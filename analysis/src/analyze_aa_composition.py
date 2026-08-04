"""
analyze_aa_composition.py
==========================
Per-site amino-acid composition of binding sites: physicochemical property
group fractions, plus a mean per-site AA fraction summary for the heatmap.

Expands the enriched CSV to residue level using collection's pdb_cif_cache,
which keys residues on (seqnum, icode).

Inputs
------
  <pdb_cif_cache.p>                   (structure lookup cache from collection)
  <nonredundant_sites_annotated_enriched.csv>

Outputs
-------
  data/site_aa_composition.csv        (one row per site: site_size + frac_* + property fractions)
  data/aa_mean_fractions.csv          (aa, site_category, n_sites, mean_fraction — heatmap input)

Usage
-----
    python analyze_aa_composition.py [--verbose]
"""

import argparse
import logging
import os
import sys
from collections import Counter

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import config  # must precede cli_utils/io_utils: injects collection/src onto sys.path
import cli_utils
import io_utils
from site_pipeline_utils import iter_site_residues
from utils import load_pickle

logger = logging.getLogger(__name__)

CATEGORIES  = config.CATEGORIES
STANDARD_AA = config.STANDARD_AA


# ─────────────────────────────────────────────────────────────────────────────
# Residue-level lookup
# ─────────────────────────────────────────────────────────────────────────────

def _parse_res_list(raw: str) -> tuple[list[tuple[str, int, str]], int]:
    """Site residues as (chain, resnum, icode) triples, plus a count of bad tokens.

    The count covers only tokens that could not be tokenised at all (no '.', or
    a non-numeric residue number).
    """
    results  = list(iter_site_residues(raw))
    n_listed = sum(1 for t in str(raw).split(';') if t.strip())
    return results, n_listed - len(results)


def _lookup_aa(seqid_to_aa: dict, chain: str, resnum: int, icode: str,
               pdb_id: str) -> str | None:
    """One-letter AA code from a pdb_cif_cache entry.

    seqid_to_aa keys are 'MODEL.CHAIN'; falls back from model 1 to any other
    model for structures with no model-1 record. Residues within a chain are
    keyed by (seqnum, icode) — residue number alone is ambiguous, since 100,
    100A and 100B all share seqid.num == 100.
    """
    chain_key = f'1.{chain}'
    chain_map = seqid_to_aa.get(chain_key)

    if chain_map is None:
        for key, cmap in seqid_to_aa.items():
            if key.endswith(f'.{chain}'):
                chain_map = cmap
                break

    if chain_map is None:
        logger.warning('%s: chain key %r not found in seqid_to_aa', pdb_id, chain_key)
        return None

    aa = chain_map.get((resnum, icode))
    if aa is None:
        logger.warning('%s: resnum %d%s not found in chain %r',
                       pdb_id, resnum, icode, chain_key)
    return aa


def build_residue_table(df: pd.DataFrame, cache: dict) -> pd.DataFrame:
    """Expand the binding-site table to one row per residue, with aa_code."""
    records = []
    n_sites       = len(df)
    n_pdb_missing = 0
    n_res_total   = 0
    n_lookup_fail = 0
    n_x_excluded  = 0
    n_unparseable = 0

    for _, row in df.iterrows():
        pdb_id = str(row['source_pdb']).upper()
        entry  = cache.get(pdb_id)

        if entry is None or entry.get('seqid_to_aa') is None:
            logger.warning('PDB %s not in cache or seqid_to_aa is None — skipping site', pdb_id)
            n_pdb_missing += 1
            continue

        seqid_to_aa  = entry['seqid_to_aa']
        site_cat     = row.get('site_category') or config.assign_site_category(row['min_dist'])
        res_list_raw = str(row['binding_site'])

        residues, n_dropped = _parse_res_list(res_list_raw)
        n_unparseable += n_dropped
        if not residues:
            logger.warning('PDB %s ligand %s: empty site_res_list', pdb_id, row['ligand_code'])
            continue

        site_residues = []
        for chain, resnum, icode in residues:
            n_res_total += 1
            aa = _lookup_aa(seqid_to_aa, chain, resnum, icode, pdb_id)
            if aa is None:
                n_lookup_fail += 1
                continue
            if aa == 'X':
                n_x_excluded += 1
                continue
            site_residues.append(aa)

        site_size = len(site_residues)
        for aa in site_residues:
            records.append({
                'ID':            row['ID'],
                'source_pdb':    pdb_id,
                'ligand_code':   row['ligand_code'],
                'site_category': site_cat,
                'site_size':     site_size,
                'aa_code':       aa,
            })

    n_listed = n_res_total + n_unparseable
    print(f'  Sites processed      : {n_sites:,}')
    print(f'  PDB entries missing  : {n_pdb_missing:,}')
    print(f'  Residues attempted   : {n_res_total:,}')
    print(f'  Unparseable tokens   : {n_unparseable:,}  '
          f'({100*n_unparseable/max(n_listed,1):.2f}% of listed)')
    print(f'  Lookup failures      : {n_lookup_fail:,}  '
          f'({100*n_lookup_fail/max(n_res_total,1):.2f}%)')
    print(f'  Modified (X) excluded: {n_x_excluded:,}  '
          f'({100*n_x_excluded/max(n_res_total,1):.2f}%)')

    logger.info('Sites: %d; PDB missing: %d; unparseable tokens: %d; '
                'lookup fail: %d; X excluded: %d',
                n_sites, n_pdb_missing, n_unparseable, n_lookup_fail, n_x_excluded)
    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────────────────────
# Per-site property fractions (consumed directly by run_pipeline.py's BUILD step)
# ─────────────────────────────────────────────────────────────────────────────

def compute_per_site_properties(df: pd.DataFrame) -> pd.DataFrame:
    """Per-site AA fractions: 20 frac_* columns plus the 4 property-group columns.

    One row per site, with site_size as the denominator used for both.
    """
    records = []
    for (sid, pdb, lig, cat), grp in df.groupby(
            ['ID', 'source_pdb', 'ligand_code', 'site_category']):
        counts = Counter(grp['aa_code'])
        total  = len(grp)

        aa_fracs   = {f'frac_{aa}': counts.get(aa, 0) / total for aa in STANDARD_AA}
        prop_fracs = {prop: sum(counts.get(a, 0) for a in aas) / total
                      for prop, aas in config.AA_PROPERTY_GROUPS.items()}

        rec = {'ID': sid, 'source_pdb': pdb, 'ligand_code': lig, 'site_category': cat,
               'site_size': total}
        rec.update(aa_fracs)
        rec.update(prop_fracs)
        records.append(rec)

    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────────────────────
# Mean per-site AA fraction (heatmap input — derived from the table above)
# ─────────────────────────────────────────────────────────────────────────────

def compute_mean_fractions(site_df: pd.DataFrame) -> pd.DataFrame:
    """Mean per-site AA fraction per category — the heatmap's input.

    Reuses compute_per_site_properties' frac_* columns rather than re-counting.
    Columns: aa, site_category, n_sites, mean_fraction.
    """
    rows = []
    for aa in STANDARD_AA:
        col = f'frac_{aa}'
        for cat in CATEGORIES:
            x = site_df.loc[site_df['site_category'] == cat, col].values
            rows.append({
                'aa':            aa,
                'site_category': cat,
                'n_sites':       len(x),
                'mean_fraction': float(x.mean()) if len(x) else 0.0,
            })
    return pd.DataFrame(rows)


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    """Expand sites to residues, then write per-site and mean-per-category AA fractions."""
    parser = argparse.ArgumentParser(
        description='Per-site amino-acid composition: heatmap summary.')
    cli_utils.add_common_args(parser)
    args = parser.parse_args()

    t0 = cli_utils.init_run('analyze_aa_composition', args, config.LOG_DIR)

    # ── [1/3] Load inputs ──────────────────────────────────────────────────────
    print('\n[1/3] Loading pdb_cif_cache and enriched CSV …')
    logger.info('Loading pdb_cif_cache from %s', config.PDB_CIF_CACHE)
    cache = load_pickle(config.PDB_CIF_CACHE)
    print(f'  Cache: {len(cache):,} PDB entries.')

    df = io_utils.load_enriched_csv(config.ENRICHED_CSV)
    print(f'  Enriched CSV: {len(df):,} rows.')
    before = len(df)
    df = df.dropna(subset=['min_dist', 'binding_site'])
    dropped = before - len(df)
    if dropped:
        print(f'  Dropped {dropped} rows with missing min_dist or site_res_list.')
        logger.warning('Dropped %d rows with missing min_dist/site_res_list', dropped)

    # ── [2/3] Expand to residue level + compute per-site fractions ─────────────
    print('\n[2/3] Expanding sites to residue level + computing per-site AA fractions …')
    residue_df  = build_residue_table(df, cache)
    per_site_all = compute_per_site_properties(residue_df)

    # Tag each site as exclusive-category (or not) — consumed by other stages
    # (e.g. descriptors, interactions); not used to filter this stage's own output.
    excl_path = os.path.join(config.DATA_DIR, 'ligand_exclusive_categories.csv')
    if os.path.exists(excl_path):
        excl_ccds = set(pd.read_csv(excl_path)['ccd_code'])
        per_site_all = io_utils.extract_ccd_code_rsplit(per_site_all)
        per_site_all['is_exclusive'] = per_site_all['ccd_code'].isin(excl_ccds)
        per_site_all = per_site_all.drop(columns=['ccd_code'])
    else:
        logger.warning('ligand_exclusive_categories.csv not found — marking all sites as non-exclusive')
        per_site_all['is_exclusive'] = False

    aa_comp_path = os.path.join(config.DATA_DIR, 'site_aa_composition.csv')
    io_utils.save_csv(per_site_all, aa_comp_path, logger)

    n_sites = len(per_site_all)
    print(f'  {n_sites:,} binding sites.')
    for cat in CATEGORIES:
        n = (per_site_all['site_category'] == cat).sum()
        print(f'    {cat:<14}: {n:,} sites')

    # ── [3/3] Mean per-site AA fractions for the heatmap ────────────────────────
    print('\n[3/3] Deriving mean per-site AA fractions for the heatmap …')
    mean_df = compute_mean_fractions(per_site_all)
    enrichment_path = os.path.join(config.DATA_DIR, 'aa_mean_fractions.csv')
    io_utils.save_csv(mean_df, enrichment_path, logger)

    cli_utils.finish_run(t0)


if __name__ == '__main__':
    main()
