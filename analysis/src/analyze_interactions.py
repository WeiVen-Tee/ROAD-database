"""
analyze_interactions.py  –  Interaction Contacts
======================================
Characterise ligand-protein contacts using pdbe-arpeggio (v1.4+).

Two analyses, each producing an aggregated per-site table of contact-type
counts/fractions and a per-atom-pair detail table:

  A  Exclusive-category ligands — appear in only one site category. Compares
     orthosteric-exclusive vs allosteric-exclusive contact profiles.
  B  Shared ligands — appear in both.

Both are computed so that every analysable site carries contact fractions in the
BUILD step's site_ligand_features.csv.

Both exclude RDKit-incompatible ligands, which requires a completed
analyze_ligand_descriptors.py run (data/ligand_descriptors.csv).

Arpeggio results are cached per site in data/arpeggio_cache/; --force bypasses
the cache. pdbe-arpeggio comes from environment.yml; --arpeggio-bin overrides
the path.

Outputs
-------
  data/interaction_contacts[_shared].csv        aggregated contact-type counts/fractions
  data/interaction_contacts_detail[_shared].csv per-atom-pair rows
  data/arpeggio_cache/                          per-site cached contact lists

Usage
-----
  # Full run (remote, CIF files in CIF_DIR):
  python analyze_interactions.py [--verbose]

  # Explicit arpeggio binary:
  python analyze_interactions.py --arpeggio-bin /path/to/pdbe-arpeggio

  # Recompute everything, bypassing the per-site Arpeggio cache:
  python analyze_interactions.py --force
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from multiprocessing import Pool, cpu_count

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import config  # must precede cli_utils/io_utils: injects collection/src onto sys.path
import cli_utils
import io_utils
from site_pipeline_utils import iter_site_residues

logger = logging.getLogger(__name__)

CATEGORIES = config.CATEGORIES

RDKIT_INCOMPAT_QUALITIES = {'failed', 'partial'}

# weak_hbond / weak_polar are merged into their base type before counting or
# plotting
_CONTACT_TYPE_MERGES = {'weak_hbond': 'hbond', 'weak_polar': 'polar'}

_SHOWN_CONTACT_TYPES = config.SHOWN_CONTACT_TYPES

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O', 'OXT'}

# Plot data: the interaction_contacts[_shared].csv files saved by this module
# serve directly as plotting input. No separate plot CSV is written.


# ─────────────────────────────────────────────────────────────────────────────
# Data loading and filtering
# ─────────────────────────────────────────────────────────────────────────────

def load_rdkit_incompat_codes(data_dir: str) -> set[str]:
    """CCD codes whose SMILES failed or only partially parsed, from ligand_descriptors.csv.

    These are excluded from both analyses. Ligand type is not a criterion.
    """
    desc_path = os.path.join(data_dir, 'ligand_descriptors.csv')
    if not os.path.exists(desc_path):
        sys.exit(
            f'ERROR: {desc_path} not found.\n'
            'Run analyze_ligand_descriptors.py first.'
        )
    desc = pd.read_csv(desc_path, low_memory=False)
    mask = desc['parse_quality'].isin(RDKIT_INCOMPAT_QUALITIES)
    codes = set(desc.loc[mask, 'ccd_code'].dropna())
    logger.info('RDKit-incompatible CCD codes: %d (from %d total)', len(codes), len(desc))
    return codes


def load_base(enriched_csv: str) -> pd.DataFrame:
    """Load the enriched sites CSV, adding 'ccd_code' and 'site_category'.

    Rows with a missing min_dist are dropped — they cannot be classified.
    """
    # dropna on min_dist first: a NaN fails the <= comparison and would be
    # silently classified allosteric. Every other stage drops these rows too.
    df = io_utils.load_enriched_csv(enriched_csv, dropna_subset=['min_dist'])
    # The CCD code is the part after the last colon in ligand_code (e.g. 'B.126:UPA' → 'UPA')
    df = io_utils.extract_ccd_code_rsplit(df)
    df = io_utils.assign_site_category_col(df, config.assign_site_category)
    return df


def load_filtered_sites(df: pd.DataFrame, incompat_codes: set[str],
                        data_dir: str) -> pd.DataFrame:
    """Analysis A's site set: exclusive-category ligands that are RDKit-compatible.

    This is the set that gets plotted and statistically tested downstream.
    """
    df = config.apply_exclusive_filter(df.copy(), data_dir, logger)
    before = len(df)
    df = df[~df['ccd_code'].isin(incompat_codes)].copy()
    logger.info('RDKit-compatible filter: %d → %d sites', before, len(df))
    return df.reset_index(drop=True)


def load_shared_sites(df: pd.DataFrame, incompat_codes: set[str]) -> pd.DataFrame:
    """Analysis B's site set: ligands appearing in both categories, RDKit-compatible.

    Shared ligands are absent from the exclusive-category CSV by definition, so
    this is the complement of load_filtered_sites within the analysable set.
    Contacts are computed for them so every analysable site carries frac_* in the
    feature table, but they are deliberately not plotted or tested.
    """
    cats_per_code = df.groupby('ccd_code')['site_category'].apply(set)
    shared_codes = set(cats_per_code[cats_per_code.apply(len) == 2].index)
    shared_codes -= incompat_codes
    df = df[df['ccd_code'].isin(shared_codes)].copy()
    logger.info('Shared-category CCD codes: %d → %d sites', len(shared_codes), len(df))
    return df.reset_index(drop=True)


# ─────────────────────────────────────────────────────────────────────────────
# CIF helpers
# ─────────────────────────────────────────────────────────────────────────────

def _cif_path(pdb_id: str, cif_dir: str) -> str:
    """Path to a PDB entry's assembly CIF, e.g. /path/to/dir/1ABC-assembly1.cif."""
    return os.path.join(cif_dir, f'{pdb_id.upper()}-assembly1.cif')


# ─────────────────────────────────────────────────────────────────────────────
# Arpeggio runner and contact parsing
# ─────────────────────────────────────────────────────────────────────────────

def _build_site_res_set(site_res_list_str: str) -> set[tuple[str, int]]:
    """Site residues as (chain, resnum) tuples, insertion code discarded.

    Arpeggio reports partners by auth_seq_id, which carries no insertion code,
    so matching has to happen on the bare number.
    """
    return {(chain, seqnum)
            for chain, seqnum, _icode in iter_site_residues(site_res_list_str)}


def _parse_ligand_code(ligand_code: str) -> tuple[str, int, str]:
    """Parse a ligand_code string into its components.

    Example: 'B.126:UPA' → chain='B', resnum=126, ccd_code='UPA'.
    A trailing insertion-code letter (if any) is stripped; the numeric part's
    sign is preserved so negative auth_seq_id values parse correctly.

    Args:
        ligand_code: String in the format 'CHAIN.RESNUM:CCDCODE'.

    Returns:
        Tuple of (chain_id, residue_number_int, ccd_code_string).
    """
    loc, ccd_code = ligand_code.rsplit(':', 1)
    chain, resnum_str = loc.split('.', 1)
    resnum_str = resnum_str.strip()
    if resnum_str and resnum_str[-1].isalpha():
        resnum_str = resnum_str[:-1]
    return chain.strip(), int(resnum_str), ccd_code.strip()


def _site_res_hash(site_res_set: set[tuple[str, int]]) -> str:
    """Short stable hash of a site's residue set (10 hex chars).

    Folded into the cache filename so that a site whose residue definition
    changed upstream misses the cache instead of reusing contacts that were
    pre-filtered against the old residue set.
    """
    key = ';'.join(f'{chain}.{resnum}' for chain, resnum in sorted(site_res_set))
    return hashlib.md5(key.encode()).hexdigest()[:10]


def _cache_path(cache_dir: str, pdb_id: str, chain: str,
                resnum: int, ccd_code: str, site_hash: str) -> str:
    """Per-site cache path, keyed on PDB/chain/resnum/CCD plus the residue-set hash."""
    key = f'{pdb_id.upper()}_{chain}_{resnum}_{ccd_code}_{site_hash}.json'
    return os.path.join(cache_dir, key)


def _run_arpeggio(cif_file: str, chain: str, resnum: int,
                  out_dir: str, arpeggio_bin: str) -> str | None:
    """Run pdbe-arpeggio on one ligand residue; returns its JSON path, or None on failure."""
    cmd = [arpeggio_bin, cif_file, '-s', f'/{chain}/{resnum}/', '-o', out_dir, '-m']
    try:
        result = subprocess.run(cmd, capture_output=True, timeout=180)
    except subprocess.TimeoutExpired:
        logger.warning('Arpeggio timed out for %s chain=%s res=%d', cif_file, chain, resnum)
        return None
    except FileNotFoundError:
        sys.exit(
            f'ERROR: arpeggio binary not found: {arpeggio_bin}\n'
            'Check that pdbe-arpeggio is installed (see environment.yml) '
            'or pass --arpeggio-bin explicitly.'
        )

    if result.returncode != 0:
        logger.warning(
            'Arpeggio failed (rc=%d) for %s chain=%s res=%d: %s',
            result.returncode, cif_file, chain, resnum,
            result.stderr.decode(errors='replace')[:200]
        )
        return None

    stem = os.path.splitext(os.path.basename(cif_file))[0]
    json_path = os.path.join(out_dir, f'{stem}.json')
    return json_path if os.path.exists(json_path) else None


def _filter_site_contacts(json_path: str, chain: str, resnum: int, ccd_code: str,
                           site_res_set: set[tuple[str, int]]) -> list[dict]:
    """Keep only INTER contacts between this ligand and a residue of this site.

    Intra-molecular contacts and contacts with non-site residues are discarded.
    The filtered list is what gets cached, and feeds both the aggregate and
    detail tables.
    """
    with open(json_path) as f:
        contacts = json.load(f)

    kept = []
    for entry in contacts:
        # Only keep contacts between two different molecules (INTER).
        # INTRA = within the same chain/molecule, which is not a ligand-protein contact.
        if entry.get('interacting_entities') != 'INTER':
            continue

        bgn = entry['bgn']
        end = entry['end']

        bgn_is_lig = (bgn['auth_asym_id'] == chain and
                      bgn['auth_seq_id'] == resnum and
                      bgn['label_comp_id'] == ccd_code)
        end_is_lig = (end['auth_asym_id'] == chain and
                      end['auth_seq_id'] == resnum and
                      end['label_comp_id'] == ccd_code)

        if bgn_is_lig:
            lig_end, partner = bgn, end
        elif end_is_lig:
            lig_end, partner = end, bgn
        else:
            continue

        if (partner['auth_asym_id'], partner['auth_seq_id']) not in site_res_set:
            continue

        kept.append({
            'lig_atom':      lig_end.get('auth_atom_id', ''),
            'partner_chain': partner['auth_asym_id'],
            'partner_resnum': partner['auth_seq_id'],
            'partner_resname': partner.get('label_comp_id', ''),
            'partner_atom':  partner.get('auth_atom_id', ''),
            'contact':       [ct.lower() for ct in entry.get('contact', [])],
            'distance':      entry.get('distance'),
            'type':          entry.get('type', ''),
        })

    return kept


def _kept_contact_tags(raw_tags: list[str]) -> list[str]:
    """Merge weak_* tags into their base type, then keep only _SHOWN_CONTACT_TYPES.

    Returned in _SHOWN_CONTACT_TYPES order; empty if nothing maps.
    """
    merged = {_CONTACT_TYPE_MERGES.get(t, t) for t in raw_tags}
    return [ct for ct in _SHOWN_CONTACT_TYPES if ct in merged]


def _aggregate_contacts(raw_contacts: list[dict]) -> dict:
    """Per-type contact counts. Contacts with no shown tag count toward nothing."""
    counts: Counter = Counter()
    for entry in raw_contacts:
        counts.update(_kept_contact_tags(entry['contact']))
    return {'n_contacts_total': sum(counts.values()), 'counts': dict(counts)}


def _extract_detail_rows(raw_contacts: list[dict], row_dict: dict,
                         ccd_code: str) -> list[dict]:
    """One detail row per atom-pair contact, with AA identity and geometry.

    Only contacts carrying a shown tag are emitted, so the detail table stays
    consistent with the aggregated counts.
    """
    base = {
        'ID':            row_dict['ID'],
        'source_pdb':    row_dict['source_pdb'],
        'ligand_code':   row_dict['ligand_code'],
        'ccd_code':      ccd_code,
        'site_category': row_dict['site_category'],
        'min_dist':      row_dict['min_dist'],
    }
    rows = []
    for entry in raw_contacts:
        kept = _kept_contact_tags(entry['contact'])
        if not kept:
            continue
        rows.append({
            **base,
            'ligand_atom':    entry['lig_atom'],
            'partner_chain':  entry['partner_chain'],
            'partner_resnum': entry['partner_resnum'],
            'partner_resname': entry['partner_resname'],
            'partner_atom':   entry['partner_atom'],
            'contact_types':  '|'.join(kept),
            'distance':       entry['distance'] if entry['distance'] is not None else float('nan'),
            'interaction_type': entry['type'],
            'is_backbone':    entry['partner_atom'] in BACKBONE_ATOMS,
        })
    return rows


# ─────────────────────────────────────────────────────────────────────────────
# Per-site worker (called by multiprocessing pool)
# ─────────────────────────────────────────────────────────────────────────────

def _process_site(args: tuple) -> dict | None:
    """Run Arpeggio for one binding site, or load it from the per-site cache.

    Takes its inputs as one tuple (row_dict, cif_dir, arpeggio_bin, cache_dir,
    force) so it can be mapped over a multiprocessing Pool.

    Returns {'summary': record, 'details': [per-contact rows]}, or None if the
    site could not be processed.
    """
    (row_dict, cif_dir, arpeggio_bin, cache_dir, force) = args

    pdb_id      = row_dict['source_pdb']
    ligand_code = row_dict['ligand_code']
    site_res    = row_dict['binding_site']

    try:
        chain, resnum, ccd_code = _parse_ligand_code(ligand_code)
    except Exception as e:
        logger.warning('%s %s: cannot parse ligand_code: %s', pdb_id, ligand_code, e)
        return None

    site_res_set = _build_site_res_set(site_res)
    site_hash    = _site_res_hash(site_res_set)

    cached = _cache_path(cache_dir, pdb_id, chain, resnum, ccd_code, site_hash)
    if not force:
        try:
            with open(cached) as f:
                data = json.load(f)
            logger.debug('Cache hit: %s', cached)
            return _build_result(row_dict, ccd_code, data)
        except Exception:
            pass

    cif_file = _cif_path(pdb_id, cif_dir)
    if not os.path.exists(cif_file):
        logger.warning('%s: CIF not found — skipped', pdb_id)
        return None

    tmp = tempfile.mkdtemp(prefix='arpeggio_')
    try:
        json_path = _run_arpeggio(cif_file, chain, resnum, tmp, arpeggio_bin)
        if json_path is None:
            return None
        raw_contacts = _filter_site_contacts(json_path, chain, resnum,
                                             ccd_code, site_res_set)
    except Exception as exc:
        logger.warning('%s %s: unexpected error: %s', pdb_id, ligand_code, exc)
        return None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    try:
        with open(cached, 'w') as f:
            json.dump(raw_contacts, f)
    except Exception as exc:
        logger.debug('Cache write failed for %s: %s', cached, exc)

    return _build_result(row_dict, ccd_code, raw_contacts)


def _build_result(row_dict: dict, ccd_code: str,
                  raw_contacts: list[dict]) -> dict:
    """Assemble {summary, details} for one site from its filtered contact list."""
    contact_data = _aggregate_contacts(raw_contacts)
    return {
        'summary': _build_record(row_dict, ccd_code, contact_data),
        'details': _extract_detail_rows(raw_contacts, row_dict, ccd_code),
    }


def _build_record(row_dict: dict, ccd_code: str, contact_data: dict) -> dict:
    """One row of interaction_contacts.csv: metadata plus n_/frac_ per contact type.

    Sites with zero detected contacts get NaN rather than 0 throughout, so they
    drop out of the statistics instead of counting as genuine zeros.
    """
    n_total = contact_data['n_contacts_total']
    counts  = contact_data['counts']

    record = {
        'ID':               row_dict['ID'],
        'source_pdb':       row_dict['source_pdb'],
        'ligand_code':      row_dict['ligand_code'],
        'ccd_code':         ccd_code,
        'site_category':    row_dict['site_category'],
        'min_dist':         row_dict['min_dist'],
        'n_contacts_total': n_total if n_total > 0 else float('nan'),
    }

    for ct in _SHOWN_CONTACT_TYPES:
        if n_total > 0:
            record[f'n_{ct}']    = counts.get(ct, 0)
            record[f'frac_{ct}'] = counts.get(ct, 0) / n_total
        else:
            record[f'n_{ct}']    = float('nan')
            record[f'frac_{ct}'] = float('nan')

    return record


# ─────────────────────────────────────────────────────────────────────────────
# Batch runner
# ─────────────────────────────────────────────────────────────────────────────

def run_arpeggio_batch(sites_df: pd.DataFrame, cif_dir: str,
                       arpeggio_bin: str, workers: int, cache_dir: str,
                       force: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run Arpeggio over every row of sites_df, in a process pool (serial if workers=1).

    Returns (contacts_df, details_df), both empty if no site produced results.
    """
    needed_cols = ['ID', 'source_pdb', 'ligand_code', 'binding_site',
                   'site_category', 'min_dist', 'ccd_code']
    avail = [c for c in needed_cols if c in sites_df.columns]
    worker_args = [
        (row, cif_dir, arpeggio_bin, cache_dir, force)
        for row in sites_df[avail].to_dict('records')
    ]

    print(f'  Running Arpeggio on {len(worker_args)} sites ({workers} workers) …')
    t0 = time.time()

    if workers > 1:
        with Pool(processes=workers) as pool:
            results = pool.map(_process_site, worker_args)
    else:
        results = [_process_site(a) for a in worker_args]

    print(f'  Arpeggio finished in {time.time() - t0:.1f}s')

    ok = [r for r in results if r is not None]
    logger.info('Arpeggio: %d succeeded, %d failed/skipped', len(ok), len(results) - len(ok))

    if not ok:
        return pd.DataFrame(), pd.DataFrame()

    contacts_df = pd.DataFrame([r['summary'] for r in ok])
    details_df  = pd.DataFrame([row for r in ok for row in r['details']])

    n_det    = contacts_df['n_contacts_total'].notna().sum()
    frac_det = n_det / len(contacts_df)
    print(f'  Contact detection rate: {n_det}/{len(contacts_df)} ({frac_det:.0%})')
    if frac_det < 0.5:
        logger.warning('Low contact detection rate (%.0f%%)', frac_det * 100)

    return contacts_df, details_df


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for analyze_interactions.py.

    Returns:
        Namespace object with attributes: arpeggio_bin, workers, verbose, force.
    """
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--arpeggio-bin', default='pdbe-arpeggio',
                   help='Path to pdbe-arpeggio executable (default: pdbe-arpeggio)')
    p.add_argument('--workers', type=int, default=max(1, cpu_count() - 1),
                   help='Parallel workers (default: nCPU-1)')
    cli_utils.add_common_args(p, force=True)
    return p.parse_args()


def _run_analysis(label: str, sites_df: pd.DataFrame,
                  cif_dir: str,
                  arpeggio_bin: str, workers: int, cache_dir: str,
                  data_dir: str,
                  suffix: str,
                  force: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run one analysis set end to end and write its two CSVs.

    `suffix` is '' for Analysis A and '_shared' for Analysis B. Statistics live
    downstream, not here.
    """
    contacts_df, details_df = run_arpeggio_batch(
        sites_df, cif_dir, arpeggio_bin, workers, cache_dir, force
    )
    if contacts_df.empty:
        logger.warning('%s: no Arpeggio results', label)
        return pd.DataFrame(), pd.DataFrame()

    def _save(df: pd.DataFrame, name: str) -> None:
        """Save df to data_dir/name and print the path."""
        path = os.path.join(data_dir, name)
        io_utils.save_csv(df, path, logger)

    _save(contacts_df, f'interaction_contacts{suffix}.csv')

    if not details_df.empty:
        _save(details_df, f'interaction_contacts_detail{suffix}.csv')

    return contacts_df, details_df


def main() -> None:
    """Entry point for analyze_interactions.py.

    Loads the enriched sites CSV, runs Analysis A (exclusive-category ligands)
    and Analysis B (shared ligands), and saves all output CSVs. Both sets are
    computed so every analysable site gets contact fractions. 
    Human-readable summaries are produced downstream.
    """
    args = parse_args()
    t0 = cli_utils.init_run('analyze_interactions', args, config.LOG_DIR)

    data_dir  = config.DATA_DIR
    cache_dir = os.path.join(data_dir, 'arpeggio_cache')
    os.makedirs(cache_dir, exist_ok=True)

    enriched_csv = config.ENRICHED_CSV

    print('Loading data …')
    incompat_codes = load_rdkit_incompat_codes(data_dir)
    base_df        = load_base(enriched_csv)

    # ── Analysis A: exclusive-category ligands (plotted and tested) ───────────
    print('\n[Analysis A] Exclusive-category ligands …')
    sites_a = load_filtered_sites(base_df, incompat_codes, data_dir)
    n_orth  = (sites_a['site_category'] == 'orthosteric').sum()
    n_allo  = (sites_a['site_category'] == 'allosteric').sum()
    print(f'  {len(sites_a)} sites: {n_orth} orthosteric-exclusive, {n_allo} allosteric-exclusive')
    logger.info('[A] Sites: %d (ortho=%d, allo=%d)', len(sites_a), n_orth, n_allo)

    if len(sites_a) == 0:
        sys.exit('No sites pass the exclusive + RDKit-compatible filter.')

    contacts_a, _details_a = _run_analysis(
        'exclusive-category ligands', sites_a,
        config.CIF_DIR, args.arpeggio_bin,
        args.workers, cache_dir, data_dir, suffix='', force=args.force
    )
    if contacts_a.empty:
        sys.exit('[A] No Arpeggio results. Check CIF files and arpeggio binary.')

    # ── Analysis B: shared ligands (feature-table coverage only) ──────────────
    print('\n[Analysis B] Shared ligands (found in both categories) …')
    sites_b  = load_shared_sites(base_df, incompat_codes)
    n_orth_b = (sites_b['site_category'] == 'orthosteric').sum()
    n_allo_b = (sites_b['site_category'] == 'allosteric').sum()
    print(f'  {len(sites_b)} sites: {n_orth_b} orthosteric, {n_allo_b} allosteric')
    logger.info('[B] Sites: %d (ortho=%d, allo=%d)', len(sites_b), n_orth_b, n_allo_b)

    if len(sites_b) > 0:
        _run_analysis(
            'shared ligands', sites_b,
            config.CIF_DIR, args.arpeggio_bin,
            args.workers, cache_dir, data_dir, suffix='_shared',
            force=args.force
        )
    else:
        logger.warning('[B] No shared-ligand sites found.')

    cli_utils.finish_run(t0)


if __name__ == '__main__':
    main()
