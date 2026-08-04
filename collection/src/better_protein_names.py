"""
Resolve per-site, chain-specific protein names from CIF entity records.

More informative than the UniProt-level protein_name, which for a polyprotein
lists every cleaved product in one string. For each site, the entity
description (`_entity.pdbx_description`) of its chains is looked up via
`_entity_poly.pdbx_strand_id`, whose author chain IDs already match
site_res_list including NCS forms like 'A-2'. Old-style 'PROTEIN (Actual Name)'
descriptions are unwrapped.

Runs between nonredundant_build_* and nonredundant_enrich_*, reading the
pre-enrichment site tables. Site keys come from the union of the annotated and
complete tables, since the two modes cluster independently and neither is a
superset of the other.

Inputs:  _nonredundant_sites_{annotated,complete}.csv, {PDB_ID}-assembly1.cif
Output:  result/protein_names/protein_names_pdb.csv — one row per site:
         source_pdb, ligand_code, protein_name_pdb, chains_in_site.

Usage:
    python better_protein_names.py [--cif-dir DIR] [-v] [-f]
"""

import argparse
import os
import re
import sys
from collections import Counter

import gemmi
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    CIF_DIR, NONREDUNDANT_CSV_ANNOTATED, NONREDUNDANT_CSV_COMPLETE,
    RESULT_PROTEIN_NAMES_DIR, LOG_DIR,
)
from utils import setup_logging

log = None  # set in main()


# ─────────────────────────────────────────────────────────────────────────────
# Site-list / CIF-resolution helpers
# ─────────────────────────────────────────────────────────────────────────────

def _parse_site_res_list(s: str) -> list[tuple[str, str]]:
    """Parse a semicolon-separated residue list string into (chain, resnum) pairs.

    The input format used in the enriched CSV is 'A.7;A.11;B.100A;...'
    where each token is 'CHAIN.RESNUM' and RESNUM may include an insertion
    code suffix (e.g. '100A').
    """
    result = []
    for token in str(s).split(';'):
        token = token.strip()
        if not token:
            continue
        chain, resnum = token.split('.', 1)
        result.append((chain.strip(), resnum.strip()))
    return result


def _build_pdb_sites(enrich_df: pd.DataFrame, all_pdbs: list[str]) -> dict[str, list[dict]]:
    """Build a per-structure index of binding sites from the enriched DataFrame.

    Each site dict has keys: ligand_code, residues (list of (chain, resnum_str) tuples).
    """
    pdb_set = set(all_pdbs)
    pdb_sites: dict[str, list[dict]] = {p: [] for p in all_pdbs}
    for _, row in enrich_df.iterrows():
        pdb_id = str(row['source_pdb']).upper()
        if pdb_id not in pdb_set:
            continue
        pdb_sites[pdb_id].append({
            'ligand_code': row['ligand_code'],
            'residues':    _parse_site_res_list(row['site_res_list']),
        })
    return pdb_sites


def _resolve_cif_path(pdb_id: str, cif_dir: str) -> str | None:
    """Path to a structure's assembly CIF, or None if it is not in cif_dir."""
    cif_file = os.path.join(cif_dir, f'{pdb_id.upper()}-assembly1.cif')
    if not os.path.exists(cif_file):
        log.warning('%s: CIF not found', pdb_id)
        return None
    return cif_file


# ─────────────────────────────────────────────────────────────────────────────
# Entity-name extraction
# ─────────────────────────────────────────────────────────────────────────────

_PROTEIN_WRAPPER_RE = re.compile(r'^PROTEIN \((.*)\)$')


def _strip_protein_wrapper(desc: str) -> str:
    """Old-style PDB entries describe entities as 'PROTEIN (Actual Name)'."""
    m = _PROTEIN_WRAPPER_RE.match(desc)
    return m.group(1) if m else desc


def _entity_names_for_structure(cif_path: str) -> dict[str, str]:
    """Author chain ID -> entity description, from _entity and _entity_poly.

    pdbx_strand_id already carries author chain IDs including NCS-expanded
    forms ('A-2'), so no renaming is needed. Empty dict if either category is
    missing or malformed.
    """
    doc = gemmi.cif.read(cif_path)
    block = doc.sole_block()
    ent_cat = block.find_mmcif_category('_entity.')
    poly_cat = block.find_mmcif_category('_entity_poly.')
    if not ent_cat or not poly_cat:
        return {}

    desc_by_entity = {
        eid: _strip_protein_wrapper(gemmi.cif.as_string(desc))
        for eid, desc in zip(ent_cat.find_column('id'), ent_cat.find_column('pdbx_description'))
    }

    chain_to_desc: dict[str, str] = {}
    for entity_id, strand_ids in zip(poly_cat.find_column('entity_id'),
                                     poly_cat.find_column('pdbx_strand_id')):
        desc = desc_by_entity.get(entity_id)
        if not desc or desc == '?':
            continue
        for chain in gemmi.cif.as_string(strand_ids).split(','):
            chain = chain.strip()
            if chain:
                chain_to_desc[chain] = desc
    return chain_to_desc


def compute_names_for_structure(
        pdb_id: str, sites: list[dict], cif_dir: str) -> list[dict]:
    """protein_name_pdb for every site in one structure.

    Returns one dict per resolved site with keys source_pdb, ligand_code,
    protein_name_pdb, chains_in_site; empty if the structure cannot be loaded.
    """
    cif_file = _resolve_cif_path(pdb_id, cif_dir)
    if cif_file is None:
        return []

    try:
        chain_to_desc = _entity_names_for_structure(cif_file)
    except Exception as exc:
        log.warning('%s: CIF parse failed: %s', pdb_id, exc)
        return []

    if not chain_to_desc:
        log.warning('%s: no entity/entity_poly data — skipped', pdb_id)
        return []

    records = []
    for site in sites:
        chain_counts = Counter(chain for chain, _ in site['residues'])
        desc_counts: Counter = Counter()
        for chain, n in chain_counts.items():
            if chain in chain_to_desc:
                desc_counts[chain_to_desc[chain]] += n
        descs_by_freq = [d for d, _ in desc_counts.most_common()]
        if not descs_by_freq:
            continue
        records.append({
            'source_pdb':       pdb_id,
            'ligand_code':      site['ligand_code'],
            'protein_name_pdb': ' / '.join(descs_by_freq),
            'chains_in_site':   ';'.join(sorted(chain_counts)),
        })
    return records


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main(cif_dir: str = CIF_DIR, verbose: bool = False,
         force: bool = False) -> None:
    """Resolve per-site, chain-specific protein names for all binding sites.

    Reads the pre-enrichment nonredundant site tables (union of the annotated
    and complete modes), groups sites by structure, and for each structure
    extracts entity descriptions per chain from its CIF file. Results are
    cached in result/protein_names/protein_names_pdb.csv; use force=True to
    recompute.
    """
    global log
    os.makedirs(RESULT_PROTEIN_NAMES_DIR, exist_ok=True)
    log = setup_logging(verbose, os.path.join(LOG_DIR, 'better_protein_names.log'))

    out_path = os.path.join(RESULT_PROTEIN_NAMES_DIR, 'protein_names_pdb.csv')
    if os.path.exists(out_path) and not force:
        print(f'  {out_path} exists — skipping (--force to recompute).')
        return

    print('\nLoading pre-enrichment nonredundant site tables …')
    complete_df  = pd.read_csv(NONREDUNDANT_CSV_COMPLETE)
    annotated_df = pd.read_csv(NONREDUNDANT_CSV_ANNOTATED)
    # annotated/complete cluster independently and can pick different cluster
    # representatives for the same (uniprot_ac, ligand) group, so neither table
    # alone is a superset of the other — take the union, deduped by site key
    # (intersecting keys carry byte-identical site_res_list either way).
    enrich_df = (pd.concat([complete_df, annotated_df], ignore_index=True)
                 .dropna(subset=['site_res_list'])
                 .drop_duplicates(subset=['source_pdb', 'ligand_code']))
    print(f'  {len(enrich_df):,} binding sites, '
          f'{enrich_df["source_pdb"].nunique():,} structures.')

    all_pdbs = sorted(enrich_df['source_pdb'].str.upper().unique())
    pdb_sites = _build_pdb_sites(enrich_df, all_pdbs)

    print(f'\nResolving protein names for {len(all_pdbs):,} structures …')
    all_records = []
    n_ok = n_fail = 0
    for i, pdb_id in enumerate(all_pdbs):
        if (i + 1) % 50 == 0 or i == 0:
            print(f'  [{i+1}/{len(all_pdbs)}] {pdb_id}', flush=True)
        records = compute_names_for_structure(pdb_id, pdb_sites[pdb_id], cif_dir)
        if records:
            all_records.extend(records); n_ok += 1
        else:
            n_fail += 1

    print(f'  OK: {n_ok:,}  failed: {n_fail:,}  site rows: {len(all_records):,}')
    if all_records:
        names_df = pd.DataFrame(all_records)
        names_df.to_csv(out_path, index=False)
        log.info('Wrote %d rows to %s', len(names_df), out_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Resolve per-site, chain-specific protein names from CIF entity records.')
    parser.add_argument('--cif-dir', default=CIF_DIR,
                        help=f'CIF directory (default: {CIF_DIR})')
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Recompute even if the output CSV already exists')
    args = parser.parse_args()
    main(args.cif_dir, args.verbose, args.force)
