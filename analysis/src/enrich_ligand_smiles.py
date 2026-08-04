"""
enrich_ligand_smiles.py
======================================
Resolve canonical SMILES for every unique ligand (CCD code) in
nonredundant_sites_annotated_enriched.csv, using the local PDB Chemical
Component Dictionary CIF file as the primary source.

Source priority:
  1. OpenEye SMILES_CANONICAL from ccd/Components-rel-alt.cif
     (explicit-bond notation, parses cleanly in RDKit)
  2. InChI from the same CIF block (roundtrip via MolFromInchi)
  3. RCSB REST API fallback (for any CCD codes absent from the CIF)

smiles_source values written to output CSV:
  'cif_smiles'    parsed directly from CIF OpenEye SMILES
  'cif_inchi'     resolved via InChI roundtrip from CIF
  'api_smiles'    parsed from RCSB API SMILES
  'api_inchi'     resolved via InChI roundtrip from RCSB API
  'partial'       partial RDKit sanitization (kekulize skipped)
  'rdkit_failed'  all RDKit attempts failed (canonical_smiles = NaN)
  'api_failed'    CIF miss + API returned nothing

Output
------
data/ligand_smiles.csv
    ccd_code, smiles_source, canonical_smiles

Usage
-----
    python enrich_ligand_smiles.py [--verbose] [--force]
"""

import argparse
import logging
import os
import sys
import time
import warnings

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import pandas as pd

import config  # must precede utils/cli_utils/io_utils: injects collection/src onto sys.path
from utils import create_session
import cli_utils
import io_utils

logger = logging.getLogger(__name__)

RCSB_CHEMCOMP_URL = 'https://data.rcsb.org/rest/v1/core/chemcomp/{ccd}'
API_SLEEP  = 0.05   # seconds between requests
API_TIMEOUT = 10    # seconds per request

# Sanity checks: (ccd_code, description, predicate(mol) -> bool)
SANITY_CHECKS = [
    ('ATP', 'ExactMolWt ≈ 507 Da',   lambda m: abs(_mw(m) - 507.18) < 1.0),
    ('ADP', 'ExactMolWt ≈ 427 Da',   lambda m: abs(_mw(m) - 427.20) < 1.0),
    ('HEM', '≥4 rings (porphyrin)',   lambda m: _n_rings(m) >= 4),
]


def _mw(mol):
    from rdkit.Chem.Descriptors import ExactMolWt
    return ExactMolWt(mol)


def _n_rings(mol):
    from rdkit.Chem import rdMolDescriptors
    return rdMolDescriptors.CalcNumRings(mol)


# ─────────────────────────────────────────────────────────────────────────────
# CIF parsing
# ─────────────────────────────────────────────────────────────────────────────

def build_cif_lookup(cif_path: str) -> dict:
    """
    Parse the full CCD CIF file once and return a lookup dict:

        ccd_code -> {
            'smiles':    str | None,   # OpenEye SMILES_CANONICAL (preferred)
            'inchi':     str | None,   # InChI string
            'chem_type': str | None,   # _chem_comp.type (e.g. 'NON-POLYMER')
            'atoms':     list of (atom_id, element, x_ideal, y_ideal, z_ideal)
                         element is title-case; coordinates may be None if '?'
            'bonds':     list of (atom_id_1, atom_id_2, value_order)
        }
    """
    import gemmi
    print(f'  Parsing CIF: {cif_path} …', flush=True)
    t0 = time.time()
    doc = gemmi.cif.read(cif_path)
    print(f'  CIF parsed in {time.time()-t0:.1f}s  ({len(doc)} blocks)', flush=True)

    lookup = {}
    for block in doc:
        ccd = block.name

        # ── SMILES / InChI ────────────────────────────────────────────────
        smiles = inchi = None
        try:
            desc_table = block.find(['_pdbx_chem_comp_descriptor.type',
                                     '_pdbx_chem_comp_descriptor.program',
                                     '_pdbx_chem_comp_descriptor.descriptor'])
            for row in desc_table:
                dtype = row[0].upper().strip('"')
                prog  = row[1].upper()
                val   = row[2].strip('"').strip()
                if not val or val == '?':
                    continue
                if dtype == 'SMILES_CANONICAL':
                    # Always prefer OpenEye SMILES; accept others only as fallback
                    if 'OPENEYE' in prog:
                        smiles = val
                    elif smiles is None:
                        smiles = val
                elif dtype == 'INCHI' and inchi is None:
                    inchi = val
        except Exception:
            pass

        # ── _chem_comp.type ────────────────────────────────────────────────
        chem_type = None
        try:
            chem_type = block.find_value('_chem_comp.type')
            if chem_type:
                chem_type = chem_type.strip('"').strip()
        except Exception:
            pass

        # ── Atom table (for 3D coords) ─────────────────────────────────────
        atoms = []
        try:
            atom_table = block.find(['_chem_comp_atom.atom_id',
                                     '_chem_comp_atom.type_symbol',
                                     '_chem_comp_atom.pdbx_model_Cartn_x_ideal',
                                     '_chem_comp_atom.pdbx_model_Cartn_y_ideal',
                                     '_chem_comp_atom.pdbx_model_Cartn_z_ideal'])
            for row in atom_table:
                element = row[1].strip('"').strip()
                if element.upper() == 'H':
                    continue
                # CIF uses all-caps (e.g. 'MG', 'AS'); RDKit needs title-case
                element = element[0].upper() + element[1:].lower()
                x = None if row[2].strip() == '?' else row[2].strip()
                y = None if row[3].strip() == '?' else row[3].strip()
                z = None if row[4].strip() == '?' else row[4].strip()
                atoms.append((row[0].strip('"'), element, x, y, z))
        except Exception:
            pass

        # ── Bond table ─────────────────────────────────────────────────────
        bonds = []
        try:
            bond_table = block.find(['_chem_comp_bond.atom_id_1',
                                     '_chem_comp_bond.atom_id_2',
                                     '_chem_comp_bond.value_order'])
            for row in bond_table:
                bonds.append((row[0].strip('"'), row[1].strip('"'),
                               row[2].strip('"').upper()))
        except Exception:
            pass

        if smiles or inchi or atoms:
            lookup[ccd] = {
                'smiles': smiles, 'inchi': inchi,
                'chem_type': chem_type,
                'atoms': atoms, 'bonds': bonds,
            }

    logger.info('CIF lookup built: %d entries in %.1fs', len(lookup), time.time()-t0)
    return lookup


# ─────────────────────────────────────────────────────────────────────────────
# SMILES sanitisation
# ─────────────────────────────────────────────────────────────────────────────

def sanitise_smiles(raw_smiles: str | None,
                    inchi: str | None = None) -> tuple[str | None, str]:
    """
    Parse and re-canonicalise a SMILES with RDKit.

    Attempt order:
      1. Standard MolFromSmiles (strict)
      2. InChI roundtrip via MolFromInchi — resolves aromatic-notation
         ambiguities that cause kekulization failures
      3. MolFromSmiles(sanitize=False) + partial sanitization (skip kekulize)

    Returns (canonical_smiles | None, source_tag).
    source_tag: 'smiles' | 'inchi' | 'partial' | 'failed'
    """
    from rdkit import Chem

    if raw_smiles:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            mol = Chem.MolFromSmiles(raw_smiles)
        if mol is not None:
            return Chem.MolToSmiles(mol, canonical=True), 'smiles'

    if inchi:
        try:
            from rdkit.Chem.inchi import MolFromInchi
            mol = MolFromInchi(inchi, sanitize=True, removeHs=True)
            if mol is not None:
                return Chem.MolToSmiles(mol, canonical=True), 'inchi'
        except Exception:
            pass

    if raw_smiles:
        try:
            mol = Chem.MolFromSmiles(raw_smiles, sanitize=False)
            if mol is not None:
                Chem.SanitizeMol(mol,
                                 Chem.SanitizeFlags.SANITIZE_ALL ^
                                 Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                # MolToSmiles still kekulizes internally, so this branch may
                # still fail for un-kekulizable molecules.  allBondsExplicit=False
                # keeps aromatic SMILES notation, which is less likely to produce
                # an invalid string than the explicit Kekulé form.
                smi = Chem.MolToSmiles(mol, canonical=True, allBondsExplicit=False)
                if smi:
                    return smi, 'partial'
        except Exception:
            pass

    return None, 'failed'


# ─────────────────────────────────────────────────────────────────────────────
# API fallback
# ─────────────────────────────────────────────────────────────────────────────

def fetch_smiles_api(ccd_code: str, session) -> tuple[str | None, str | None]:
    """Fetch SMILES and InChI from RCSB REST API. Returns (smiles, inchi)."""
    url = RCSB_CHEMCOMP_URL.format(ccd=ccd_code)
    try:
        resp = session.get(url, timeout=API_TIMEOUT)
        if resp.status_code == 404:
            return None, None
        resp.raise_for_status()
        data = resp.json()
        smiles = inchi = None
        for source_key in ('rcsb_chem_comp_descriptor', 'pdbx_chem_comp_descriptor'):
            descs = data.get(source_key, [])
            if not isinstance(descs, list):
                descs = [descs]
            for desc in descs:
                dtype = desc.get('type', '').upper()
                val   = desc.get('descriptor', '').strip()
                if not val:
                    continue
                if smiles is None and 'SMILES' in dtype:
                    smiles = val
                if inchi is None and dtype == 'INCHI':
                    inchi = val
        return smiles, inchi
    except Exception as exc:
        logger.warning('%s: API fetch error — %s', ccd_code, exc)
        return None, None


# ─────────────────────────────────────────────────────────────────────────────
# Core resolution
# ─────────────────────────────────────────────────────────────────────────────

def resolve_smiles(ccd_codes: list[str],
                   cif_lookup: dict,
                   session) -> pd.DataFrame:
    """
    Resolve SMILES for each CCD code from CIF (primary) or API (fallback).

    Returns DataFrame: ccd_code, smiles_source, canonical_smiles
    """
    records = []
    n_total = len(ccd_codes)
    counts = {'cif_smiles': 0, 'cif_inchi': 0,
              'api_smiles': 0, 'api_inchi': 0,
              'partial': 0, 'rdkit_failed': 0, 'api_failed': 0}

    for i, ccd in enumerate(ccd_codes, 1):
        if i % 500 == 0 or i == n_total:
            print(f'  {i:,}/{n_total:,} …', flush=True)

        entry = cif_lookup.get(ccd)
        raw_smiles = entry['smiles'] if entry else None
        inchi      = entry['inchi']  if entry else None
        prefix     = 'cif'

        if raw_smiles is None and inchi is None:
            # CIF miss — try API
            raw_smiles, inchi = fetch_smiles_api(ccd, session)
            time.sleep(API_SLEEP)
            prefix = 'api'
            if raw_smiles is None and inchi is None:
                records.append({'ccd_code': ccd, 'smiles_source': 'api_failed',
                                'canonical_smiles': None})
                counts['api_failed'] += 1
                logger.debug('%s: not in CIF and API returned nothing', ccd)
                continue

        canonical, parse_src = sanitise_smiles(raw_smiles, inchi)

        if canonical is None:
            source = 'rdkit_failed'
            logger.warning('%s: all parse attempts failed', ccd)
        elif parse_src == 'partial':
            source = 'partial'
            logger.warning('%s: partial sanitization (kekulize skipped)', ccd)
        else:
            source = f'{prefix}_{parse_src}'   # e.g. 'cif_smiles', 'api_inchi'

        counts[source] = counts.get(source, 0) + 1
        records.append({'ccd_code': ccd, 'smiles_source': source,
                        'canonical_smiles': canonical})

    n_resolved = sum(v for k, v in counts.items()
                     if k not in ('rdkit_failed', 'api_failed'))
    print(f'\n  Resolved from CIF SMILES : {counts["cif_smiles"]:,}')
    print(f'  Resolved from CIF InChI  : {counts["cif_inchi"]:,}')
    print(f'  Resolved from API SMILES : {counts["api_smiles"]:,}')
    print(f'  Resolved from API InChI  : {counts["api_inchi"]:,}')
    print(f'  Partial (no kekulize)    : {counts["partial"]:,}')
    print(f'  RDKit failed             : {counts["rdkit_failed"]:,}')
    print(f'  API failed (CIF miss)    : {counts["api_failed"]:,}')
    print(f'  Total coverage           : {100*n_resolved/max(n_total,1):.1f}%')
    logger.info('Resolved %d/%d; cif_smi=%d cif_inchi=%d api=%d '
                'partial=%d rdkit_fail=%d api_fail=%d',
                n_resolved, n_total,
                counts['cif_smiles'], counts['cif_inchi'],
                counts['api_smiles'] + counts['api_inchi'],
                counts['partial'], counts['rdkit_failed'], counts['api_failed'])
    return pd.DataFrame(records)


# ─────────────────────────────────────────────────────────────────────────────
# Sanity checks
# ─────────────────────────────────────────────────────────────────────────────

def run_sanity_checks(smiles_df: pd.DataFrame) -> None:
    from rdkit import Chem
    print('\n  Sanity checks:')
    for ccd, description, predicate in SANITY_CHECKS:
        row = smiles_df[smiles_df['ccd_code'] == ccd]
        if row.empty:
            print(f'    {ccd}: SKIP (not in dataset)')
            continue
        smi = row.iloc[0]['canonical_smiles']
        if not smi or pd.isna(smi):
            print(f'    {ccd}: FAIL (no SMILES resolved)')
            continue
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            print(f'    {ccd}: FAIL (RDKit parse error on stored SMILES)')
            continue
        status = 'PASS' if predicate(mol) else 'FAIL'
        print(f'    {ccd} [{description}]: {status}')


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description='Resolve canonical SMILES for all ligand CCD codes.')
    cli_utils.add_common_args(parser, force=True)
    args = parser.parse_args()

    t0 = cli_utils.init_run('enrich_ligand_smiles', args, config.LOG_DIR)

    # ── Load CSV ───────────────────────────────────────────────────────────────
    print('\n[1/4] Loading enriched CSV …')
    df = io_utils.load_enriched_csv(config.ENRICHED_CSV)
    print(f'  Loaded {len(df):,} rows.')

    # ── Extract unique CCD codes ───────────────────────────────────────────────
    print('\n[2/4] Extracting unique CCD codes …')
    df = df.dropna(subset=['ligand_code']) # drop rows where ligand_code is NaN
    df = io_utils.extract_ccd_code_rsplit(df)
    unique_ccds = sorted(df['ccd_code'].dropna().unique().tolist())
    print(f'  Unique CCD codes: {len(unique_ccds):,}')

    # ── Output path ────────────────────────────────────────────────────────────
    out_path = os.path.join(config.DATA_DIR, 'ligand_smiles.csv')

    if io_utils.output_is_cached(out_path, args.force):
        print(f'\n[3/4] Output already exists — loading {out_path}')
        print('      (use --force to re-run)')
        smiles_df = pd.read_csv(out_path)
    else:
        # ── Parse CIF ─────────────────────────────────────────────────────────
        print('\n[3/4] Building CIF lookup …')
        if not os.path.exists(config.CCD_CIF):
            print(f'  ERROR: CIF not found at {config.CCD_CIF}')
            sys.exit(1)
        cif_lookup = build_cif_lookup(config.CCD_CIF)

        session = create_session()
        smiles_df = resolve_smiles(unique_ccds, cif_lookup, session)
        smiles_df.to_csv(out_path, index=False)
        print(f'\n  Saved → {out_path}  ({len(smiles_df):,} rows)')
        logger.info('Wrote %d rows to %s', len(smiles_df), out_path)

    # ── Sanity checks ─────────────────────────────────────────────────────────
    print('\n[4/4] Sanity checks …')
    run_sanity_checks(smiles_df)

    cli_utils.finish_run(t0)


if __name__ == '__main__':
    main()
