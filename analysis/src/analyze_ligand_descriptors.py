"""
analyze_ligand_descriptors.py
============================================
Compute RDKit molecular descriptors for each unique ligand (CCD code), tabulated
by binding-site category (orthosteric vs allosteric).

Descriptors are emitted at the ligand level (avoids pseudoreplication): only
ligands exclusive to one category (orthosteric or allosteric, never both) are
used.  Site counts and exclusive-category assignments are loaded from Step 0
outputs (ligand_site_counts.csv, ligand_exclusive_categories.csv), with an inline
fallback if those files are absent (standalone use without the pipeline).

Inputs
------
  data/ligand_smiles.csv             (from enrich_ligand_smiles.py)
  data/ligand_site_counts.csv        (from run_analysis_pipeline.py Step 0)
  data/ligand_exclusive_categories.csv (from run_analysis_pipeline.py Step 0)
  nonredundant_sites_annotated_enriched.csv

Outputs
-------
  data/ligand_descriptors.csv          (per unique CCD code)
  data/ligand_uniqueness_report.txt    (exclusive-ligand counts)
  data/plot_ligand_descriptors.csv     (long format: ccd_code, exclusive_category, metric, value)
  data/plot_ligand_npr.csv             (wide format: NPR1/NPR2 columns for scatter plots)

Usage
-----
    python analyze_ligand_descriptors.py [--verbose]
"""

import argparse
import logging
import multiprocessing as mp
import os
import sys
import time

import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, Descriptors, Fragments, rdMolDescriptors, QED
rdBase.DisableLog('rdApp.*')   # silence RDKit C++ console output globally
from rdkit.Chem import Descriptors3D
from rdkit.Geometry import rdGeometry

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import config  # must precede cli_utils/io_utils: injects collection/src onto sys.path
import cli_utils
import io_utils
from enrich_ligand_smiles import build_cif_lookup

logger = logging.getLogger(__name__)

CATEGORIES = config.CATEGORIES

EXCLUDE_SMILES_SOURCES = {'rdkit_failed', 'api_failed', 'partial'}
EXCLUDE_PARSE_QUALITY  = {'partial', 'failed'}


# ─────────────────────────────────────────────────────────────────────────────
# Worker pool shared state (set once per process via initializer)
# ─────────────────────────────────────────────────────────────────────────────

_worker_cif_lookup: dict = {}


def _worker_init(cif_lookup: dict) -> None:
    """Pool initializer: stash the CIF lookup so _worker_compute needn't be passed it."""
    global _worker_cif_lookup
    _worker_cif_lookup = cif_lookup
    rdBase.DisableLog('rdApp.*')  # re-silence in each worker process


def _worker_compute(args) -> dict:
    """compute_descriptors for one ligand; args is (ccd_code, smiles, smiles_source)."""
    ccd_code, smiles, smiles_source = args
    return compute_descriptors(ccd_code, smiles, _worker_cif_lookup,
                               smiles_source=smiles_source)




# ─────────────────────────────────────────────────────────────────────────────
# Descriptor computation
# ─────────────────────────────────────────────────────────────────────────────

_METAL_SYMBOLS = frozenset([
    'Li', 'Na', 'K',  'Rb', 'Cs', 'Be', 'Mg', 'Ca', 'Sr', 'Ba',
    'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Al', 'Ga', 'In', 'Sn', 'Tl', 'Pb', 'Bi',
    'La', 'Ce', 'Nd', 'Eu', 'Gd', 'Tb', 'Yb', 'Lu',
])

_EMBED_KWARGS = dict(randomSeed=42, useExpTorsionAnglePrefs=True,
                    useBasicKnowledge=True, ETversion=2, maxAttempts=50)

_BOND_TYPE_MAP = {
    'SING': Chem.BondType.SINGLE,
    'DOUB': Chem.BondType.DOUBLE,
    'TRIP': Chem.BondType.TRIPLE,
    'AROM': Chem.BondType.AROMATIC,
}


def _embed_3d(ccd: str, mol, cif_lookup: dict):
    """Give `mol` a 3D conformer, preferring the CCD's ideal coordinates.

    Needed for NPR1/NPR2 and volume. The CIF path rebuilds connectivity from the
    CCD bond table so atom order matches the coordinate list, then assigns
    positions directly — no force field, and reproducible. Falls back to ETKDG
    when CIF coordinates are missing or incomplete.

    Returns (mol_with_conformer, 'cif'|'etkdg') or (None, 'failed').
    """
    entry = cif_lookup.get(ccd, {})
    atom_rows = entry.get('atoms', [])
    bond_rows  = entry.get('bonds', [])

    if atom_rows and bond_rows:
        try:
            cif_rwmol = Chem.RWMol()
            atom_id_to_idx = {}
            for atom_id, element, x, y, z in atom_rows:
                atom_id_to_idx[atom_id] = cif_rwmol.AddAtom(Chem.Atom(element))
            for a1, a2, order in bond_rows:
                i, j = atom_id_to_idx.get(a1), atom_id_to_idx.get(a2)
                if i is not None and j is not None:
                    cif_rwmol.AddBond(i, j, _BOND_TYPE_MAP.get(order, Chem.BondType.SINGLE))
            Chem.SanitizeMol(cif_rwmol)

            # Assign CCD ideal coords directly — atom_id order is deterministic
            # and matches the CIF atom table
            n_atoms = cif_rwmol.GetNumAtoms()
            has_coords = sum(1 for _, _, x, y, z in atom_rows if x is not None)
            if has_coords / max(n_atoms, 1) >= 0.8:
                conf = Chem.Conformer(n_atoms)
                for i, (atom_id, _, x, y, z) in enumerate(atom_rows):
                    idx = atom_id_to_idx[atom_id]
                    if x is not None:
                        conf.SetAtomPosition(idx, rdGeometry.Point3D(float(x), float(y), float(z)))
                    else:
                        # Atom missing ideal coords: place at origin as a safe placeholder
                        conf.SetAtomPosition(idx, rdGeometry.Point3D(0.0, 0.0, 0.0))
                cif_rwmol.AddConformer(conf, assignId=True)
                return cif_rwmol.GetMol(), 'cif'
        except Exception:
            pass

    mol_h = Chem.AddHs(mol)
    res = AllChem.EmbedMolecule(mol_h, **_EMBED_KWARGS)
    return (mol_h, 'etkdg') if res == 0 else (None, 'failed')


def _num_ring_systems(mol) -> int:
    """Count ring systems (SSSR rings merged when they share an atom)."""
    atom_rings = mol.GetRingInfo().AtomRings()
    if not atom_rings:
        return 0
    parent = {}

    def find(x):
        root = x
        while parent.get(root, root) != root:
            root = parent[root]
        while parent.get(x, x) != root:
            parent[x], x = root, parent.get(x, x)
        return root

    def union(x, y):
        parent[find(x)] = find(y)

    for ring in atom_rings:
        for atom in ring[1:]:
            union(ring[0], atom)

    ring_atoms = {a for ring in atom_rings for a in ring}
    return len({find(a) for a in ring_atoms})


def compute_descriptors(ccd_code: str, smiles: str, cif_lookup: dict,
                        smiles_source: str = '') -> dict:
    """All RDKit descriptors for one CCD code.

    Parsing tries strict sanitisation, then partial (no kekulization); which one
    succeeded is recorded in parse_quality. The 3D descriptors (NPR1, NPR2,
    MolVolume) come from _embed_3d's conformer and stay NaN if embedding failed.

    Always returns at least ccd_code and parse_quality.
    """
    import warnings
    rec = {'ccd_code': ccd_code, 'smiles_source': smiles_source}

    parse_quality = 'ok'
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        try:
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is not None:
                Chem.SanitizeMol(mol,
                                 Chem.SanitizeFlags.SANITIZE_ALL ^
                                 Chem.SanitizeFlags.SANITIZE_KEKULIZE)
                parse_quality = 'partial'
                logger.warning('%s [%s]: partial sanitization (kekulize skipped) — '
                               'descriptors may be unreliable', ccd_code, smiles_source)
        except Exception:
            mol = None
    if mol is None:
        rec['parse_quality'] = 'failed'
        logger.warning('%s [%s]: could not parse SMILES — excluded from analysis',
                       ccd_code, smiles_source)
        return rec

    rec['parse_quality'] = parse_quality

    try:
        rec['ExactMolWt']     = Descriptors.ExactMolWt(mol)
        rec['HeavyAtomCount'] = mol.GetNumHeavyAtoms()

        rec['MolLogP'] = Descriptors.MolLogP(mol)

        rec['TPSA']          = Descriptors.TPSA(mol)
        rec['LabuteASA']     = rdMolDescriptors.CalcLabuteASA(mol)
        rec['NumHDonors']    = rdMolDescriptors.CalcNumHBD(mol)
        rec['NumHAcceptors'] = rdMolDescriptors.CalcNumHBA(mol)

        n_heavy = mol.GetNumHeavyAtoms() or 1
        rec['NumAromaticRings']      = rdMolDescriptors.CalcNumAromaticRings(mol)
        n_arom_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
        rec['NumAromaticAtoms']      = n_arom_atoms
        rec['FractionAromaticAtoms'] = n_arom_atoms / n_heavy

        rec['FractionCSP3'] = rdMolDescriptors.CalcFractionCSP3(mol)

        n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
        rec['NumRotatableBonds']    = n_rot
        rec['RotatableBondFraction'] = n_rot / n_heavy

        rec['RingCount'] = rdMolDescriptors.CalcNumRings(mol)
        rec['NumAliphaticRings']      = rdMolDescriptors.CalcNumAliphaticRings(mol)
        rec['NumSaturatedRings']      = rdMolDescriptors.CalcNumSaturatedRings(mol)
        rec['NumAromaticHeterocycles'] = rdMolDescriptors.CalcNumAromaticHeterocycles(mol)
        rec['NumRingSystems']         = _num_ring_systems(mol)

        ring_sizes = [len(r) for r in mol.GetRingInfo().AtomRings()]
        rec['largest_ring_size'] = max(ring_sizes, default=0)
        rec['NumRings5'] = ring_sizes.count(5)
        rec['NumRings6'] = ring_sizes.count(6)

        rec['fr_sulfonamd']  = Fragments.fr_sulfonamd(mol)
        rec['fr_piperzine']  = Fragments.fr_piperzine(mol)
        rec['fr_methoxy']    = Fragments.fr_methoxy(mol)
        rec['fr_halogen']    = Fragments.fr_halogen(mol)
        rec['fr_aniline']    = Fragments.fr_aniline(mol)
        rec['fr_nitro_arom'] = Fragments.fr_nitro_arom(mol)
        rec['fr_ether']      = Fragments.fr_ether(mol)
        rec['fr_Ar_COO']     = Fragments.fr_Ar_COO(mol)
        rec['fr_amide']      = Fragments.fr_amide(mol)

        rec['QED'] = QED.qed(mol)
        mw   = rec['ExactMolWt']
        logp = rec['MolLogP']
        hbd  = rec['NumHDonors']
        hba  = rec['NumHAcceptors']
        rec['LipinskiViolations'] = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)

        rec['has_metal'] = bool({a.GetSymbol() for a in mol.GetAtoms()} & _METAL_SYMBOLS)

        # ── 3D shape (CIF coordMap → ETKDG fallback) ──────────────────────────
        mol_3d, embed_src = _embed_3d(ccd_code, mol, cif_lookup)
        rec['embed_source'] = embed_src
        if mol_3d is not None:
            try:
                rec['NPR1'] = Descriptors3D.NPR1(mol_3d)
                rec['NPR2'] = Descriptors3D.NPR2(mol_3d)
            except Exception as exc:
                logger.debug('%s: NPR descriptors failed — %s', ccd_code, exc)
                rec['NPR1'] = float('nan')
                rec['NPR2'] = float('nan')
            try:
                rec['MolVolume'] = AllChem.ComputeMolVolume(mol_3d)
            except Exception as exc:
                logger.debug('%s: MolVolume failed — %s', ccd_code, exc)
                rec['MolVolume'] = float('nan')
        else:
            rec['NPR1']      = float('nan')
            rec['NPR2']      = float('nan')
            rec['MolVolume'] = float('nan')

    except Exception as exc:
        logger.warning('%s: descriptor computation error — %s', ccd_code, exc)

    return rec


def _flush_cache(new_records: list, cached_df: pd.DataFrame, cache_path: str) -> None:
    """Rewrite the cache CSV as cached_df + new_records, via write-then-rename.

    The .tmp-then-os.replace dance keeps a crash mid-write from truncating an
    existing cache.
    """
    new_df = pd.DataFrame(new_records)
    combined = pd.concat([cached_df, new_df], ignore_index=True) if not cached_df.empty else new_df
    tmp = cache_path + '.tmp'
    combined.to_csv(tmp, index=False)
    os.replace(tmp, cache_path)


# Sentinel object used to tell build_descriptor_table to skip disk caching.
# Use a unique object() instead of None or a string so that callers
# can safely pass cache_path=_SKIP_CACHE without accidentally matching
# any real path or falsy value.
_SKIP_CACHE = object()


def build_descriptor_table(smiles_df: pd.DataFrame, cif_lookup: dict,
                           cache_path=None) -> pd.DataFrame:
    """Descriptors for every CCD code with a valid SMILES, one row each.

    Ligands with no SMILES still appear, with parse_quality='failed' and NaN
    descriptors. Parallelised across all cores above 200 ligands.

    Results are cached by ccd_code in data/descriptor_cache.csv and flushed
    every 1000 new rows; a cache missing any currently-produced column is
    treated as stale and discarded. Pass _SKIP_CACHE to disable caching (--force).
    """
    if cache_path is None:
        cache_path = os.path.join(config.DATA_DIR, 'descriptor_cache.csv')
    use_cache = cache_path is not _SKIP_CACHE
    valid = smiles_df.dropna(subset=['canonical_smiles']).copy()
    invalid = smiles_df[smiles_df['canonical_smiles'].isna()]
    if 'smiles_source' not in valid.columns:
        valid = valid.assign(smiles_source='')

    expected_cols = set(ALL_DESCRIPTOR_COLS)
    cached_df = pd.DataFrame()
    hit_ccds: set = set()
    if use_cache and os.path.exists(cache_path):
        header_cols = set(pd.read_csv(cache_path, nrows=0).columns)
        if expected_cols.issubset(header_cols):
            cached_df = pd.read_csv(cache_path, low_memory=False)
            hit_ccds = set(cached_df['ccd_code'])
        else:
            missing = expected_cols - header_cols
            print(f'  Cache stale ({len(missing)} new descriptor columns) — recomputing all …')
            logger.info('Cache stale: missing columns %s', sorted(missing))

    todo = valid[~valid['ccd_code'].isin(hit_ccds)]
    n_cached = len(valid) - len(todo)
    if n_cached:
        print(f'  Cache: {n_cached:,} CCD codes loaded, {len(todo):,} to compute …')
        logger.info('Cache hit: %d/%d CCD codes; %d to compute',
                    n_cached, len(valid), len(todo))

    print(f'  Computing descriptors for {len(todo):,} CCD codes '
          f'({len(invalid):,} skipped — no SMILES) …')
    logger.info('Descriptor computation: %d attempted, %d skipped (no SMILES)',
                len(todo), len(invalid))

    records: list[dict] = []
    n_total = len(todo)
    n_3d_cif = n_3d_etkdg = n_3d_fail = 0
    n_workers = max(1, os.cpu_count() or 1)
    use_parallel = n_total >= 200
    # flushed_df tracks what is on disk; `records` holds only the batch since the last
    # flush. Each flush rewrites the whole CSV to keep a single consistent column order.
    flushed_df = cached_df

    tasks = list(zip(todo['ccd_code'], todo['canonical_smiles'], todo['smiles_source']))

    if use_parallel:
        print(f'  Using {n_workers} parallel workers …')
    t_start = time.time()

    def _collect(rec: dict, done: int) -> None:
        nonlocal n_3d_cif, n_3d_etkdg, n_3d_fail, flushed_df
        records.append(rec)
        src = rec.get('embed_source', 'failed')
        if src == 'cif':
            n_3d_cif += 1
        elif src == 'etkdg':
            n_3d_etkdg += 1
        else:
            n_3d_fail += 1
        if use_cache and done % 1000 == 0:
            _flush_cache(records, flushed_df, cache_path)
            flushed_df = (pd.concat([flushed_df, pd.DataFrame(records)], ignore_index=True)
                          if not flushed_df.empty else pd.DataFrame(records))
            records.clear()
            logger.info('Cache flushed at %d records', n_cached + done)
        if done % 100 == 0 or done == n_total:
            elapsed = time.time() - t_start
            rate = done / elapsed if elapsed > 0 else 0
            eta = (n_total - done) / rate if rate > 0 else 0
            print(f'\r  {done:,}/{n_total:,}  '
                  f'{rate:.0f} lig/s  ETA {eta/60:.1f} min   ',
                  end='', flush=True)

    if n_total > 0:
        if use_parallel:
            chunksize = max(1, n_total // (n_workers * 10))
            with mp.Pool(n_workers, initializer=_worker_init, initargs=(cif_lookup,)) as pool:
                for done, rec in enumerate(
                        pool.imap_unordered(_worker_compute, tasks, chunksize=chunksize), 1):
                    _collect(rec, done)
        else:
            for done, (ccd_code, smiles, smiles_source) in enumerate(tasks, 1):
                _collect(compute_descriptors(ccd_code, smiles, cif_lookup,
                                             smiles_source=smiles_source), done)
        print()  # newline after progress line

    if use_cache and records:
        _flush_cache(records, flushed_df, cache_path)
        logger.info('Cache written: %d new rows → %s', n_total, cache_path)

    # Build full table in memory from flushed_df (everything periodically flushed)
    # and any remaining records not yet flushed.
    new_df = pd.DataFrame(records) if records else pd.DataFrame()
    if not flushed_df.empty and not new_df.empty:
        full_df = pd.concat([flushed_df, new_df], ignore_index=True)
    elif not new_df.empty:
        full_df = new_df
    else:
        full_df = flushed_df
    desc_df = full_df[full_df['ccd_code'].isin(valid['ccd_code'])].copy()

    inv_src = invalid['smiles_source'] if 'smiles_source' in invalid.columns else [''] * len(invalid)
    invalid_rows = [{'ccd_code': ccd, 'smiles_source': src, 'parse_quality': 'failed'}
                    for ccd, src in zip(invalid['ccd_code'], inv_src)]
    if invalid_rows:
        desc_df = pd.concat([desc_df, pd.DataFrame(invalid_rows)], ignore_index=True)

    pct_3d_fail = 100 * n_3d_fail / max(n_total, 1)
    if n_total > 0:
        print(f'  3D embedding: CIF coords={n_3d_cif:,}  ETKDG fallback={n_3d_etkdg:,}  '
              f'failed={n_3d_fail:,} ({pct_3d_fail:.1f}%)')
        logger.info('3D embedding: cif=%d etkdg=%d failed=%d (%.1f%%)',
                    n_3d_cif, n_3d_etkdg, n_3d_fail, pct_3d_fail)
    return desc_df


# ─────────────────────────────────────────────────────────────────────────────
# Descriptor metadata
# ─────────────────────────────────────────────────────────────────────────────

# Descriptors stored in cache (raw, before derived columns)
_RAW_DESCRIPTORS = [
    'ExactMolWt', 'HeavyAtomCount',
    'MolLogP',
    'TPSA', 'LabuteASA', 'NumHDonors', 'NumHAcceptors',
    'NumAromaticRings', 'NumAromaticAtoms', 'FractionAromaticAtoms',
    'FractionCSP3',
    'NumRotatableBonds', 'RotatableBondFraction',
    'RingCount', 'NumAliphaticRings', 'NumSaturatedRings',
    'NumAromaticHeterocycles', 'NumRingSystems',
    'largest_ring_size', 'NumRings5', 'NumRings6',
    'fr_sulfonamd', 'fr_piperzine', 'fr_methoxy', 'fr_halogen',
    'fr_aniline', 'fr_nitro_arom', 'fr_ether', 'fr_Ar_COO', 'fr_amide',
    'QED', 'LipinskiViolations',
    'NPR1', 'NPR2', 'MolVolume',
    'has_metal',
]

# Complete set of columns compute_descriptors() can return.
ALL_DESCRIPTOR_COLS = [
    'ccd_code', 'smiles_source', 'parse_quality', 'embed_source',
    *_RAW_DESCRIPTORS,
]

ALL_DESCRIPTORS = [
    # Size, lipophilicity, surface area, volume
    'ExactMolWt', 'HeavyAtomCount', 'MolLogP', 'logP_per_HA',
    'LabuteASA', 'MolVolume', 'ASA_per_Volume',
    # Ring complexity and flexibility
    'RingCount', 'RingCount_per_HA',
    'NumAromaticRings', 'NumAromaticRings_per_HA',
    'NumAliphaticRings', 'NumAliphaticRings_per_HA',
    'NumSaturatedRings', 'NumSaturatedRings_per_HA',
    'NumAromaticHeterocycles', 'NumAromaticHeterocycles_per_HA',
    'NumRingSystems', 'NumRingSystems_per_HA',
    'NumRings5', 'NumRings5_per_HA', 'NumRings6', 'NumRings6_per_HA',
    'largest_ring_size',
    'NumRotatableBonds', 'NumRotatableBonds_per_HA',
    'NumAromaticAtoms', 'FractionAromaticAtoms', 'FractionCSP3', 'QED',
    # Functional groups (RDKit Fragments SMARTS counters)
    'fr_sulfonamd', 'fr_piperzine', 'fr_methoxy', 'fr_halogen',
    'fr_aniline', 'fr_nitro_arom', 'fr_ether', 'fr_Ar_COO', 'fr_amide',
    # Polarity
    'TPSA', 'TPSA_per_HA', 'NumHDonors', 'HBD_per_HA', 'NumHAcceptors', 'HBA_per_HA',
    # Shape (NPR scatter — wide format in plot_ligand_npr.csv)
    'NPR1', 'NPR2',
]



# ─────────────────────────────────────────────────────────────────────────────
# Plot data export
# ─────────────────────────────────────────────────────────────────────────────

def write_plot_csv(desc_df: pd.DataFrame, data_dir: str) -> str:
    """Write the two plot-ready CSVs.

    plot_ligand_descriptors.csv is long (ccd_code, exclusive_category, metric,
    value) over ALL_DESCRIPTORS minus the NPR columns; plot_ligand_npr.csv stays
    wide because the scatter needs both NPR axes at once.
    """
    # Violin descriptors: all except the NPR columns (those get their own scatter CSV)
    violin_descs = [d for d in ALL_DESCRIPTORS
                    if d in desc_df.columns and not d.startswith('NPR')]

    # Melt wide → long so each row is one (ligand, metric, value) observation.
    # Long format is easier for plotting tools to slice by metric name.
    id_vars = ['ccd_code', 'exclusive_category']
    long_df = desc_df[id_vars + violin_descs].melt(
        id_vars=id_vars, var_name='metric', value_name='value'
    )
    long_df = long_df.dropna(subset=['value'])

    out_violin = os.path.join(data_dir, 'plot_ligand_descriptors.csv')
    long_df.to_csv(out_violin, index=False)
    print(f'  Plot data (violin) → {out_violin}')

    # NPR scatter: keep wide format since the scatter uses both axes simultaneously
    npr_cols = [c for c in ['ccd_code', 'exclusive_category', 'NPR1', 'NPR2']
                if c in desc_df.columns]
    npr_subset = [c for c in ['NPR1', 'NPR2'] if c in desc_df.columns]
    npr_df = desc_df[npr_cols].dropna(subset=npr_subset, how='any') if npr_subset else desc_df[npr_cols]
    out_npr = os.path.join(data_dir, 'plot_ligand_npr.csv')
    npr_df.to_csv(out_npr, index=False)
    print(f'  Plot data (NPR scatter) → {out_npr}')

    return out_violin


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    """Run the full ligand-descriptor pipeline: load inputs, compute descriptors, run stats, export.

    Steps:
        1. Load the enriched site CSV and derive ``ccd_code`` /
           ``site_category`` columns.
        2. Load ``ligand_smiles.csv`` (produced by enrich_ligand_smiles.py).
        3. Load the CCD CIF lookup for ideal-coordinate 3D embedding.
        4. Compute RDKit descriptors for every CCD code (with disk cache).
        5. Load site counts + exclusive categories from Step 0 outputs
           (or compute inline if running standalone); merge into descriptor
           table; apply SMILES/quality exclusions; restrict analysis set to
           exclusive ligands; write uniqueness report.
        6. Add heavy-atom-count-normalised derived columns.
        7. Save ``ligand_descriptors.csv`` (exclusive_category dropped —
           per-analysis context, not a descriptor).
        8. Write plot-ready CSVs (``plot_ligand_descriptors.csv``,
           ``plot_ligand_npr.csv``).
    """
    parser = argparse.ArgumentParser(
        description='Compute ligand descriptors by site category.')
    cli_utils.add_common_args(parser, force=True)
    args = parser.parse_args()

    t0 = cli_utils.init_run('analyze_ligand_descriptors', args, config.LOG_DIR)

    # ── [1/8] Load inputs ──────────────────────────────────────────────────────
    print('\n[1/8] Loading enriched CSV …')
    df = io_utils.load_enriched_csv(config.ENRICHED_CSV, dropna_subset=['ligand_code', 'min_dist'])
    # ligand_code format is "{chain}.{resnum}:{ccd}"; rsplit at the last colon
    # so CCD codes that happen to contain a colon are not accidentally split.
    df = io_utils.extract_ccd_code_rsplit(df)
    df = io_utils.assign_site_category_col(df, config.assign_site_category)
    print(f'  {len(df):,} sites loaded.')
    for cat in CATEGORIES:
        n = (df['site_category'] == cat).sum()
        print(f'    {cat:<14}: {n:,} sites')

    # ── [2/8] Load SMILES ──────────────────────────────────────────────────────
    print('\n[2/8] Loading ligand_smiles.csv …')
    smiles_path = os.path.join(config.DATA_DIR, 'ligand_smiles.csv')
    if not os.path.exists(smiles_path):
        print(f'  ERROR: {smiles_path} not found. Run enrich_ligand_smiles.py first.')
        sys.exit(1)
    smiles_df = pd.read_csv(smiles_path)
    n_with_smiles = smiles_df['canonical_smiles'].notna().sum()
    print(f'  {len(smiles_df):,} unique CCD codes; '
          f'{n_with_smiles:,} with valid SMILES.')
    logger.info('%d unique CCD codes; %d with valid SMILES (%d skipped — no SMILES)',
                len(smiles_df), n_with_smiles, len(smiles_df) - n_with_smiles)

    # ── [3/8] Load CIF lookup ──────────────────────────────────────────────────
    print('\n[3/8] Loading CCD CIF lookup …')
    t_cif = time.time()
    cif_lookup = build_cif_lookup(config.CCD_CIF)
    print(f'  CIF lookup: {len(cif_lookup):,} entries loaded in {time.time()-t_cif:.1f}s')
    logger.info('CIF lookup: %d entries', len(cif_lookup))

    # ── [4/8] Compute descriptors ──────────────────────────────────────────────
    print('\n[4/8] Computing RDKit descriptors …')
    cache_path = _SKIP_CACHE if args.force else os.path.join(config.DATA_DIR, 'descriptor_cache.csv')
    desc_df = build_descriptor_table(smiles_df, cif_lookup, cache_path=cache_path)
    logger.info('Descriptor table: %d rows × %d cols', *desc_df.shape)

    # ── [5/8] Load site counts + exclusivity from Step 0 (or compute inline) ───
    print('\n[5/8] Loading site counts and exclusive categories …')
    counts_path = os.path.join(config.DATA_DIR, 'ligand_site_counts.csv')
    excl_path   = os.path.join(config.DATA_DIR, 'ligand_exclusive_categories.csv')

    if os.path.exists(counts_path) and os.path.exists(excl_path):
        site_counts = pd.read_csv(counts_path)
        excl_cats   = pd.read_csv(excl_path)
        print(f'  Loaded from Step 0: {len(excl_cats):,} exclusive CCD codes')
        logger.info('Loaded Step 0 outputs: %d site-count rows, %d exclusive CCD codes',
                    len(site_counts), len(excl_cats))
    else:
        print('  [Step 0 files not found] Computing exclusive categories inline …')
        logger.warning('Step 0 outputs missing — computing inline (run_analysis_pipeline.py Step 0 '
                       'should normally produce these first)')
        counts_pivot = (df.groupby(['ccd_code', 'site_category'])
                          .size().unstack(fill_value=0).reset_index())
        for cat in CATEGORIES:
            if cat not in counts_pivot.columns:
                counts_pivot[cat] = 0

        # excl_mask/excl_cats work on the raw category columns (e.g. 'orthosteric');
        # site_counts is just those same columns renamed for the output CSV.
        excl_mask = (counts_pivot[CATEGORIES] > 0).sum(axis=1) == 1
        excl_cats = (counts_pivot[excl_mask]
                     .assign(exclusive_category=lambda x: x[CATEGORIES].idxmax(axis=1))
                     [['ccd_code', 'exclusive_category']])
        site_counts = counts_pivot[['ccd_code'] + CATEGORIES].rename(
            columns={cat: f'n_sites_{cat}' for cat in CATEGORIES})

    desc_df = desc_df.merge(site_counts, on='ccd_code', how='left')
    for col in ['n_sites_orthosteric', 'n_sites_allosteric']:
        if col in desc_df.columns:
            desc_df[col] = desc_df[col].fillna(0).astype(int)
    desc_df = desc_df.merge(excl_cats, on='ccd_code', how='left')

    # Drop ligands with low-confidence SMILES (failed/partial lookups),
    # low-confidence parsing (failed/partial sanitization), or metals/organometallics
    # (whose logP, TPSA, HBD/HBA values are chemically meaningless).
    mask_src   = desc_df['smiles_source'].isin(EXCLUDE_SMILES_SOURCES)
    mask_qual  = desc_df['parse_quality'].isin(EXCLUDE_PARSE_QUALITY)
    mask_metal = (desc_df['has_metal'].fillna(False)
                  if 'has_metal' in desc_df.columns
                  else pd.Series(False, index=desc_df.index))
    mask_bad   = mask_src | mask_qual | mask_metal
    n_excl = int(mask_bad.sum())
    if n_excl:
        n_src, n_qual, n_metal = int(mask_src.sum()), int(mask_qual.sum()), int(mask_metal.sum())
        print(f'  Excluded {n_excl:,} low-confidence/inorganic ligands '
              f'(smiles_source={n_src}, parse_quality={n_qual}, metal={n_metal})')
        logger.info('Excluded %d ligands: smiles_source=%d parse_quality=%d metal=%d',
                    n_excl, n_src, n_qual, n_metal)

    # desc_df_conf: full confident set (kept for the uniqueness report below).
    # desc_df_analysis: further restricted to ligands exclusive to one site category
    desc_df_conf    = desc_df[~mask_bad].copy()
    desc_df_analysis = desc_df_conf[desc_df_conf['exclusive_category'].notna()].copy()

    n_ortho = int((desc_df_analysis['exclusive_category'] == 'orthosteric').sum())
    n_allo  = int((desc_df_analysis['exclusive_category'] == 'allosteric').sum())
    # Compute true overlap (appears in ≥2 categories) separately from "not in any site".
    count_cols = [f'n_sites_{cat}' for cat in CATEGORIES if f'n_sites_{cat}' in desc_df_conf.columns]
    if count_cols:
        site_presence = (desc_df_conf[count_cols] > 0).sum(axis=1)
        n_overlap = int((site_presence > 1).sum())
        n_no_site = int((site_presence == 0).sum())
    else:
        n_overlap = int(desc_df_conf['exclusive_category'].isna().sum())
        n_no_site = 0
    print(f'  Orthosteric-only ligands : {n_ortho:,}')
    print(f'  Allosteric-only ligands  : {n_allo:,}')
    print(f'  Overlap (found in both)  : {n_overlap:,}  [excluded from all downstream analyses]')
    if n_no_site:
        print(f'  Not in any site          : {n_no_site:,}  [also excluded]')

    report_path = os.path.join(config.DATA_DIR, 'ligand_uniqueness_report.txt')
    with open(report_path, 'w') as fh:
        fh.write('Ligand category uniqueness report\n')
        fh.write('=' * 50 + '\n')
        fh.write(f'Orthosteric-only ligands : {n_ortho:,}\n')
        fh.write(f'Allosteric-only ligands  : {n_allo:,}\n')
        fh.write(f'Overlap (found in both)  : {n_overlap:,}  '
                 f'[excluded from all downstream analyses]\n')
        if n_no_site:
            fh.write(f'Not in any site          : {n_no_site:,}  [also excluded]\n')
    print(f'  Saved → {report_path}')
    logger.info('Uniqueness report written to %s', report_path)

    # ── [6/8] Add derived (HA-normalised) columns ─────────────────────────────
    print('\n[6/8] Computing derived (HA-normalised) columns …')
    desc_df_analysis = config.add_derived_columns(desc_df_analysis)

    # ── [7/8] Save output tables ───────────────────────────────────────────────
    print('\n[7/8] Saving output tables …')
    desc_out = os.path.join(config.DATA_DIR, 'ligand_descriptors.csv')
    desc_df.drop(columns=['exclusive_category']).to_csv(desc_out, index=False)
    print(f'  Saved → {desc_out}  ({len(desc_df):,} rows, full; '
          f'{len(desc_df_analysis):,} in analysis set)')
    logger.info('Wrote %s: %d rows total, %d in analysis set',
                os.path.basename(desc_out), len(desc_df), len(desc_df_analysis))

    # ── [8/8] Writing plot data ────────────────────────────────────────────────
    print('\n[8/8] Writing plot data …')
    write_plot_csv(desc_df_analysis, config.DATA_DIR)

    cli_utils.finish_run(t0)


if __name__ == '__main__':
    main()
