"""
analysis/src/config.py
Central configuration for the analysis module.

Injects collection/src/ onto sys.path so that utilities (setup_logging,
load_pickle, cif2seq, …) and collection constants (CIF_DIR, MIN_CHAIN_LENGTH)
can be imported directly.
"""

import os
import sys

# ── Path setup ────────────────────────────────────────────────────────────────
_SRC_DIR      = os.path.dirname(os.path.abspath(__file__))
_ANALYSIS_DIR = os.path.dirname(_SRC_DIR)
_DATASET_DIR  = os.path.join(os.path.dirname(_ANALYSIS_DIR), 'collection')
_DATASET_SRC  = os.path.join(_DATASET_DIR, 'src')

if _DATASET_SRC not in sys.path:
    sys.path.insert(0, _DATASET_SRC)

# ── Shared constants imported from collection ─────────────────────────────────
# Both files are named 'config.py', so a plain `from config import …` would
# resolve to this file itself.  Use importlib to load by absolute path instead.
import importlib.util as _ilu

_ds_cfg_path = os.path.join(_DATASET_SRC, 'config.py')
if not os.path.isfile(_ds_cfg_path):
    # spec_from_file_location succeeds for a nonexistent path and only fails at
    # exec_module, deep inside importlib — check here so the message is useful.
    raise RuntimeError(
        f'analysis requires a sibling collection/ checkout: expected {_ds_cfg_path}. '
        f'Keep collection/ and analysis/ under one parent directory (see README).')
_ds_cfg_spec = _ilu.spec_from_file_location('collection_config', _ds_cfg_path)
_ds_cfg = _ilu.module_from_spec(_ds_cfg_spec)
_ds_cfg_spec.loader.exec_module(_ds_cfg)

# CIF_DIR  – path to directory containing {PDB_ID}-assembly1.cif files
# MIN_CHAIN_LENGTH – minimum residues for a chain (used inside cif2seq)
# NUCLEIC_ACIDS – residue/monomer names treated as nucleic acid (not protein)
#
# These three are also a re-export contract, not just a convenience. Because this
# module is itself named 'config', collection's utils.py resolves its own
# `from config import MIN_CHAIN_LENGTH` to *this* file whenever analysis imports it.
# Any name collection/src/utils.py imports from config must therefore be re-exported
# here, or analysis breaks with an unrelated-looking ImportError.
CIF_DIR         = _ds_cfg.CIF_DIR
MIN_CHAIN_LENGTH = _ds_cfg.MIN_CHAIN_LENGTH
NUCLEIC_ACIDS    = _ds_cfg.NUCLEIC_ACIDS

# ── Input dataset (user-supplied copy) ────────────────────────────────────────
INPUT_DATASET_DIR = os.path.join(_ANALYSIS_DIR, 'input_dataset')
ENRICHED_CSV = os.path.join(INPUT_DATASET_DIR,
                             'nonredundant_sites_annotated_enriched.csv')

# ── Paths inside collection/result ────────────────────────────────────────────
RESULT_NONREDUNDANT_DIR = os.path.join(_DATASET_DIR, 'result', 'nonredundant')
PDB_CIF_CACHE = os.path.join(RESULT_NONREDUNDANT_DIR, 'pdb_cif_cache.p')
ENRICHED_PKL  = os.path.join(RESULT_NONREDUNDANT_DIR,
                              'nonredundant_sites_annotated_enriched.p')
# SCOP2/SCOP2B domain-segment mappings (per PDB chain) — collection's
# map_pdb_nextgen_annotations.py output. Used by the CLUSTER stage's SCOP2
# purity check (numeric SF IDs only; names are fetched separately from the
# EBI SCOP2 API, see cluster_protein_similarity.py).
SCOP_MAPPING_CSV = os.path.join(_ds_cfg.RESULT_MAPPING_DIR, 'mapped_scop.csv')
# CATH domain-segment mappings (per PDB chain) — same collection pipeline step as
# SCOP_MAPPING_CSV. Used by the CLUSTER stage's CATH topology-level fallback
# annotation (names are fetched separately from CATH's release site, see
# cluster_protein_similarity.py).
CATH_MAPPING_CSV = os.path.join(_ds_cfg.RESULT_MAPPING_DIR, 'mapped_cath.csv')

# ── analysis-local paths ───────────────────────────────────────────────────
ANALYSIS_DIR = _ANALYSIS_DIR
DATA_DIR   = os.path.join(_ANALYSIS_DIR, 'data')
LOG_DIR    = os.path.join(_ANALYSIS_DIR, 'logs')
CCD_DIR    = os.path.join(_ANALYSIS_DIR, 'ccd')
CCD_CIF    = os.path.join(CCD_DIR, 'Components-rel-alt.cif')

# Ensure local directories exist when the module is first imported
for _d in (DATA_DIR, LOG_DIR, CCD_DIR):
    os.makedirs(_d, exist_ok=True)

# ── Site-category distance thresholds ────────────────────────────────────────
#   orthosteric : min_dist <= ORTHOSTERIC_MAX_DIST
#   allosteric  : min_dist >  ORTHOSTERIC_MAX_DIST
ORTHOSTERIC_MAX_DIST = 5.0

CATEGORIES = ['orthosteric', 'allosteric']

# Arpeggio contact types that are counted, plotted and written to the output CSVs.
# weak_hbond/weak_polar are merged into their base type first (see
# analyze_interactions._CONTACT_TYPE_MERGES); everything else is dropped.
SHOWN_CONTACT_TYPES = ['hydrophobic', 'ionic', 'aromatic', 'hbond', 'polar']


def assign_site_category(min_dist: float) -> str | None:
    """Classify a site by its distance to the nearest functional residue.

    Returns None when min_dist is unknown. A bare `min_dist <= threshold` would
    evaluate NaN as False and label those sites 'allosteric'; collection's
    enriched CSVs carry a null min_dist for every site with no measured distance
    (none in the 'annotated' table, 35,982 of 63,081 in 'complete').
    """
    if min_dist is None or min_dist != min_dist:   # None or NaN
        return None
    return 'orthosteric' if min_dist <= ORTHOSTERIC_MAX_DIST else 'allosteric'


def apply_exclusive_filter(df, data_dir: str, logger):
    """Filter df to ligands exclusive to one site category.

    Reads ligand_exclusive_categories.csv produced by analyze_ligand_descriptors.py.  Adds a
    ccd_code column (parsed from ligand_code) if not already present.
    Returns df unchanged and logs a warning when the file is absent.
    """
    import pandas as pd
    excl_path = os.path.join(data_dir, 'ligand_exclusive_categories.csv')
    if not os.path.exists(excl_path):
        logger.warning('ligand_exclusive_categories.csv missing — run analyze_ligand_descriptors.py first; '
                       'no exclusive filter applied')
        return df
    import io_utils
    excl_ccds = set(pd.read_csv(excl_path)['ccd_code'].dropna())
    if 'ccd_code' not in df.columns:
        df = io_utils.extract_ccd_code_rsplit(df.copy())
    before = len(df)
    df = df[df['ccd_code'].isin(excl_ccds)].copy()
    logger.info('Exclusive filter applied: %d → %d sites', before, len(df))
    print(f'  [exclusive filter] {before - len(df):,} sites removed; {len(df):,} remain.')
    return df

# ── Derived descriptor columns ────────────────────────────────────────────────
# Heavy-atom-normalised variants: (source column, derived column). Lives here so
# analyze_ligand_descriptors (plot CSVs) and the BUILD step (site_ligand_features)
# derive them identically — they used to have separate, drifting copies.
DERIVED_PER_HA_COLUMNS = [
    ('MolLogP',                 'logP_per_HA'),
    ('TPSA',                    'TPSA_per_HA'),
    ('NumHDonors',              'HBD_per_HA'),
    ('NumHAcceptors',           'HBA_per_HA'),
    ('RingCount',               'RingCount_per_HA'),
    ('NumAromaticRings',        'NumAromaticRings_per_HA'),
    ('NumRotatableBonds',       'NumRotatableBonds_per_HA'),
    ('NumAliphaticRings',       'NumAliphaticRings_per_HA'),
    ('NumSaturatedRings',       'NumSaturatedRings_per_HA'),
    ('NumAromaticHeterocycles', 'NumAromaticHeterocycles_per_HA'),
    ('NumRingSystems',          'NumRingSystems_per_HA'),
    ('NumRings5',               'NumRings5_per_HA'),
    ('NumRings6',               'NumRings6_per_HA'),
    ('fr_sulfonamd',            'fr_sulfonamd_per_HA'),
    ('fr_piperzine',            'fr_piperzine_per_HA'),
    ('fr_methoxy',              'fr_methoxy_per_HA'),
    ('fr_halogen',              'fr_halogen_per_HA'),
    ('fr_aniline',              'fr_aniline_per_HA'),
    ('fr_nitro_arom',           'fr_nitro_arom_per_HA'),
    ('fr_ether',                'fr_ether_per_HA'),
    ('fr_Ar_COO',               'fr_Ar_COO_per_HA'),
    ('fr_amide',                'fr_amide_per_HA'),
]


def add_derived_columns(df):
    """Add per-heavy-atom normalised columns and ASA_per_Volume, in place.

    Zero denominators become NaN rather than inf. Source columns that aren't
    present are skipped, so this is safe on a partial descriptor table.
    """
    if 'HeavyAtomCount' not in df.columns:
        return df
    ha = df['HeavyAtomCount'].replace(0, float('nan'))
    for src, dst in DERIVED_PER_HA_COLUMNS:
        if src in df.columns:
            df[dst] = df[src] / ha
    if 'LabuteASA' in df.columns and 'MolVolume' in df.columns:
        df['ASA_per_Volume'] = df['LabuteASA'] / df['MolVolume'].replace(0, float('nan'))
    return df


# ── Statistical parameters ────────────────────────────────────────────────────
# Hypothesis testing (raw per-test Mann-Whitney U) lives downstream.
SIZE_CONFOUND_RHO = 0.3    # |Spearman ρ| above which site-size confounding is flagged

# ── RCSB download URL ─────────────────────────────────────────────────────────
# Used by site_pipeline_utils.resolve_cif_path when a CIF is absent from CIF_DIR
# and the caller supplies a download_dir to cache it into.
RCSB_ASSEMBLY_URL = 'https://files.rcsb.org/download/{pdb_id}-assembly1.cif'

# ── Amino-acid reference data ─────────────────────────────────────────────────

# Standard 20 one-letter codes (excludes 'X' / modified residues)
STANDARD_AA = list('ACDEFGHIKLMNPQRSTVWY')

# ── Columns dropped from site_ligand_features.csv at the BUILD step ──────────
# Either: merge artifacts, descriptors unused downstream, redundant median
# variants, or per-AA fraction columns not consumed downstream.
BUILD_EXCLUDE_COLUMNS = [
    'ccd_code', 'smiles_source', 'parse_quality',
    'NumAromaticAtoms', 'LipinskiViolations', 'embed_source',
    'RotatableBondFraction',
    'n_sites_orthosteric', 'n_sites_allosteric',
    'modal_category', 'exclusive_category',
    'n_res_with_sasa', 'n_res_total_x', 'n_res_total_y',
    'n_res_with_coords', 'chains_in_site',
    'n_universe_residues',
] + [f'frac_{aa}' for aa in STANDARD_AA]

# ── Column order for site_ligand_features.csv, right of 'is_exclusive' ───────
# Grouped for readability: ligand descriptors first, then site geometry and AA
# composition, with each raw/_per_HA sibling pair kept adjacent and columns
# belonging to neither group in a trailing block.
SITE_LIGAND_FEATURE_COLUMN_ORDER = [
    # Block A: ligand descriptors, raw/_per_HA siblings paired
    'ExactMolWt', 'HeavyAtomCount', 'MolLogP', 'logP_per_HA', 'LabuteASA', 'MolVolume',
    'ASA_per_Volume', 'NumRotatableBonds', 'NumRotatableBonds_per_HA',
    'FractionAromaticAtoms', 'FractionCSP3', 'QED', 'TPSA', 'TPSA_per_HA',
    'NumHDonors', 'HBD_per_HA', 'NumHAcceptors', 'HBA_per_HA',
    'RingCount', 'RingCount_per_HA',
    'NumAromaticRings', 'NumAromaticRings_per_HA',
    'NumAliphaticRings', 'NumAliphaticRings_per_HA',
    'NumSaturatedRings', 'NumSaturatedRings_per_HA',
    'NumAromaticHeterocycles', 'NumAromaticHeterocycles_per_HA',
    'NumRingSystems', 'NumRingSystems_per_HA',
    'NumRings5', 'NumRings5_per_HA',
    'NumRings6', 'NumRings6_per_HA',
    # Block B: site geometry + AA composition
    'site_size', 'mean_num_neighbors', 'rSASA_apo', 'rSASA_bound',
    'frac_surface_apo', 'frac_intermediate_apo',
    'Charged', 'Polar uncharged',
    'Hydrophobic aliphatic', 'Hydrophobic aromatic',
    # Block C: orphans (no plotted sibling), original relative order
    'largest_ring_size',
    'fr_sulfonamd', 'fr_sulfonamd_per_HA', 'fr_piperzine', 'fr_piperzine_per_HA',
    'fr_methoxy', 'fr_methoxy_per_HA', 'fr_halogen', 'fr_halogen_per_HA',
    'fr_aniline', 'fr_aniline_per_HA', 'fr_nitro_arom', 'fr_nitro_arom_per_HA',
    'fr_ether', 'fr_ether_per_HA', 'fr_Ar_COO', 'fr_Ar_COO_per_HA',
    'fr_amide', 'fr_amide_per_HA',
    'has_metal', 'NPR1', 'NPR2',
    'frac_hydrophobic', 'frac_ionic', 'frac_aromatic', 'frac_hbond', 'frac_polar',
]

# Kyte-Doolittle hydrophobicity scale
KD_SCALE = {
    'A':  1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C':  2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I':  4.5,
    'L':  3.8, 'K': -3.9, 'M':  1.9, 'F':  2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V':  4.2,
}
  
# Physicochemical property groups
AA_PROPERTY_GROUPS = {
    'Charged':               ['K', 'R', 'D', 'E'],
    'Polar uncharged':       ['S', 'T', 'N', 'Q', 'H', 'Y', 'G', 'A'],
    'Hydrophobic aliphatic': ['V', 'I', 'L', 'M', 'P', 'C'],
    'Hydrophobic aromatic':  ['F', 'W'],
}

# Reverse map: one-letter code → property group name
AA_TO_GROUP = {
    aa: grp
    for grp, aas in AA_PROPERTY_GROUPS.items()
    for aa in aas
}
