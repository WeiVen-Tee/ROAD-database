"""
Central configuration for the collection pipeline.

Holds every output path under collection/result/, the location of the local
PDB assembly CIF cache, and the tunable constants each step reads (contact
cutoffs, clustering thresholds, service URLs). Grouped by the step that uses
them; anything shared by several steps sits in the block at the top.

CIF_DIR is the only machine-specific value and is never hardcoded — see the
note above its definition.
"""

import os

SRC_DIR  = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SRC_DIR)

RESULT_DIR = os.path.join(BASE_DIR, 'result')
DATA_DIR   = os.path.join(BASE_DIR, 'data')
LOG_DIR    = os.path.join(BASE_DIR, 'log')

# Result subdirectories (semantic names)
RESULT_MAPPING_DIR      = os.path.join(RESULT_DIR, 'mapping')
RESULT_LIGANDS_DIR      = os.path.join(RESULT_DIR, 'ligands')
RESULT_BINDING_DIR      = os.path.join(RESULT_DIR, 'binding_sites')
RESULT_FUNCTIONAL_DIR   = os.path.join(RESULT_DIR, 'functional')
RESULT_DISTANCES_DIR    = os.path.join(RESULT_DIR, 'distances')
RESULT_NONREDUNDANT_DIR = os.path.join(RESULT_DIR, 'nonredundant')

# Per-step log subdirectories
LOG_FUNCTIONAL_DIR = os.path.join(LOG_DIR, 'functional')

# Directory containing downloaded PDB assembly CIF files (*-assembly1.cif).
# Machine-specific, so it's not hardcoded here. Set it either via the
# STRUCTURE_CIF_DIR environment variable, or by creating a gitignored
# config_local.py next to this file with a line like:
#     CIF_DIR = '/path/to/your/PDB/assembly/cifs/'
CIF_DIR = os.environ.get('STRUCTURE_CIF_DIR', '')
if not CIF_DIR:
    try:
        from config_local import CIF_DIR  # noqa: F811
    except ImportError:
        pass  # CIF_DIR stays '' — fine for steps that don't touch CIF files

# ── nextgen_annotations ───────────────────────────────────────────────────────
NEXTGEN_BASE_URL = 'https://files-nextgen.pdbj.org/pdb_nextgen/data/entries/divided'
NEXTGEN_OUT_DIR  = os.path.join(RESULT_MAPPING_DIR, 'nextgen_annotations')
NEXTGEN_TARGET_CATEGORIES = {
    '_pdbx_sifts_unp_segments',
    '_pdbx_sifts_xref_db_segments',
}
NEXTGEN_CATEGORY_CSV = {
    '_pdbx_sifts_unp_segments':     'pdbx_sifts_unp_segments.csv',
    '_pdbx_sifts_xref_db_segments': 'pdbx_sifts_xref_db_segments.csv',
}

# ── uniprot_names ─────────────────────────────────────────────────────────────
UNIPROT_SEARCH_URL = 'https://rest.uniprot.org/uniprotkb/search'
UNIPROT_CHUNK_SIZE = 100
UNIPROT_WORKERS    = 10

# ── ligand_blacklist ──────────────────────────────────────────────────────────
CCD_PATH    = os.path.join(DATA_DIR, 'Components-rel-alt.cif')
BIOLIP_PATH = os.path.join(DATA_DIR, 'possible_artifact_molecule_list.txt')
BADLIG_ADHOC: set = {
    'UNK', 'UNL', 'DUM',
    'BMA', 'MAN', 'OXY', 'GLA', 'GLC', 'GLD',
    'BGC', 'XYS', 'XYP', 'GTX',
    'CLR', 'IOB', 'HEC', 'HEM', 'SF4', 'C8E', 'FES',
}
MOL_WT_MIN = 150.0
MOL_WT_MAX = 1000.0

# ── site_extraction ───────────────────────────────────────────────────────────
CONTACT_CUTOFF     = 4.5   # Angstroms; ligand–protein residue contact threshold
ION_LIGAND_CUTOFF  = 4.5   # Angstroms; ion proximity to organic ligand atoms
# Note: there is deliberately no STANDARD_AA set here. Whether a residue counts as
# protein is decided per structure by its label_seq (polymer membership) in
# find_resi_cif.py::_extract_atoms() and utils.py::cif2seq(), not by its CCD code —
# a name-only test would misclassify both modified in-chain residues and
# free-standing amino acids bound as ligands.
NUCLEIC_ACIDS = {
    'A', 'C', 'G', 'T', 'U',
    'DA', 'DC', 'DG', 'DT',
    'ADE', 'CYT', 'GUA', 'THY', 'URA',
    '+A', '+C', '+G', '+T', '+U',
    'I', '5IU',
}

# ── site_filter ───────────────────────────────────────────────────────────────
MIN_BINDING_RESIDUES = 3   # min residues to retain a binding site entry

# ── uniprot_features ──────────────────────────────────────────────────────────
UNIPROT_FETCH_TIMEOUT = 30   # seconds; initial per-request timeout (remote servers need >5s)
UNIPROT_RETRY_TIMEOUT = 60   # seconds; longer timeout for the retry pass

# ── functional_site_mapping ───────────────────────────────────────────────────
MIN_CHAIN_LENGTH      = 30     # minimum residues for a chain to be considered
USE_SIFTS_MATCHING    = True   # True: SIFTS arithmetic as primary; False: fuzzy-only
FUNC_SITE_WINDOW_SIZE = 10     # residues on each side of functional site in fuzzy window

# ── nonredundant_build_* / nonredundant_enrich_* ─────────────────────────────
EDIT_DIST_THRESHOLD  = 0.2   # normalised Levenshtein ≤ this → merge into one cluster
MIN_SITE_RESIDUES    = 3     # drop sites with fewer than this many AA residues in aa-tuple
RCSB_META_WORKERS    = 20
RCSB_META_TIMEOUT    = 10

# ── distances ─────────────────────────────────────────────────────────────────
DIST_MAX_WORKERS         = os.cpu_count() or 4
DIST_CHECKPOINT_INTERVAL = 100

# Site-category distance threshold, mirrored from analysis/src/config.py
# (kept in sync manually — analysis/src/config.py is the canonical definition).
#   orthosteric : min_dist <= ORTHOSTERIC_MAX_DIST
#   allosteric  : min_dist >  ORTHOSTERIC_MAX_DIST
ORTHOSTERIC_MAX_DIST = 5.0

# ── better_protein_names ──────────────────────────────────────────────────────
RESULT_PROTEIN_NAMES_DIR = os.path.join(RESULT_DIR, 'protein_names')
NONREDUNDANT_CSV_ANNOTATED = os.path.join(RESULT_NONREDUNDANT_DIR,
                                           '_nonredundant_sites_annotated.csv')
NONREDUNDANT_CSV_COMPLETE = os.path.join(RESULT_NONREDUNDANT_DIR,
                                          '_nonredundant_sites_complete.csv')
