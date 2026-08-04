"""
Match UniProt functional site positions onto PDB structures.

Workflow:
  1. Filter to PDBs that have both a detected ligand (pdb_with_ligand) and a mapped
     functional site (unp2site); extract ±FUNC_SITE_WINDOW_SIZE residue windows from
     the canonical UniProt sequence (seq_dict) around each site position.
  2. For each PDB, parse its assembly CIF (cif2seq) and map windows to chain residues:
     primary via SIFTS label_seq arithmetic (pdbx_sifts_unp_segments.csv); fallback via
     fuzzy Levenshtein matching on chains SIFTS does not resolve.
  3. Save matched (resi_idx, 'model.chain.resnum') tuples → func_site_dict.p.

Usage:   python map_functional_sites.py [--cif-dir DIR] [-v]
Inputs:  result/mapping/pdb2unp.p, pdbx_sifts_unp_segments.csv, unp2site.p, seq_dict.p, CIF files
Outputs: result/functional/func_site_dict.p
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (RESULT_MAPPING_DIR, RESULT_BINDING_DIR, RESULT_FUNCTIONAL_DIR,
                    CIF_DIR, LOG_DIR, USE_SIFTS_MATCHING, FUNC_SITE_WINDOW_SIZE)
from utils import save_pickle, load_pickle, cif2seq, setup_logging

from rapidfuzz.distance import Levenshtein

log = None  # set in main()


# ── SIFTS lookup ──────────────────────────────────────────────────────────────

def _build_sifts_lookup(csv_path: str) -> dict:
    """(pdb_id, unp_acc) -> [(entity_id, seq_id_start, unp_start, unp_end), ...].

    seq_id_* are 1-based label_seq positions (entity numbering, including
    unmodeled residues); unp_* are 1-based canonical UniProt positions. Only
    best_mapping == 'y' rows are kept.
    """
    lookup: dict = {}
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            if row['best_mapping'] != 'y':
                continue
            key = (row['pdb_id'].upper(), row['unp_acc'])
            lookup.setdefault(key, []).append((
                row['entity_id'],
                int(row['seq_id_start']),
                int(row['unp_start']),
                int(row['unp_end']),
            ))
    return lookup


def _residue_id(chain_key: str, chain, resi_idx: int) -> str:
    """Residue ID for one chain position, e.g. '1.A.100A' (icode appended)."""
    return f'{chain_key}.{chain.auth_nums[resi_idx]}{chain.icodes[resi_idx]}'


def _sifts_map_position(pdb_id: str, unp_acc: str, uniprot_pos: int,
                         sifts_lookup: dict, pdb_chains: dict):
    """Map a UniProt position to PDB residues via SIFTS label_seq arithmetic:

        target_label_seq = seq_id_start + (uniprot_pos - unp_start)

    Every chain is searched and filtered by entity_id, so all symmetry copies of
    a homo-oligomer are covered in one pass. Returns one
    (resi_idx, 'model.chain.resnum') pair per matching chain, or [].
    """
    hits = []
    # For multi-segment constructs the first covering segment wins per chain.
    seen_chains = set()
    for entity_id, seq_start, unp_start, unp_end in sifts_lookup.get((pdb_id, unp_acc), []):
        if not (unp_start <= uniprot_pos <= unp_end):
            continue
        target_label_seq = seq_start + (uniprot_pos - unp_start)
        for chain_key, chain in pdb_chains.items():
            if chain_key in seen_chains:
                continue
            try:
                resi_idx = chain.label_seqs.index(target_label_seq)
            except ValueError:
                continue
            if chain.entity_ids[resi_idx] != entity_id:
                continue  # same label_seq in an unrelated entity — not a real hit
            hits.append((resi_idx, _residue_id(chain_key, chain, resi_idx)))
            seen_chains.add(chain_key)
    return hits


# ── Fuzzy window matching (fallback) ─────────────────────────────────────────

def _get_site_window(sequence: str, one_based_index: int, window_size: int = 10):
    """Sequence window centred on a 1-based position.

    Returns (window_str, offset_of_site_within_window), or None if the position
    is out of range. The window is clipped at sequence boundaries.
    """
    idx = one_based_index - 1
    if idx < 0 or idx >= len(sequence):
        return None
    before = sequence[max(0, idx - window_size):idx]
    after  = sequence[idx + 1:idx + 1 + window_size]
    return before + sequence[idx] + after, len(before)


def _find_best_match(protein_seq: str, query: str, max_edits: int = 4):
    """Substring position of `query` in `protein_seq` with the lowest edit distance.

    Exact-length windows are tried first, then alternately longer and shorter
    ones to absorb indels.

    Args:
        protein_seq: One-letter sequence of the PDB chain.
        query:       Window string from _get_site_window().
        max_edits:   Maximum allowable edit distance (inclusive).

    Returns:
        (start_index, matched_substring, edit_distance) of the best match,
        or None if no match within max_edits was found.
    """
    best_match = None
    best_dist  = max_edits + 1
    qlen = len(query)
    protein_len = len(protein_seq)

    # Try exact window length first, then alternate between longer and shorter.
    window_order = [qlen]
    for offset in range(1, max_edits + 1):
        window_order.append(qlen + offset)
        if qlen - offset >= 1:
            window_order.append(qlen - offset)

    for wsize in window_order:
        for i in range(max(0, protein_len - wsize + 1)):
            d = Levenshtein.distance(query, protein_seq[i:i + wsize],
                                     score_cutoff=best_dist - 1)
            if d < best_dist:
                best_dist  = d
                best_match = (i, protein_seq[i:i + wsize], d)
                if d == 0:
                    return best_match  # Exact match — no point searching further.
    return best_match


# ── Combined matching ─────────────────────────────────────────────────────────

def _match_sites(sites: list, pdb_chains: dict, pdb_id: str, sifts_lookup: dict,
                 use_sifts: bool = True) -> list:
    """Map functional site positions to PDB residues across all chains.

    sites is a list of (unp_ac, uniprot_pos, query_window, site_offset). SIFTS
    arithmetic is tried first, per chain; chains it cannot resolve fall back to
    fuzzy window matching. Pass use_sifts=False for fuzzy matching only.

    The fuzzy path is approximate: when the matched window differs in length from
    the query, the site can land one residue off.
    """
    out = []
    for unp_ac, uniprot_pos, query, site_offset in sites:
        sifts_chains_hit = set()
        if use_sifts:
            hits = _sifts_map_position(pdb_id, unp_ac, uniprot_pos, sifts_lookup, pdb_chains)
            for resi_idx, res_id in hits:
                out.append((resi_idx, res_id))
                sifts_chains_hit.add(res_id.rsplit('.', 1)[0])  # 'model.chain' prefix
        # Fuzzy: search every chain SIFTS didn't already resolve, collect all matches
        for chain_id, chain in pdb_chains.items():
            if chain_id in sifts_chains_hit:
                continue
            m = _find_best_match(chain.seq, query)
            if m:
                resi_idx = m[0] + site_offset
                if 0 <= resi_idx < len(chain.auth_nums):
                    out.append((resi_idx, _residue_id(chain_id, chain, resi_idx)))
    return out


# ── Main ──────────────────────────────────────────────────────────────────────

def main(cif_dir: str = CIF_DIR, verbose: bool = False) -> None:
    """Map UniProt functional site positions onto PDB chain sequences and save to func_site_dict.p.

    Pipeline:
      1. Load inputs: pdb2unp, unp2site (active/DNA-binding site positions), seq_dict
         (canonical UniProt sequences), pdb_with_ligand, and SIFTS segment CSV.
      2. Restrict to PDBs that have both a detected ligand AND a known functional site
         on one of their UniProt-mapped proteins.
      3. For each qualifying PDB, build a list of (unp_ac, position, query_window, offset)
         tuples by extracting ±FUNC_SITE_WINDOW_SIZE residues around each site position
         from the canonical UniProt sequence.
      4. Match those windows against the actual PDB chain sequences (from cif2seq) using
         either SIFTS arithmetic or fuzzy Levenshtein matching (controlled by
         USE_SIFTS_MATCHING in config.py).
      5. Save the result as func_site_dict.p: {pdb_id: [(resi_idx, 'model.chain.resnum'), ...]}.
    """
    global log
    os.makedirs(RESULT_FUNCTIONAL_DIR, exist_ok=True)
    log = setup_logging(verbose, os.path.join(LOG_DIR, 'functional', 'map_functional_sites.log'))

    pdb2unp         = load_pickle(os.path.join(RESULT_MAPPING_DIR, 'pdb2unp.p'))
    unp2site        = load_pickle(os.path.join(RESULT_FUNCTIONAL_DIR, 'unp2site.p'))
    seq_dict        = load_pickle(os.path.join(RESULT_FUNCTIONAL_DIR, 'seq_dict.p'))
    pdb_with_ligand = set(load_pickle(os.path.join(RESULT_BINDING_DIR, 'pdb_with_ligand.p')))

    sifts_csv = os.path.join(RESULT_MAPPING_DIR, 'nextgen_annotations', 'pdbx_sifts_unp_segments.csv')
    log.info(f"Loading SIFTS segments from {sifts_csv}")
    sifts_lookup = _build_sifts_lookup(sifts_csv)
    log.info(f"SIFTS lookup: {len(sifts_lookup)} (pdb, unp) pairs")

    # Identify PDBs with both a binding-site ligand and a known active/DNA-binding site.
    unp_with_sites = set(unp2site.uniprot_id.unique())
    # Invert pdb2unp to get a UniProt → [PDB IDs] mapping for efficient lookup.
    unp2pdb: dict = {}
    for pdb, unps in pdb2unp.items():
        for unp in unps:
            unp2pdb.setdefault(unp, []).append(pdb)

    pdbs_with_sites = {pdb for unp in unp_with_sites for pdb in unp2pdb.get(unp, [])}
    target_pdbs = list(pdb_with_ligand & pdbs_with_sites)
    log.info(f"PDBs with ligand: {len(pdb_with_ligand)}  with sites: {len(pdbs_with_sites)}"
             f"  intersection: {len(target_pdbs)}")
    print(f"Matching functional sites for {len(target_pdbs):,} PDBs...")

    # Pre-build per-PDB site window lists before the CIF-parsing loop to separate
    # pandas groupby work from file I/O.
    site_groups = unp2site.groupby('uniprot_id')['position']
    pdb2site: dict = {}
    for pdb in target_pdbs:
        windows = []
        for unp in pdb2unp.get(pdb, []):
            if unp in unp_with_sites and unp in seq_dict:
                # Deduplicate positions: range-type features (DNA binding) are stored
                # one row per residue, so overlapping features repeat positions, and
                # each duplicate would otherwise re-run the full fuzzy search.
                for pos_int in sorted({int(p) for p in site_groups.get_group(unp)}):
                    result = _get_site_window(seq_dict[unp], pos_int, FUNC_SITE_WINDOW_SIZE)
                    if result:
                        windows.append((unp, pos_int, result[0], result[1]))
        pdb2site[pdb] = windows

    # Match each PDB's site windows against its CIF chain sequences.
    func_site_dict: dict = {}
    for idx, pdb in enumerate(target_pdbs, 1):
        sites = pdb2site.get(pdb, [])
        if not sites:
            continue
        try:
            chains = cif2seq(pdb, os.path.join(cif_dir, f'{pdb}-assembly1.cif'))
            results = _match_sites(sites, chains, pdb, sifts_lookup, use_sifts=USE_SIFTS_MATCHING)
            if results:
                func_site_dict[pdb] = results
        except Exception as e:
            log.warning(f'{pdb}: {e}', exc_info=True)
        if idx % 100 == 0:
            print(f"  [{idx}/{len(target_pdbs)}] processed")

    save_pickle(func_site_dict, os.path.join(RESULT_FUNCTIONAL_DIR, 'func_site_dict.p'))
    log.info(f"func_site_dict: {len(func_site_dict)} PDBs with matched sites")
    print(f"Done. func_site_dict: {len(func_site_dict):,} PDBs → {RESULT_FUNCTIONAL_DIR}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cif-dir", default=CIF_DIR)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(args.cif_dir, args.verbose)
