"""
Build chain-level and PDB-level PDB → UniProt AC mappings from SIFTS data.

Reads mapped_uniprot.csv and produces:
    pdb_chain_to_unp.p      — dict[(pdb_id, chain_label), frozenset[str]]  (best mappings only)
    pdb_chain_to_unp_all.p  — dict[(pdb_id, chain_label), frozenset[str]]  (every mapping)
    pdb2unp.p               — dict[pdb_id, set[str]]

Usage:
    python build_chain_unp_map.py [-v]
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import RESULT_MAPPING_DIR, LOG_DIR
from utils import save_pickle, setup_logging

log = None  # set in main()


def main(verbose: bool = False) -> None:
    """Build the PDB→UniProt maps from SIFTS-derived mapped_uniprot.csv.

    Writes three pickles:
      pdb_chain_to_unp.p      (pdb_id, chain) -> frozenset of accessions, SIFTS
                              best mappings only. A chain can still map to
                              several accessions (isoforms, chimeras).
      pdb_chain_to_unp_all.p  the same built from every row; used only as a
                              fallback for chains with no best mapping, so
                              restricting the primary map loses nothing.
      pdb2unp.p               pdb_id -> set of all accessions, for coarse
                              filtering and as the last-resort fallback.
    """
    global log
    log = setup_logging(verbose, os.path.join(LOG_DIR, 'mapping', 'build_chain_unp_map.log'))
    os.makedirs(RESULT_MAPPING_DIR, exist_ok=True)

    uniprot_csv = os.path.join(RESULT_MAPPING_DIR, 'mapped_uniprot.csv')
    log.info(f"Reading {uniprot_csv}")

    # Use defaultdict(set) during construction so we can add multiple UniProt ACs
    # per (pdb_id, chain) without checking for key existence each time.
    chain_best: dict[tuple, set] = defaultdict(set)
    chain_all:  dict[tuple, set] = defaultdict(set)
    pdb2unp: dict[str, set] = defaultdict(set)

    with open(uniprot_csv, newline='') as f:
        for row in csv.DictReader(f):
            pdb_id  = row['pdb_id']
            chain   = row['asym_id']   # asym_id is the author chain label (e.g. 'A').
            unp_acc = row['unp_acc']
            chain_all[(pdb_id, chain)].add(unp_acc)
            if row.get('best_mapping') == 'y':
                chain_best[(pdb_id, chain)].add(unp_acc)
            pdb2unp[pdb_id].add(unp_acc)

    # Convert mutable sets to frozensets for the chain-level maps; plain dict for PDB-level.
    chain_map     = {k: frozenset(v) for k, v in chain_best.items()}
    chain_map_all = {k: frozenset(v) for k, v in chain_all.items()}
    pdb_map       = dict(pdb2unp)

    save_pickle(chain_map,     os.path.join(RESULT_MAPPING_DIR, 'pdb_chain_to_unp.p'))
    save_pickle(chain_map_all, os.path.join(RESULT_MAPPING_DIR, 'pdb_chain_to_unp_all.p'))
    save_pickle(pdb_map,       os.path.join(RESULT_MAPPING_DIR, 'pdb2unp.p'))

    n_pdbs = len(pdb_map)
    only_non_best = len(set(chain_map_all) - set(chain_map))
    log.info(f"pdb_chain_to_unp (best): {len(chain_map):,} (pdb, chain) entries across {n_pdbs:,} PDBs")
    log.info(f"pdb_chain_to_unp_all:    {len(chain_map_all):,} entries "
             f"({only_non_best:,} chains have no best mapping)")
    log.info(f"pdb2unp: {len(pdb_map):,} PDB IDs")
    print(f"pdb_chain_to_unp: {len(chain_map):,} best entries "
          f"({len(chain_map_all):,} total)  |  pdb2unp: {len(pdb_map):,} PDB IDs")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(args.verbose)
