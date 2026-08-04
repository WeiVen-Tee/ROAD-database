"""
Filter ligand binding site entries from site_res.txt.

Keeps entries with >= MIN_BINDING_RESIDUES contacting residues; drops ion-only
entries (logging them). Produces a flat dict.

Also filters site_ion_db.p (nearby-ion annotations produced by find_resi_cif.py)
to match the keys kept in site_res_db.p.

Usage:
    python filter_binding_sites.py [-v]

Inputs:  result/binding_sites/site_res.txt, site_ion_db.p, result/ligands/ions.p
Outputs: result/binding_sites/site_res_db.p, site_ion_db.p, pdb_with_ligand.p,
         dropped_ions.log
"""

import argparse
import os
import sys

from natsort import natsorted

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import RESULT_LIGANDS_DIR, RESULT_BINDING_DIR, MIN_BINDING_RESIDUES, LOG_DIR
from utils import save_pickle, load_pickle, setup_logging

log = None  # set in main()


def main(verbose: bool = False) -> None:
    """Filter site_res.txt into site_res_db.p.

    Keeps entries with at least MIN_BINDING_RESIDUES contacting residues and
    drops those whose ligand is a crystallographic ion (logged to
    dropped_ions.log). site_ion_db.p is restricted to the same keys, and
    pdb_with_ligand.p lists the PDBs that kept at least one site.
    """
    global log
    log = setup_logging(verbose, os.path.join(LOG_DIR, 'filter_binding_sites.log'))
    os.makedirs(RESULT_BINDING_DIR, exist_ok=True)

    ions = set(load_pickle(os.path.join(RESULT_LIGANDS_DIR, 'ions.p')))

    # Parse site_res.txt — the file alternates between a ligand-ID line (contains ':')
    # and a space-separated residue list on the very next line.
    raw: dict[str, list[str]] = {}
    with open(os.path.join(RESULT_BINDING_DIR, 'site_res.txt')) as fin:
        for line in fin:
            if ':' in line:
                lig = line.strip()
                try:
                    # The residue list is always on the line immediately after the ligand line.
                    res = next(fin).strip().split()
                except StopIteration:
                    log.warning(f"EOF after ligand line (truncated file?): {lig}")
                    break
                if len(res) >= MIN_BINDING_RESIDUES:
                            raw[lig] = natsorted(res)   # A.9 before A.10
    log.info(f"Parsed {len(raw)} ligand entries (>= {MIN_BINDING_RESIDUES} residues)")

    # Drop ion-only singletons — ions like ZN, MG are not small-molecule ligands;
    # keeping them would add noise to downstream distance calculations.
    site_res_db: dict[str, list[str]] = {}
    dropped = 0
    log_path = os.path.join(RESULT_BINDING_DIR, 'dropped_ions.log')
    with open(log_path, 'w') as flog:
        for lig_key, res_list in raw.items():
            # The last colon-separated field of the key is the 3-letter ligand code.
            lig_name = lig_key.split(':')[-1]
            if lig_name in ions:
                flog.write(lig_key + '\n')
                dropped += 1
            else:
                site_res_db[lig_key] = res_list

    log.info(f"Dropped {dropped} ion singletons  |  Kept {len(site_res_db)} entries")

    save_pickle(site_res_db, os.path.join(RESULT_BINDING_DIR, 'site_res_db.p'))

    # Filter site_ion_db.p to match the kept keys in site_res_db.
    ion_db_path = os.path.join(RESULT_BINDING_DIR, 'site_ion_db.p')
    raw_ion_db: dict = {}
    if os.path.exists(ion_db_path):
        raw_ion_db = load_pickle(ion_db_path)
    # For each kept ligand key, carry over its nearby-ion list (default to []).
    site_ion_db = {k: raw_ion_db.get(k, []) for k in site_res_db}
    save_pickle(site_ion_db, ion_db_path)
    log.info(f"Filtered site_ion_db: {len(site_ion_db)} entries → {ion_db_path}")

    # The PDB ID is the first colon-separated field of every key.
    pdb_with_ligand = sorted({k.split(':')[0] for k in site_res_db})
    save_pickle(pdb_with_ligand, os.path.join(RESULT_BINDING_DIR, 'pdb_with_ligand.p'))

    print(f"Binding sites: {len(site_res_db):,}  |  PDBs with ligand: {len(pdb_with_ligand):,}"
          f"  |  Dropped ions: {dropped:,}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(args.verbose)
