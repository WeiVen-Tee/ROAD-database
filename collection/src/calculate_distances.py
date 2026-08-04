"""
Calculate distances between functional site residues and ligand binding site residues.

For each PDB with both a functional site annotation and a binding site, computes
min, median, and avg heavy-atom pairwise distances (Å), rounded to 2 decimal places.
Parallel across structures; checkpoints every N structures.

Usage:
    python calculate_distances.py [--cif-dir DIR] [-v]

Inputs:  result/ligands/ions.p, result/binding_sites/site_res_db.p,
         result/functional/func_site_dict.p, CIF files
Outputs: result/distances/distance_results.p
         Values: dict[pdb_id, dict[lig_key, (min_dist, median_dist, avg_dist)]]
"""

import argparse
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, FrozenSet, List, Tuple

import numpy as np
import gemmi
from scipy.spatial.distance import cdist

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    RESULT_LIGANDS_DIR, RESULT_BINDING_DIR, RESULT_FUNCTIONAL_DIR,
    RESULT_DISTANCES_DIR, CIF_DIR, LOG_DIR,
    DIST_MAX_WORKERS, DIST_CHECKPOINT_INTERVAL,
)
from utils import save_pickle, load_pickle, setup_logging

log = None  # set in main()

# Worker-process globals (loaded once via initializer, not re-pickled per task)
_worker_binding_site_dict: Dict = {}
_worker_ions: FrozenSet[str] = frozenset()


def _worker_init(binding_site_dict: Dict, ions: FrozenSet[str]) -> None:
    """Store the large shared inputs as worker globals, transferred once per process."""
    global _worker_binding_site_dict, _worker_ions
    _worker_binding_site_dict = binding_site_dict
    _worker_ions = ions


def _extract_coords(structure_path: str, residue_ids: List[str]) -> Dict[str, np.ndarray]:
    """Read a CIF structure file and return 3D atomic coordinates for a set of residues.

    Only heavy (non-hydrogen) atoms are collected

    Args:
        residue_ids:    List of residue IDs in 'MODEL.CHAIN.RESNUM[ICODE]' format.

    Returns:
        Dict mapping each residue ID to a float32 numpy array of shape (N_atoms, 3).
        Residues with no atoms found are absent from the returned dict.
    """
    if not os.path.exists(structure_path):
        return {}

    # Pre-build a lookup from (model_idx, chain, resnum, icode) → [residue_id, ...]
    # so that traversing the structure is O(atoms) rather than O(atoms × residues).
    lookup: Dict = {}
    for res_id in residue_ids:
        parts = res_id.split('.')
        if len(parts) != 3:
            continue
        # gemmi model indices are 0-based; our IDs use 1-based model numbers.
        model_num = int(parts[0]) - 1
        chain_id  = parts[1]
        res_str   = parts[2]
        try:
            # Simple case: pure integer residue number, no insertion code.
            res_num, icode = int(res_str), ''
        except ValueError:
            # Insertion code: find where the digits end (e.g. '10A' → num=10, icode='A').
            for i, ch in enumerate(res_str):
                if not ch.isdigit() and ch != '-':
                    res_num, icode = int(res_str[:i]), res_str[i:]
                    break
            else:
                continue
        lookup.setdefault((model_num, chain_id, res_num, icode), []).append(res_id)

    if not lookup:
        return {}

    try:
        structure = gemmi.read_structure(structure_path)
    except Exception as e:
        print(f"  [WARN] Could not read structure {structure_path}: {e}", flush=True)
        return {}

    coords: Dict[str, list] = {rid: [] for rid in residue_ids}
    for model_idx, model in enumerate(structure):
        for chain in model:
            for residue in chain:
                key = (model_idx, chain.name, residue.seqid.num, residue.seqid.icode.strip())
                if key not in lookup:
                    continue
                for atom in residue:
                    if atom.element.name == 'H':
                        continue  # Skip hydrogen atoms
                    pos = atom.pos
                    for rid in lookup[key]:
                        coords[rid].append([pos.x, pos.y, pos.z])

    # Only include residues that actually had atoms in the file.
    return {rid: np.array(c, dtype=np.float32) for rid, c in coords.items() if c}


def _is_ion(lig_id: str, ions: FrozenSet[str]) -> bool:
    """True if this ligand key's 3-letter code is a crystallographic ion.

    The code is the last colon-separated field of 'PDB:MODEL.CHAIN:RESNUM:RESNAME'.
    """
    return lig_id.split(':')[-1] in ions


def _process_pdb(pdb: str, cif_dir: str, functional_residues: List[str],
                 ligand_ids: List[str]) -> Tuple[str, Dict]:
    """Functional-site to binding-site distances for every ligand in one PDB.

    Runs in a worker process, reading the shared inputs from the globals set by
    _worker_init(). All needed coordinates are extracted in a single pass over
    the CIF.

    Returns (pdb, {lig_id: (min, median, avg) in A}); a ligand whose residues
    have no extractable coordinates maps to (None, None, None).
    """
    binding_site_dict = _worker_binding_site_dict
    ions = _worker_ions
    cif_path = os.path.join(cif_dir, f"{pdb}-assembly1.cif")
    if not os.path.exists(cif_path) or not functional_residues:
        return pdb, {}

    # Collect all residue IDs we will need in a single set to parse the CIF once.
    all_residues = set(functional_residues)
    for lig_id in ligand_ids:
        if _is_ion(lig_id, ions):
            continue
        all_residues.update(binding_site_dict[lig_id])

    all_coords = _extract_coords(cif_path, list(all_residues))

    func_coords = {r: all_coords[r] for r in functional_residues if r in all_coords}
    if not func_coords:
        print(f"  [WARN] {pdb}: no functional-residue coords found in CIF", flush=True)
        return pdb, {}

    func_atoms = np.vstack(list(func_coords.values()))
    results: Dict = {}

    for lig_id in ligand_ids:
        if _is_ion(lig_id, ions):
            continue
        binding_res = binding_site_dict[lig_id]
        bind_coords = {r: all_coords[r] for r in binding_res if r in all_coords}
        if not bind_coords:
            print(f"  [WARN] {pdb}/{lig_id}: no binding-residue coords found", flush=True)
            results[lig_id] = (None, None, None)
            continue
        bind_atoms = np.vstack(list(bind_coords.values()))
        dists = cdist(func_atoms, bind_atoms, metric='euclidean')
        results[lig_id] = (
            round(float(np.min(dists)),    2),
            round(float(np.median(dists)), 2),
            round(float(np.mean(dists)),   2),
        )

    return pdb, results


def _run_batch(cif_dir: str, func_site_dict: Dict, binding_site_dict: Dict,
               ions: FrozenSet[str], checkpoint_file: str,
               max_workers: int = DIST_MAX_WORKERS,
               checkpoint_interval: int = DIST_CHECKPOINT_INTERVAL) -> Dict:
    """Run _process_pdb across a process pool, checkpointing as it goes.

    Only PDBs present in both func_site_dict and binding_site_dict are processed.
    Progress is written to checkpoint_file every checkpoint_interval structures
    so an interrupted run can resume; pass '' to disable.

    Returns {pdb_id: {lig_key: (min_dist, median_dist, avg_dist)}}.
    """
    # Group all ligand keys by their PDB ID for efficient dispatch.
    ligands_by_pdb: Dict[str, List[str]] = {}
    for lig_id in binding_site_dict:
        ligands_by_pdb.setdefault(lig_id.split(':')[0], []).append(lig_id)

    # Only process PDBs that have both functional-site and binding-site annotations.
    common_pdbs = sorted(set(func_site_dict) & set(ligands_by_pdb))
    total = len(common_pdbs)

    results: Dict = {}
    processed: set = set()
    if checkpoint_file and os.path.exists(checkpoint_file):
        try:
            with open(checkpoint_file, 'rb') as f:
                ckpt = pickle.load(f)
            results   = ckpt.get('results', {})
            processed = ckpt.get('processed', set())
            print(f"Resumed from checkpoint: {len(processed):,} done, "
                  f"{total - len(processed):,} remaining")
        except Exception as e:
            print(f"Could not load checkpoint: {e}. Starting from scratch.")

    remaining = [p for p in common_pdbs if p not in processed]
    if not remaining:
        print("All structures already processed.")
        return results

    # Each work item is a tuple of arguments passed to _process_pdb.
    # The large shared dicts (binding_site_dict, ions) are NOT included here —
    # they were passed to _worker_init and live in each worker's global namespace.
    work = [(p, cif_dir, [r for _, r in func_site_dict[p]], ligands_by_pdb[p])
            for p in remaining]

    print(f"Processing {total:,} structures ({max_workers} workers)...")
    initial_done = len(processed)
    completed    = initial_done
    last_ckpt = completed
    start = time.time()

    with ProcessPoolExecutor(max_workers=max_workers,
                             initializer=_worker_init,
                             initargs=(binding_site_dict, ions)) as executor:
        futures = {executor.submit(_process_pdb, *item): item[0] for item in work}
        for future in as_completed(futures):
            try:
                pdb, pdb_results = future.result()
                if pdb_results:
                    results[pdb] = pdb_results
                processed.add(pdb)
            except Exception as e:
                print(f"  Error {futures[future]}: {e}")
            completed += 1
            if completed % 100 == 0 or completed == total:
                elapsed = time.time() - start
                rate = (completed - initial_done) / elapsed if elapsed > 0 else 0
                eta  = (total - completed) / rate if rate > 0 else 0
                print(f"  [{completed:>6,}/{total:,}] {100*completed/total:.1f}% | "
                      f"{rate:.1f}/s | ETA {eta/60:.1f} min")
            # Write checkpoint when enough new results have accumulated.
            if checkpoint_file and completed - last_ckpt >= checkpoint_interval:
                with open(checkpoint_file, 'wb') as f:
                    pickle.dump({'results': results, 'processed': processed,
                                 'timestamp': time.time()}, f)
                last_ckpt = completed

    # Final checkpoint marks the run as fully complete.
    if checkpoint_file:
        with open(checkpoint_file, 'wb') as f:
            pickle.dump({'results': results, 'processed': processed,
                         'timestamp': time.time(), 'completed': True}, f)

    elapsed = time.time() - start
    print(f"Completed: {len(results):,} structures in {elapsed/60:.1f} min")
    return results


def main(cif_dir: str = CIF_DIR, verbose: bool = False) -> None:
    """Orchestrate the parallel distance calculation and save results to disk.

    Loads the three required inputs (ions, binding sites, functional sites), runs
    the parallel batch, persists the result, and writes a diagnostic file listing
    binding sites that were not processed (e.g. because the PDB had no functional
    site annotation or the CIF file was missing).
    """
    global log
    os.makedirs(RESULT_DISTANCES_DIR, exist_ok=True)
    log = setup_logging(verbose, os.path.join(LOG_DIR, 'calculate_distances.log'))

    # Use frozenset for ions so membership tests in worker processes are O(1).
    ions           = frozenset(load_pickle(os.path.join(RESULT_LIGANDS_DIR, 'ions.p')))
    site_res_db    = load_pickle(os.path.join(RESULT_BINDING_DIR, 'site_res_db.p'))
    func_site_dict = load_pickle(os.path.join(RESULT_FUNCTIONAL_DIR, 'func_site_dict.p'))
    log.info(f"ions={len(ions)}  binding_sites={len(site_res_db)}  func_sites={len(func_site_dict)}")

    distance_results = _run_batch(
        cif_dir, func_site_dict, site_res_db, ions,
        checkpoint_file=os.path.join(RESULT_DISTANCES_DIR, 'distance_results_checkpoint.p'),
    )

    out_path = os.path.join(RESULT_DISTANCES_DIR, 'distance_results.p')
    save_pickle(distance_results, out_path)
    log.info(f"Saved distance_results.p: {len(distance_results):,} PDBs")

    total_entries = sum(len(v) for v in distance_results.values())
    # Count ligand entries where min_dist is not None (coordinates were found in the CIF).
    # Example scenario when min_dist is None -- binding residues not present in the assembly CIF
    with_coords   = sum(1 for v in distance_results.values()
                        for (mn, _med, _avg) in v.values() if mn is not None)
    print(f"PDBs: {len(distance_results):,}  entries: {total_entries:,}  with coords: {with_coords:,}")

    # Write a diagnostic file listing all binding-site ligands that have no distance result,
    # useful for debugging missing CIF files or absent functional-site annotations.
    missing_path = os.path.join(RESULT_DISTANCES_DIR, 'missing_from_results.txt')
    with open(missing_path, 'w') as fh:
        for lig_id in sorted(site_res_db):
            if lig_id.split(':')[0] not in distance_results:
                fh.write(f'{lig_id}\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cif-dir", default=CIF_DIR)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(args.cif_dir, args.verbose)
