"""
Detect ligand-residue contacts in all *-assembly1.cif files in a directory.

For each CIF file, extracts protein residues and ligand atoms, finds contacts
within CONTACT_CUTOFF Angstroms using spatial binning + NumPy, and writes results
to a text file. Nucleic acid chains are excluded. Ion-only structures are skipped
unless a real ligand is also present.

Output format (result/binding_sites/site_res.txt):
    PDB:MODEL.CHAIN:RESNUM:RESNAME
    CHAIN.RESNUM CHAIN.RESNUM ...

Usage:
    python find_resi_cif.py <cif_dir> [-o OUTPUT] [-j WORKERS] [-v]
"""

import argparse
import glob
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import gemmi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    RESULT_LIGANDS_DIR, RESULT_BINDING_DIR, LOG_DIR,
    CONTACT_CUTOFF, ION_LIGAND_CUTOFF, NUCLEIC_ACIDS,
)
from utils import load_pickle, save_pickle, setup_logging

log = None  # set in main(); workers use their own state


def _is_nucleic_acid_chain(chain) -> bool:
    """True if the chain is a DNA, RNA, or DNA/RNA hybrid polymer.

    Returns False on any exception so edge-case chains are kept, not dropped.
    """
    try:
        polymer = chain.get_polymer()
        if not polymer:
            return False
        na_types = {gemmi.PolymerType.Dna, gemmi.PolymerType.Rna, gemmi.PolymerType.DnaRnaHybrid}
        return polymer.type in na_types
    except Exception:
        return False


def _extract_atoms(cif_path: str, bad_ligands: set, ions: set):
    """Classify a CIF's heavy atoms into (residues, ligands, ion_ligands).

    Each is a dict of 'MODEL.CHAIN:RESNUM:RESNAME' -> [[x, y, z], ...]; all three
    are empty if the file cannot be parsed.

      residues     polymer-chain residues (the protein side), including in-chain
                   modified residues and PTMs such as SEP, TPO, MSE
      ligands      non-excluded, non-ion HETATM molecules, including free-standing
                   amino acids bound as ligands
      ion_ligands  inorganic ions, kept for nearby-ion annotation

    Membership of a polymer chain (label_seq), not the CCD code, decides protein
    vs ligand. Water, nucleic-acid residues and whole nucleic-acid chains are
    excluded; alternative conformations are collapsed to the first conformer.
    """
    try:
        structure = gemmi.read_structure(cif_path)
        structure.remove_alternative_conformations() # uses only the "A" altloc
    except Exception:
        return {}, {}, {}

    residues: dict[str, list] = {}
    ligands: dict[str, list] = {}
    ion_ligands: dict[str, list] = {}

    for model_idx, model in enumerate(structure):
        model_id = model_idx + 1
        for chain in model:
            if _is_nucleic_acid_chain(chain):
                continue
            for residue in chain:
                resname = residue.name.strip()
                if resname == 'HOH' or resname in NUCLEIC_ACIDS:
                    continue

                # Build the residue key, appending insertion code when present.
                resnum = str(residue.seqid.num)
                if residue.seqid.icode.strip():
                    resnum += residue.seqid.icode.strip()
                key = f"{model_id}.{chain.name}:{resnum}:{resname}"

                # Filter on element, not atom name: a name test would also drop
                # heavy atoms named H* (mercury's atom in CCD 'HG' is named HG).
                coords = [[a.pos.x, a.pos.y, a.pos.z]
                          for a in residue if a.element.name != 'H']
                if not coords:
                    continue

                # Polymer membership, not CCD code, decides protein vs ligand: an
                # in-chain modified residue (SEP, TPO, MSE) has an integer
                # label_seq, a free-standing amino acid does not.
                in_polymer = residue.label_seq is not None
                if in_polymer:
                    residues[key] = coords
                elif resname not in bad_ligands:
                    if resname in ions:
                        ion_ligands[key] = coords
                    else:
                        ligands[key] = coords

    return residues, ligands, ion_ligands


def _find_contacts(ligands: dict, residues: dict, cutoff: float = CONTACT_CUTOFF) -> dict:
    """Residue keys within cutoff Å of each ligand, via a 3D grid of cell size cutoff.

    Ligands with no contacts are omitted.
    """
    cutoff_sq = cutoff * cutoff
    grid: dict = defaultdict(list)
    for res_key, atoms in residues.items():
        for coord in atoms:
            cell = (int(coord[0] / cutoff), int(coord[1] / cutoff), int(coord[2] / cutoff))
            grid[cell].append((res_key, coord))

    connections: dict = {}
    for lig_key, lig_atoms in ligands.items():
        contacted: set[str] = set()
        for latom in lig_atoms:
            cx, cy, cz = int(latom[0] / cutoff), int(latom[1] / cutoff), int(latom[2] / cutoff)
            nearby_atoms, nearby_keys = [], []
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for rk, rc in grid.get((cx+dx, cy+dy, cz+dz), []):
                            nearby_atoms.append(rc)
                            nearby_keys.append(rk)
            if not nearby_atoms:
                continue
            diff = np.array(nearby_atoms, dtype=np.float32) - np.array(latom, dtype=np.float32)
            contacts = np.sum(diff * diff, axis=1) < cutoff_sq
            for idx in np.where(contacts)[0]:
                contacted.add(nearby_keys[idx])
        if contacted:
            connections[lig_key] = contacted
    return connections


def _process_cif(args: tuple) -> tuple[list[str], dict[str, list[str]]]:
    """Worker: process one CIF file.

    Returns:
        (lines, ion_map) where lines is the site_res.txt output for this PDB
        and ion_map maps each full ligand key (PDB:model.chain:resnum:resname)
        to a sorted list of nearby ion strings ("{chain}.{resnum}:{resname}").
    """
    cif_path, bad_ligands, ions = args
    fname = os.path.basename(cif_path)
    m = re.match(r'^([A-Za-z0-9]{4})-assembly1\.cif$', fname, re.IGNORECASE)
    if not m:
        return [], {}
    pdb_id = m.group(1).upper()

    residues, ligands, ion_ligands = _extract_atoms(cif_path, bad_ligands, ions)
    all_ligands = {**ligands, **ion_ligands} if ligands else {}
    if not all_ligands:
        return [], {}

    connections = _find_contacts(all_ligands, residues)

    # Find which ions are within ION_LIGAND_CUTOFF of each organic ligand atom.
    # Reuse _find_contacts treating ion_ligands as the "residues" side.
    ion_connections: dict = {}
    if ligands and ion_ligands:
        ion_connections = _find_contacts(ligands, ion_ligands, cutoff=ION_LIGAND_CUTOFF)

    lines: list[str] = []
    ion_map: dict[str, list[str]] = {}
    for lig_key in sorted(connections):
        full_key = f"{pdb_id}:{lig_key}"
        lines.append(f"{full_key}\n")
        res_ids = []
        for res in connections[lig_key]:
            chain_model, resnum, _ = res.split(':')
            res_ids.append(f"{chain_model}.{resnum}")
        lines.append(' '.join(res_ids) + '\n')

        # Build ion annotation for this ligand (only for organic ligands).
        nearby = []
        for ion_key in sorted(ion_connections.get(lig_key, ())):
            ion_parts = ion_key.split(':')   # ["1.A", "100", "ZN"]
            chain = ion_parts[0].split('.')[1]
            nearby.append(f"{chain}.{ion_parts[1]}:{ion_parts[2]}")
        ion_map[full_key] = nearby

    return lines, ion_map


def main(cif_dir: str, output_file: str = None, workers: int = None,
         verbose: bool = False) -> None:
    """Process all *-assembly1.cif files in cif_dir and write contact results to site_res.txt.

    Each CIF file is processed in its own worker process (CPU-bound work).
    Results from all workers are collected into an in-memory dict keyed by PDB ID,
    then written to site_res.txt in sorted order so the output is deterministic
    regardless of the order workers finish.

    A parallel output (site_ion_db.p) records, for each organic ligand, which
    metal ions are within ION_LIGAND_CUTOFF Å.  This is used in later steps to
    annotate binding sites with nearby-ion information.

    Args:
        workers:     Number of worker processes (default: all available CPUs).
    """
    global log
    if output_file is None:
        output_file = os.path.join(RESULT_BINDING_DIR, 'site_res.txt')
    log = setup_logging(verbose, os.path.join(LOG_DIR, 'find_resi_cif.log'))
    os.makedirs(RESULT_BINDING_DIR, exist_ok=True)

    bad_ligands = load_pickle(os.path.join(RESULT_LIGANDS_DIR, 'excluded_lig.p'))
    ions        = load_pickle(os.path.join(RESULT_LIGANDS_DIR, 'ions.p'))
    bad_ligands = set(bad_ligands)
    ions        = set(ions)

    cif_files = sorted(glob.glob(os.path.join(cif_dir, '*-assembly1.cif')))
    total = len(cif_files)
    if not total:
        log.warning(f"No *-assembly1.cif files found in: {cif_dir}")
        return

    log.info(f"Processing {total} CIF files with {workers or os.cpu_count()} workers")
    print(f"Finding ligand contacts in {total} CIF files...")

    # Collect all results keyed by PDB ID for deterministic output order.
    all_results: dict[str, list[str]] = {}
    all_ion_maps: dict[str, list[str]] = {}
    completed = 0
    # Pass bad_ligands and ions as part of each work tuple — ProcessPoolExecutor
    # pickles them once per task, which is acceptable for these moderate-sized sets.
    work = [(f, bad_ligands, ions) for f in cif_files]

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_process_cif, item): item[0] for item in work}
        for future in as_completed(futures):
            cif_path = futures[future]
            try:
                lines, ion_map = future.result()
                if lines:
                    pdb_id = os.path.basename(cif_path)[:4].upper()
                    all_results[pdb_id] = lines
                    all_ion_maps.update(ion_map)
            except Exception as e:
                log.warning(f"Error processing {cif_path}: {e}")
            completed += 1
            if completed % 1000 == 0 or completed == total:
                print(f"  [{completed}/{total}] processed  ({len(all_results)} with contacts)")

    # Write site_res.txt in sorted PDB order for determinism across runs.
    with open(output_file, 'w') as fout:
        for pdb_id in sorted(all_results):
            fout.writelines(all_results[pdb_id])

    # Write ion proximity data as a pickle parallel to site_res.txt.
    ion_db_path = os.path.join(RESULT_BINDING_DIR, 'site_ion_db.p')
    save_pickle(all_ion_maps, ion_db_path)

    log.info(f"Wrote {len(all_results)} PDB entries to {output_file}")
    log.info(f"Wrote {len(all_ion_maps)} ion entries to {ion_db_path}")
    print(f"Done. {len(all_results)} PDBs with contacts → {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cif_dir", help="Directory containing *-assembly1.cif files")
    parser.add_argument("-o", "--output", default=None,
                        help="Output file (default: result/binding_sites/site_res.txt)")
    parser.add_argument("-j", "--workers", type=int, default=None,
                        help="Number of parallel workers (default: all CPUs)")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(args.cif_dir, args.output, args.workers, args.verbose)
