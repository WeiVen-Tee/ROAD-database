"""
analyze_site.py
================================================
Binding-site structural properties:

  1. Curvature  — neighbourhood density against all protein residues (proxy)
  2. SASA       — relative solvent-accessible surface area (bound vs apo)

FreeSASA is run with --radii protor and --rsa, so values come back as relative
fractions. frac_surface = fraction of site residues with rSASA_apo > 0.25;
frac_intermediate = fraction with 0 < rSASA_apo <= 0.25.

Inputs
------
  nonredundant_sites_annotated_enriched.csv
  {CIF_DIR}/{PDB_ID}-assembly1.cif

Outputs
-------
  data/sasa_results.csv
  data/curvature_results.csv

Usage
-----
    python analyze_site.py [--verbose] [--force]
"""

import argparse
import logging
import os
import string
import subprocess
import sys
from collections import defaultdict

import gemmi
import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import config  # must precede cli_utils/io_utils: injects collection/src onto sys.path
import cli_utils
import io_utils
from site_pipeline_utils import resolve_cif_path, _build_pdb_sites

logger = logging.getLogger(__name__)

# Curvature parameters
_CURVATURE_THRESHOLD = 4.5   # Å, heavy-atom distance threshold
_MIN_COORD_COVERAGE  = 0.6   # minimum fraction of residues with resolved coords
_H = gemmi.Element('H')
_D = gemmi.Element('D')



# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────

def _rename_chains(st) -> dict[str, str]:
    """Rename multi-character chain IDs to single characters, in place.

    PDB format allows only one character per chain ID, but CIF assemblies use
    names like 'AA' or '1.A'. Returns the original -> new mapping for every
    chain. Raises ValueError if there aren't enough free single-character IDs.
    """
    chain_names = [chain.name for model in st for chain in model]
    singles = list(string.ascii_uppercase + string.ascii_lowercase + string.digits)
    used = set(c for c in chain_names if len(c) == 1)
    available = [c for c in singles if c not in used]
    multi = [n for n in sorted(set(chain_names)) if len(n) > 1]
    if len(multi) > len(available):
        raise ValueError(
            f'{len(multi)} multi-character chains but only {len(available)} free '
            f'single-character IDs — cannot emit valid PDB'
        )
    mapping: dict[str, str] = {}
    idx = 0
    for name in sorted(set(chain_names)):
        if len(name) > 1:
            mapping[name] = available[idx]; idx += 1
        else:
            mapping[name] = name
    for model in st:
        for chain in model:
            chain.name = mapping.get(chain.name, chain.name)
    return mapping


# ─────────────────────────────────────────────────────────────────────────────
# Curvature (neighbourhood density)
# ─────────────────────────────────────────────────────────────────────────────

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


def _extract_protein_residues(st) -> dict[tuple[str, int, str], np.ndarray]:
    """Heavy-atom coords per polymer residue, keyed (chain, seqnum, icode).

    The residue universe curvature neighbours are counted against. Excludes
    water, nucleic acids, and anything not covalently in a polymer chain.
    """
    residues: dict[tuple[str, int, str], np.ndarray] = {}
    for chain in st[0]:
        if _is_nucleic_acid_chain(chain):
            continue
        for res in chain:
            resname = res.name.strip()
            if resname == 'HOH' or resname in config.NUCLEIC_ACIDS:
                continue
            if res.label_seq is None:
                continue
            coords = np.array([
                [atom.pos.x, atom.pos.y, atom.pos.z]
                for atom in res
                if atom.element not in (_H, _D)  # skip hydrogens and deuterium
            ])
            if len(coords) == 0:
                continue
            residues[(chain.name, res.seqid.num, res.seqid.icode.strip())] = coords
    return residues


def _build_spatial_grid(residues: dict[tuple, np.ndarray], cutoff: float) -> dict:
    """Bin heavy atoms into a 3D grid of cell size == cutoff."""
    grid: dict = defaultdict(list)
    for key, coords in residues.items():
        for coord in coords:
            cell = (int(coord[0] / cutoff), int(coord[1] / cutoff), int(coord[2] / cutoff))
            grid[cell].append((key, coord))
    return grid


def _query_contacts(query: dict[tuple, np.ndarray], grid: dict, cutoff: float) -> dict:
    """Residues within cutoff of each query residue, via a 27-cell scan of `grid`.

    Self-contacts are excluded; keys with no contacts are omitted.
    """
    cutoff_sq = cutoff * cutoff
    connections: dict = {}
    for qkey, qcoords in query.items():
        contacted: set = set()
        for qatom in qcoords:
            cx, cy, cz = int(qatom[0] / cutoff), int(qatom[1] / cutoff), int(qatom[2] / cutoff)
            nearby_atoms, nearby_keys = [], []
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for uk, uc in grid.get((cx + dx, cy + dy, cz + dz), []):
                            if uk == qkey:
                                continue
                            nearby_atoms.append(uc)
                            nearby_keys.append(uk)
            if not nearby_atoms:
                continue
            diff = np.array(nearby_atoms) - np.array(qatom)
            within = np.sum(diff * diff, axis=1) <= cutoff_sq
            for idx in np.where(within)[0]:
                contacted.add(nearby_keys[idx])
        if contacted:
            connections[qkey] = contacted
    return connections


def compute_curvature_for_structure(
        pdb_id: str, sites: list[dict], cif_dir: str,
        threshold: float) -> list[dict]:
    """Neighbourhood density for every site in one structure, as a curvature proxy.

    Counts, per site residue, how many other protein residues in the assembly
    have a heavy atom within `threshold` A. A deeply buried pocket sits in a
    densely packed region and so scores higher than a shallow surface patch.

    Returns one dict per site with keys source_pdb, ligand_code, site_category,
    n_res_total, n_res_with_coords, n_universe_residues, mean_num_neighbors;
    empty if the structure cannot be loaded.
    """
    cif_file = resolve_cif_path(pdb_id, cif_dir, None, logger, config.RCSB_ASSEMBLY_URL)
    if cif_file is None:
        return []

    try:
        st = gemmi.read_structure(cif_file)
        st.remove_alternative_conformations()
        chain_mapping = _rename_chains(st)
    except Exception as exc:
        logger.warning('%s: CIF parse failed: %s', pdb_id, exc)
        return []

    # All protein residues in the assembly, and spatial
    # grid are built once per structure and reused across every site below.
    protein_residues = _extract_protein_residues(st)
    grid = _build_spatial_grid(protein_residues, threshold)

    records = []
    for site in sites:
        residues    = site['residues']
        n_res_total = len(residues)
        site_coords: dict[tuple[str, int, str], np.ndarray] = {}
        for orig_chain, seqnum, icode in residues:
            key = (chain_mapping.get(orig_chain, orig_chain), seqnum, icode)
            coords = protein_residues.get(key)
            if coords is not None:
                site_coords[key] = coords

        n_with_coords = len(site_coords)

        if n_with_coords < 2:
            continue
        # Skip sites where too few residues have coordinates
        if n_with_coords / n_res_total < _MIN_COORD_COVERAGE:
            logger.warning('%s %s: coord coverage %.0f%% < 60%% — skipped',
                           pdb_id, site['ligand_code'],
                           100 * n_with_coords / n_res_total)
            continue

        contacts = _query_contacts(site_coords, grid, threshold)
        counts = [len(contacts.get(key, ())) for key in site_coords]

        records.append({
            'ID':                  site['ID'],
            'source_pdb':          pdb_id,
            'ligand_code':         site['ligand_code'],
            'site_category':       site['site_category'],
            'n_res_total':         n_res_total,
            'n_res_with_coords':   n_with_coords,
            'n_universe_residues': len(protein_residues),
            'mean_num_neighbors':  float(np.mean(counts)),
        })
    return records


# ─────────────────────────────────────────────────────────────────────────────
# SASA
# ─────────────────────────────────────────────────────────────────────────────

def _make_pdb_strings(cif_file: str) -> tuple[str, str, dict[str, str]]:
    """One CIF assembly -> (bound PDB, apo PDB, chain_mapping).

    Bound keeps HETATM records so ligands occlude protein surface; apo has
    ligands and waters removed. Both share one chain-name mapping so residue
    lookups stay consistent.
    """
    st_bound = gemmi.read_structure(cif_file)
    st_bound.remove_alternative_conformations()
    # Apply the same chain renaming to both copies so residue IDs match.
    chain_mapping = _rename_chains(st_bound)
    st_apo = st_bound.clone()
    # Remove all ligands and water molecules to create the apo (unliganded)
    # condition.  SASA computed on this structure reflects exposure without any
    # bound small molecule.
    st_apo.remove_ligands_and_waters()
    return st_bound.make_pdb_string(), st_apo.make_pdb_string(), chain_mapping


def _run_freesasa(pdb_str: str, include_hetatm: bool, freesasa_bin: str) -> str:
    """Run freesasa on a PDB string (via stdin) and return its RSA stdout.

    include_hetatm adds --hetatm so ligand atoms count toward occlusion — use it
    for the bound condition only.
    """
    cmd = [freesasa_bin, '--rsa', '--radii', 'protor']
    if include_hetatm:
        cmd.append('--hetatm')
    result = subprocess.run(cmd, input=pdb_str.encode(), capture_output=True, check=True)
    return result.stdout.decode()


def _parse_rsa(rsa_text: str) -> dict[tuple[str, str], float]:
    """freesasa RSA stdout -> {(chain, resnum): rSASA fraction in [0,1]}.

    Unparseable percentages are omitted.
    """
    result: dict[tuple[str, str], float] = {}
    for line in rsa_text.splitlines():
        if not line.startswith('RES'):
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        chain_field = parts[2]
        if len(chain_field) == 1:
            chain = chain_field
            resno = parts[3]
            rel_idx = 5
        else:
            chain = chain_field[0]
            resno = chain_field[1:]
            rel_idx = 4
        if rel_idx >= len(parts):
            continue
        try:
            # Cap at 100 before dividing to guard against rare freesasa rounding
            # artefacts that can produce values like 100.02%.
            rsa_pct = min(float(parts[rel_idx]), 100.0)
        except ValueError:
            continue
        # Convert from percentage (0–100) to fraction (0–1).
        result[(chain, resno)] = rsa_pct / 100.0
    return result


def _lookup_rsasa(rsa_map: dict, chain: str, resnum: str) -> float | None:
    """rSASA for one residue, keyed by (chain, residue number with any icode).

    Falls back to an insertion-coded key ('10' -> '10A') only when exactly one
    candidate exists, since freesasa and the CIF sometimes disagree about whether
    the code is part of the residue number. Residue 10 and residue 10A are
    different residues, so an ambiguous match is dropped rather than guessed.
    """
    val = rsa_map.get((chain, resnum))
    if val is not None:
        return val
    if not resnum.isdigit():
        return None
    candidates = [v for (c, rno), v in rsa_map.items()
                  if c == chain and len(rno) == len(resnum) + 1
                  and rno.startswith(resnum) and rno[-1].isalpha()]
    return candidates[0] if len(candidates) == 1 else None


def compute_sasa_for_structure(
        pdb_id: str, sites: list[dict], cif_dir: str,
        freesasa_bin: str) -> list[dict]:
    """Per-site SASA metrics for one structure, from one bound and one apo freesasa run.

    Returns one dict per site with keys source_pdb, ligand_code, site_category,
    n_res_with_sasa, n_res_total, rSASA_bound, rSASA_apo, frac_surface_apo,
    frac_intermediate_apo; empty if the structure cannot be processed.
    """
    cif_file = resolve_cif_path(pdb_id, cif_dir, None, logger, config.RCSB_ASSEMBLY_URL)
    if cif_file is None:
        return []

    try:
        bound_str, apo_str, chain_mapping = _make_pdb_strings(cif_file)
    except Exception as exc:
        logger.warning('%s: CIF→PDB conversion failed: %s', pdb_id, exc)
        return []

    try:
        # Bound condition: include HETATM so the ligand occludes protein surface.
        rsa_bound_text = _run_freesasa(bound_str, include_hetatm=True,
                                       freesasa_bin=freesasa_bin)
        # Apo condition: ligands already removed from apo_str; no HETATM needed.
        rsa_apo_text   = _run_freesasa(apo_str,   include_hetatm=False,
                                       freesasa_bin=freesasa_bin)
    except subprocess.CalledProcessError as exc:
        logger.warning('%s: freesasa failed: %s', pdb_id, exc)
        return []

    bound_map = _parse_rsa(rsa_bound_text)
    apo_map   = _parse_rsa(rsa_apo_text)

    records = []
    for site in sites:
        rsasa_bound_vals = []
        rsasa_apo_vals   = []

        for chain, seqnum, icode in site['residues']:
            mapped_chain = chain_mapping.get(chain, chain)
            resnum       = f'{seqnum}{icode}'   # freesasa keys carry the icode inline
            rsasa_bound  = _lookup_rsasa(bound_map, mapped_chain, resnum)
            rsasa_apo    = _lookup_rsasa(apo_map,   mapped_chain, resnum)
            if rsasa_bound is None or rsasa_apo is None:
                continue
            rsasa_bound_vals.append(rsasa_bound)
            rsasa_apo_vals.append(rsasa_apo)

        n = len(rsasa_bound_vals)
        n_res_total = len(site['residues'])
        # Same inclusion rule as compute_curvature_for_structure, so the two
        # metrics describe the same set of sites rather than silently averaging
        # over different fractions of one.
        if n < 2:
            logger.warning('%s %s: %d rSASA residue(s) — site skipped',
                           pdb_id, site['ligand_code'], n)
            continue
        if n / n_res_total < _MIN_COORD_COVERAGE:
            logger.warning('%s %s: rSASA coverage %.0f%% < %.0f%% — skipped',
                           pdb_id, site['ligand_code'],
                           100 * n / n_res_total, 100 * _MIN_COORD_COVERAGE)
            continue

        b = np.array(rsasa_bound_vals)
        a = np.array(rsasa_apo_vals)

        records.append({
            'ID':                  site['ID'],
            'source_pdb':          pdb_id,
            'ligand_code':         site['ligand_code'],
            'site_category':       site['site_category'],
            'n_res_with_sasa':     n,
            'n_res_total':         n_res_total,
            'rSASA_bound':         float(b.mean()),
            'rSASA_apo':           float(a.mean()),
            'frac_surface_apo':    float((a > 0.25).sum() / n),
            'frac_intermediate_apo': float(((a > 0) & (a <= 0.25)).sum() / n),
        })
    return records


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    """Compute neighbourhood-density curvature and bound/apo SASA for every site.

    Results are cached in data/curvature_results.csv and data/sasa_results.csv;
    --force recomputes. Statistical comparison between site categories happens in
    downstream, not here.
    """
    parser = argparse.ArgumentParser(
        description='Combined site analysis: size, curvature, SASA.')
    parser.add_argument('--cif-dir', default=config.CIF_DIR,
                        help=f'CIF directory (default: {config.CIF_DIR})')
    parser.add_argument('--freesasa-bin', default='freesasa',
                        help='Path to freesasa executable (default: freesasa on PATH)')
    parser.add_argument('--threshold', type=float, default=_CURVATURE_THRESHOLD,
                        help=f'Curvature distance threshold in Å (default {_CURVATURE_THRESHOLD})')
    cli_utils.add_common_args(parser, force=True)
    args = parser.parse_args()

    t0 = cli_utils.init_run('analyze_site', args, config.LOG_DIR)

    # ── Load enriched CSV ──────────────────────────────────────────────────────
    print('\nLoading enriched CSV …')
    enrich_df = io_utils.load_enriched_csv(config.ENRICHED_CSV, dropna_subset=['min_dist', 'binding_site'])
    enrich_df = io_utils.assign_site_category_col(enrich_df, config.assign_site_category)
    print(f'  {len(enrich_df):,} binding sites, '
          f'{enrich_df["source_pdb"].nunique():,} structures.')

    # ── 1. Curvature ───────────────────────────────────────────────────────────
    print('\n[1/2] Computing curvature …')
    curv_path = os.path.join(config.DATA_DIR, 'curvature_results.csv')

    if io_utils.output_is_cached(curv_path, args.force):
        print('  curvature_results.csv exists — loading cache (--force to recompute).')
        curv_df = pd.read_csv(curv_path)
        print(f'  {len(curv_df):,} cached site rows.')
    else:
        all_pdbs = sorted(enrich_df['source_pdb'].str.upper().unique())
        pdb_sites = _build_pdb_sites(enrich_df, all_pdbs)
        all_records = []
        n_ok = n_fail = 0
        print(f'  Processing {len(all_pdbs):,} structures (threshold={args.threshold} Å) …')
        for i, pdb_id in enumerate(all_pdbs):
            if (i + 1) % 50 == 0 or i == 0:
                print(f'    [{i+1}/{len(all_pdbs)}] {pdb_id}', flush=True)
            records = compute_curvature_for_structure(
                pdb_id, pdb_sites[pdb_id], args.cif_dir, args.threshold)
            if records:
                all_records.extend(records); n_ok += 1
            else:
                n_fail += 1

        print(f'  OK: {n_ok:,}  failed: {n_fail:,}  site rows: {len(all_records):,}')
        if all_records:
            curv_df = pd.DataFrame(all_records)
            curv_df.to_csv(curv_path, index=False)
            print(f'  Saved → {curv_path}')
        else:
            curv_df = pd.DataFrame()

    # ── 2. SASA ────────────────────────────────────────────────────────────────
    print('\n[2/2] Computing SASA …')
    sasa_path = os.path.join(config.DATA_DIR, 'sasa_results.csv')

    if io_utils.output_is_cached(sasa_path, args.force):
        sasa_df = pd.read_csv(sasa_path)
        print(f'  Loaded {len(sasa_df):,} cached site rows.')
    else:
        all_pdbs = sorted(enrich_df['source_pdb'].str.upper().unique())
        pdb_sites = _build_pdb_sites(enrich_df, all_pdbs)

        print(f'  Running freesasa on {len(all_pdbs):,} structures …')
        all_sasa_records = []
        n_ok = n_fail = 0
        for i, pdb_id in enumerate(all_pdbs):
            if (i + 1) % 50 == 0 or i == 0:
                print(f'    [{i+1}/{len(all_pdbs)}] {pdb_id} …', flush=True)
            records = compute_sasa_for_structure(
                pdb_id, pdb_sites[pdb_id], args.cif_dir, args.freesasa_bin)
            if records:
                all_sasa_records.extend(records); n_ok += 1
            else:
                n_fail += 1

        print(f'  OK: {n_ok:,}  failed: {n_fail:,}  site rows: {len(all_sasa_records):,}')
        if all_sasa_records:
            sasa_df = pd.DataFrame(all_sasa_records)
            sasa_df.to_csv(sasa_path, index=False)
            print(f'  Saved → {sasa_path}')
        else:
            sasa_df = pd.DataFrame()

    cli_utils.finish_run(t0)


if __name__ == '__main__':
    main()
