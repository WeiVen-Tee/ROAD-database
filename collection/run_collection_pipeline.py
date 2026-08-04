"""
Data collection pipeline orchestrator.

Runs the full pipeline from SIFTS annotation extraction through non-redundant
binding site output. Steps are executed in dependency order; already-completed
steps (output files exist) are skipped.

Usage:
    python run_collection_pipeline.py [--cif-dir DIR] [--steps STEP ...] [-v]

Steps:
    nextgen_annotations              extract_nextgen_annotations
    uniprot_mapping                  map_pdb_nextgen_annotations
    chain_unp_map                    build_chain_unp_map
    uniprot_names                    fetch_uniprot_names
    ligand_blacklist                 build_badlig
    site_extraction                  find_resi_cif
    site_filter                      filter_binding_sites
    uniprot_features                 fetch_uniprot_features
    functional_site_mapping          map_functional_sites
    distances                        calculate_distances
    nonredundant_build_annotated     build_nonredundant_sites (annotated)
    nonredundant_build_complete      build_nonredundant_sites (complete)
    protein_names                    better_protein_names
    nonredundant_enrich_annotated    enrich_nonredundant_sites (annotated)
    nonredundant_enrich_complete     enrich_nonredundant_sites (complete)
    nonredundant_pdb_copy            copy_nonredundant_pdb
"""

import argparse
import os
import sys
import time

# Ensure src/ is on path
SRC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src')
sys.path.insert(0, SRC_DIR)

from config import (
    CIF_DIR, LOG_DIR,
    RESULT_MAPPING_DIR, RESULT_LIGANDS_DIR, RESULT_BINDING_DIR,
    RESULT_FUNCTIONAL_DIR, RESULT_DISTANCES_DIR, RESULT_NONREDUNDANT_DIR,
    NEXTGEN_OUT_DIR,
)
from utils import setup_logging
log = None  # set in main()


def _exists(*paths: str, force: bool = False) -> bool:
    """Return True only if all paths exist AND force-rerun is not requested."""
    return not force and all(os.path.exists(p) for p in paths)


def _run_step(label: str, fn, *args, **kwargs) -> None:
    """Call a step's main(), logging its timing.

    Exceptions are logged and re-raised so the pipeline stops rather than
    silently continuing past a failed step.
    """
    log.info(f"=== START {label} ===")
    t0 = time.time()
    try:
        fn(*args, **kwargs)
        log.info(f"=== DONE  {label} ({time.time()-t0:.1f}s) ===")
    except Exception as e:
        log.error(f"=== FAILED {label}: {e} ===")
        raise


def main(cif_dir: str = CIF_DIR, steps: list = None, verbose: bool = False,
         force: bool = False) -> None:
    """Run the pipeline, or the subset named in `steps` (None means all).

    Steps run in dependency order and are skipped when their outputs already
    exist, unless force=True. Each step is imported lazily so heavy optional
    dependencies are only loaded when that step actually runs.
    """
    global log
    os.makedirs(LOG_DIR, exist_ok=True)
    log = setup_logging(verbose, os.path.join(LOG_DIR, 'pipeline.log'))
    if force:
        log.info("Force re-run enabled — all step outputs will be regenerated")

    all_steps = [
        'nextgen_annotations', 'uniprot_mapping', 'chain_unp_map', 'uniprot_names',
        'ligand_blacklist', 'site_extraction', 'site_filter', 'uniprot_features',
        'functional_site_mapping', 'distances', 'nonredundant_build_annotated',
        'nonredundant_build_complete', 'protein_names',
        'nonredundant_enrich_annotated', 'nonredundant_enrich_complete',
        'nonredundant_pdb_copy',
    ]
    run = set(steps) if steps else set(all_steps)
    unknown = run - set(all_steps)
    if unknown:
        raise ValueError(f"Unknown step(s): {sorted(unknown)}. Valid steps: {all_steps}")

    def skip(*paths: str) -> bool:
        """Return True if all outputs exist and force-rerun is not requested."""
        return _exists(*paths, force=force)

    # Lazy imports (only import what is needed)
    if 'nextgen_annotations' in run:
        nextgen_csv  = os.path.join(NEXTGEN_OUT_DIR, 'pdbx_sifts_unp_segments.csv')
        nextgen_csv2 = os.path.join(NEXTGEN_OUT_DIR, 'pdbx_sifts_xref_db_segments.csv')
        if skip(nextgen_csv, nextgen_csv2):
            log.info("nextgen_annotations: Skipping (output exists)")
            print("nextgen_annotations: extract_nextgen_annotations — skipped (output exists)")
        else:
            import extract_nextgen_annotations
            # force is forwarded: this step appends to its CSVs.
            _run_step('nextgen_annotations (extract_nextgen_annotations.py)',
                      extract_nextgen_annotations.main, cif_dir,
                      verbose=verbose, force=force)

    if 'uniprot_mapping' in run:
        out  = os.path.join(RESULT_MAPPING_DIR, 'mapped_uniprot.csv')
        out2 = os.path.join(RESULT_MAPPING_DIR, 'mapped_cath.csv')
        if skip(out, out2):
            log.info("uniprot_mapping: Skipping (output exists)")
            print("uniprot_mapping: map_pdb_nextgen_annotations — skipped (output exists)")
        else:
            import map_pdb_nextgen_annotations
            _run_step('uniprot_mapping (map_pdb_nextgen_annotations.py)',
                      map_pdb_nextgen_annotations.main, cif_dir, verbose=verbose)

    if 'chain_unp_map' in run:
        out  = os.path.join(RESULT_MAPPING_DIR, 'pdb_chain_to_unp.p')
        out2 = os.path.join(RESULT_MAPPING_DIR, 'pdb2unp.p')
        out3 = os.path.join(RESULT_MAPPING_DIR, 'pdb_chain_to_unp_all.p')
        if skip(out, out2, out3):
            log.info("chain_unp_map: Skipping (output exists)")
            print("chain_unp_map: build_chain_unp_map — skipped (output exists)")
        else:
            import build_chain_unp_map
            _run_step('chain_unp_map (build_chain_unp_map.py)',
                      build_chain_unp_map.main, verbose=verbose)

    if 'uniprot_names' in run:
        out = os.path.join(RESULT_MAPPING_DIR, 'uniprot_names.p')
        if skip(out):
            log.info("uniprot_names: Skipping (output exists)")
            print("uniprot_names: fetch_uniprot_names — skipped (output exists)")
        else:
            import fetch_uniprot_names
            _run_step('uniprot_names (fetch_uniprot_names.py)',
                      fetch_uniprot_names.main, verbose=verbose)

    if 'ligand_blacklist' in run:
        out = os.path.join(RESULT_LIGANDS_DIR, 'excluded_lig.p')
        if skip(out):
            log.info("ligand_blacklist: Skipping (output exists)")
            print("ligand_blacklist: build_badlig — skipped (output exists)")
        else:
            import build_badlig
            _run_step('ligand_blacklist (build_badlig.py)', build_badlig.main, verbose=verbose)

    if 'site_extraction' in run:
        out = os.path.join(RESULT_BINDING_DIR, 'site_res.txt')
        if skip(out):
            log.info("site_extraction: Skipping (output exists)")
            print("site_extraction: find_resi_cif — skipped (output exists)")
        else:
            import find_resi_cif
            _run_step('site_extraction (find_resi_cif.py)',
                      find_resi_cif.main, cif_dir, verbose=verbose)

    if 'site_filter' in run:
        out = os.path.join(RESULT_BINDING_DIR, 'site_res_db.p')
        if skip(out):
            log.info("site_filter: Skipping (output exists)")
            print("site_filter: filter_binding_sites — skipped (output exists)")
        else:
            import filter_binding_sites
            _run_step('site_filter (filter_binding_sites.py)',
                      filter_binding_sites.main, verbose=verbose)

    if 'uniprot_features' in run:
        out = os.path.join(RESULT_FUNCTIONAL_DIR, 'unp2site.p')
        if skip(out):
            log.info("uniprot_features: Skipping (output exists)")
            print("uniprot_features: fetch_uniprot_features — skipped (output exists)")
        else:
            import fetch_uniprot_features
            # force is forwarded: this step keeps its own seq_dict.p cache.
            _run_step('uniprot_features (fetch_uniprot_features.py)',
                      fetch_uniprot_features.main, verbose=verbose, force=force)

    if 'functional_site_mapping' in run:
        out = os.path.join(RESULT_FUNCTIONAL_DIR, 'func_site_dict.p')
        if skip(out):
            log.info("functional_site_mapping: Skipping (output exists)")
            print("functional_site_mapping: map_functional_sites — skipped (output exists)")
        else:
            import map_functional_sites
            _run_step('functional_site_mapping (map_functional_sites.py)',
                      map_functional_sites.main, cif_dir, verbose=verbose)

    if 'distances' in run:
        out = os.path.join(RESULT_DISTANCES_DIR, 'distance_results.p')
        ckpt = os.path.join(RESULT_DISTANCES_DIR, 'distance_results_checkpoint.p')
        # This step uses an incremental checkpoint file so it can resume after a crash.
        # When force=True the checkpoint must be deleted, otherwise the step would
        # immediately resume from (potentially stale) prior progress.
        if force and os.path.exists(ckpt):
            os.remove(ckpt)
            log.info("distances: Removed stale checkpoint (force mode)")
        if skip(out):
            log.info("distances: Skipping (output exists)")
            print("distances: calculate_distances — skipped (output exists)")
        else:
            import calculate_distances
            _run_step('distances (calculate_distances.py)',
                      calculate_distances.main, cif_dir, verbose=verbose)

    if 'nonredundant_build_annotated' in run:
        out = os.path.join(RESULT_NONREDUNDANT_DIR, '_nonredundant_sites_annotated.p')
        if skip(out):
            log.info("nonredundant_build_annotated: Skipping (output exists)")
            print("nonredundant_build_annotated: build_nonredundant_sites (annotated) — skipped")
        else:
            import build_nonredundant_sites
            # force is forwarded: this step keeps its own pdb_cif_cache.p.
            _run_step('nonredundant_build_annotated (build_nonredundant_sites.py)',
                      build_nonredundant_sites.main, 'annotated', cif_dir,
                      verbose=verbose, force=force)

    if 'nonredundant_build_complete' in run:
        out = os.path.join(RESULT_NONREDUNDANT_DIR, '_nonredundant_sites_complete.p')
        if skip(out):
            log.info("nonredundant_build_complete: Skipping (output exists)")
            print("nonredundant_build_complete: build_nonredundant_sites (complete) — skipped")
        else:
            import build_nonredundant_sites
            _run_step('nonredundant_build_complete (build_nonredundant_sites.py)',
                      build_nonredundant_sites.main, 'complete', cif_dir,
                      verbose=verbose, force=force)

    if 'protein_names' in run:
        from config import RESULT_PROTEIN_NAMES_DIR
        # Needs both pre-enrichment site tables (their union of keys), not just
        # 'complete' — the two build modes cluster independently and can pick
        # different representative sites for the same (uniprot_ac, ligand) group.
        annotated_in = os.path.join(RESULT_NONREDUNDANT_DIR, '_nonredundant_sites_annotated.csv')
        complete_in  = os.path.join(RESULT_NONREDUNDANT_DIR, '_nonredundant_sites_complete.csv')
        out = os.path.join(RESULT_PROTEIN_NAMES_DIR, 'protein_names_pdb.csv')
        if skip(out):
            log.info("protein_names: Skipping (output exists)")
            print("protein_names: better_protein_names — skipped (output exists)")
        elif not (os.path.exists(annotated_in) and os.path.exists(complete_in)):
            log.info("protein_names: Skipping (missing inputs — run nonredundant_build_* first)")
            print("protein_names: better_protein_names — skipped (missing inputs)")
        else:
            import better_protein_names
            # force is forwarded: this step re-checks its own output file.
            _run_step('protein_names (better_protein_names.py)',
                      better_protein_names.main, cif_dir=cif_dir, verbose=verbose,
                      force=force)

    if 'nonredundant_enrich_annotated' in run:
        inp = os.path.join(RESULT_NONREDUNDANT_DIR, '_nonredundant_sites_annotated.p')
        out = os.path.join(RESULT_NONREDUNDANT_DIR, 'nonredundant_sites_annotated_enriched.p')
        if skip(out):
            log.info("nonredundant_enrich_annotated: Skipping (output exists)")
            print("nonredundant_enrich_annotated: enrich_nonredundant_sites (annotated) — skipped")
        else:
            import enrich_nonredundant_sites
            _run_step('nonredundant_enrich_annotated (enrich_nonredundant_sites.py)',
                      enrich_nonredundant_sites.main, inp, verbose=verbose)

    if 'nonredundant_enrich_complete' in run:
        inp = os.path.join(RESULT_NONREDUNDANT_DIR, '_nonredundant_sites_complete.p')
        out = os.path.join(RESULT_NONREDUNDANT_DIR, 'nonredundant_sites_complete_enriched.p')
        if skip(out):
            log.info("nonredundant_enrich_complete: Skipping (output exists)")
            print("nonredundant_enrich_complete: enrich_nonredundant_sites (complete) — skipped")
        else:
            import enrich_nonredundant_sites
            _run_step('nonredundant_enrich_complete (enrich_nonredundant_sites.py)',
                      enrich_nonredundant_sites.main, inp, verbose=verbose)

    if 'nonredundant_pdb_copy' in run:
        inp = os.path.join(RESULT_NONREDUNDANT_DIR, 'nonredundant_sites_annotated_enriched.csv')
        out = os.path.join(RESULT_NONREDUNDANT_DIR, 'PDB.tar.gz')
        out_allo = os.path.join(RESULT_NONREDUNDANT_DIR, 'PDB_allosteric.tar.gz')
        if skip(out, out_allo):
            log.info("nonredundant_pdb_copy: Skipping (output exists)")
            print("nonredundant_pdb_copy: copy_nonredundant_pdb — skipped (output exists)")
        elif not os.path.exists(inp):
            log.info("nonredundant_pdb_copy: Skipping (missing input — "
                     "run nonredundant_enrich_annotated first)")
            print("nonredundant_pdb_copy: copy_nonredundant_pdb — skipped (missing input)")
        else:
            import copy_nonredundant_pdb
            # No force to forward: each archive is written whole or not at all, so the
            # skip() check above is the only re-run decision.
            _run_step('nonredundant_pdb_copy (copy_nonredundant_pdb.py)',
                      copy_nonredundant_pdb.main, cif_dir, verbose=verbose)

    log.info("Pipeline complete.")
    print("\nPipeline complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cif-dir", default=CIF_DIR)
    parser.add_argument("--steps", nargs="+", metavar="STEP",
                        help="Run specific steps only "
                             "(e.g. --steps site_extraction site_filter distances)")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-f", "--force", action="store_true",
                        help="Force re-run all steps, ignoring existing outputs")
    args = parser.parse_args()
    main(args.cif_dir, args.steps, args.verbose, args.force)
