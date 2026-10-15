"""
run_analysis_pipeline.py — Orchestrator for the analysis pipeline stages.

Execution graph (stage: depends on)
-----------------------------------
  exclusive_categories : —
  smiles               : —
  site                 : —
  descriptors          : smiles
  aa_composition       : exclusive_categories
  recap_fragments      : smiles, exclusive_categories
  interactions         : descriptors, exclusive_categories
  build                : — (best-effort; merges whatever is on disk)
  cluster              : build

--only closes over these dependencies, so `--only cluster` also runs `build`.
`build` is the exception: its merges are all best-effort, so it never pulls the
feature stages in — it merges whatever is on disk and reports what is missing.


See README.md for the full stage/script/output table.

Usage
-----
    python run_analysis_pipeline.py [--force] [--verbose]
                           [--only smiles aa_composition ...]
                           [--cif-dir PATH] [--freesasa-bin PATH]
                           [--foldseek-bin PATH]
"""

import argparse
import os
import sys
import time
import subprocess
from concurrent.futures import ThreadPoolExecutor, Future
from typing import List, Optional

import pandas as pd

SRC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src')
sys.path.insert(0, SRC_DIR)

import config  # triggers collection path injection
import io_utils

DATA_DIR = config.DATA_DIR


# ── exclusive_categories: compute ligand_exclusive_categories.csv ───────────

def compute_exclusive_categories(force: bool) -> None:
    """
    Extraction of ligand_exclusive_categories.csv and ligand_site_counts.csv.

    Only requires the enriched sites CSV — no SMILES, no RDKit.
    Both files are produced together and consumed by analyze_ligand_descriptors.py.
    """

    out_path    = os.path.join(DATA_DIR, 'ligand_exclusive_categories.csv')
    counts_path = os.path.join(DATA_DIR, 'ligand_site_counts.csv')
    if os.path.exists(out_path) and os.path.exists(counts_path) and not force:
        print('[exclusive_categories] SKIP  ligand_exclusive_categories + ligand_site_counts — files exist')
        return

    t0 = time.time()
    print('[exclusive_categories] RUN   compute_exclusive_categories …')

    df = pd.read_csv(config.ENRICHED_CSV,
                     usecols=['ligand_code', 'min_dist', 'site_category'])
    df = io_utils.extract_ccd_code_rsplit(df)
    df = io_utils.assign_site_category_col(df, config.assign_site_category)

    # Count appearances per (ccd_code, site_category). groupby drops rows whose
    # category is unknown, so sites with no measured distance don't vote.
    counts = (df.groupby(['ccd_code', 'site_category'])
                .size()
                .unstack(fill_value=0))

    # Write n_sites for ALL ligands (before exclusive mask)
    counts_out = counts.reset_index()
    counts_out = counts_out.rename(columns={c: f'n_sites_{c}'
                                            for c in counts_out.columns if c != 'ccd_code'})
    counts_out.to_csv(counts_path, index=False)

    # Keep only ligands that appear in exactly one category
    exclusive_mask = (counts > 0).sum(axis=1) == 1
    exclusive = counts[exclusive_mask]

    if exclusive.empty:
        print('[exclusive_categories] WARNING: no exclusive ligands found — writing empty file')
        pd.DataFrame(columns=['ccd_code', 'exclusive_category']).to_csv(out_path, index=False)
        n_excl = 0
    else:
        result = exclusive.idxmax(axis=1).reset_index()
        result.columns = ['ccd_code', 'exclusive_category']
        result.to_csv(out_path, index=False)
        n_excl = len(result)

    print(f'[exclusive_categories] DONE  {n_excl} exclusive ligands + {len(counts_out)} site-count rows '
          f'({time.time()-t0:.1f}s)')


# ── Stage runner ───────────────────────────────────────────────────────────────

ALL_STAGES = ['exclusive_categories', 'smiles', 'descriptors', 'aa_composition', 'site',
              'interactions', 'recap_fragments', 'build', 'cluster']

# Hard input dependencies: stage -> stages whose output files it reads.
# `build` is deliberately absent — every one of its merges is guarded, so it runs
# against whatever exists and reports what it couldn't find (see
# build_site_ligand_features). Selecting it with --only will NOT re-run the
# feature stages; pass them explicitly, or run the whole pipeline.
_DEPS = {
    'exclusive_categories': [],
    'smiles':               [],
    'descriptors':          ['smiles'],
    'aa_composition':       ['exclusive_categories'],
    'site':                 [],
    'interactions':         ['descriptors', 'exclusive_categories'],
    'recap_fragments':      ['smiles', 'exclusive_categories'],
    'build':                [],
    'cluster':              ['build'],
}


def expand_stages(requested: List[str]) -> List[str]:
    """Close `requested` over _DEPS, returned in ALL_STAGES order.

    --only names the stages you want; their inputs are pulled in automatically,
    so `--only cluster` also runs `build`.
    """
    seen: set = set()
    stack = list(requested)
    while stack:
        stage = stack.pop()
        if stage in seen:
            continue
        seen.add(stage)
        stack.extend(_DEPS.get(stage, []))
    return [s for s in ALL_STAGES if s in seen]

_SCRIPT = {
    'smiles':          'enrich_ligand_smiles.py',
    'descriptors':     'analyze_ligand_descriptors.py',
    'aa_composition':  'analyze_aa_composition.py',
    'site':            'analyze_site.py',
    'interactions':    'analyze_interactions.py',
    'recap_fragments': 'analyze_recap_fragments.py',
    'cluster':         'cluster_protein_similarity.py',
    # 'build' is not a subprocess — handled by build_site_ligand_features() directly
}

# Stages that accept --force (build is handled separately via build_site_ligand_features)
_ACCEPTS_FORCE = {'smiles', 'descriptors', 'site', 'interactions', 'recap_fragments', 'cluster'}
# Stages that accept --verbose
_ACCEPTS_VERBOSE = {'smiles', 'descriptors', 'aa_composition', 'site',
                     'interactions', 'recap_fragments', 'cluster'}


# Every file a stage writes. All must exist for the stage to count as cached —
# checking only the first would let a run interrupted mid-stage look complete,
# and several stages write their headline file before their secondary ones.
_STAGE_OUTPUTS = {
    'smiles':          ['ligand_smiles.csv'],
    'descriptors':     ['ligand_descriptors.csv', 'plot_ligand_descriptors.csv',
                        'plot_ligand_npr.csv'],
    'aa_composition':  ['aa_mean_fractions.csv', 'site_aa_composition.csv'],
    'site':            ['sasa_results.csv', 'curvature_results.csv'],
    'interactions':    ['interaction_contacts.csv', 'interaction_contacts_shared.csv'],
    'recap_fragments': ['recap_fragments.csv', 'recap_fragment_associations.csv'],
    'build':           ['site_ligand_features.csv'],
    'cluster':         ['clusters_foldseek_all.csv', 'clusters_foldseek_allosteric.csv',
                        'largest_clusters.csv'],
}


def _stage_outputs(stage: str) -> list:
    """Absolute paths of every file `stage` writes; the first is its headline output."""
    return [os.path.join(DATA_DIR, f) for f in _STAGE_OUTPUTS[stage]]


def _primary_output(stage: str) -> str:
    return _stage_outputs(stage)[0]


def run_stage(stage: str, args: argparse.Namespace,
              extra_args: Optional[List[str]] = None) -> None:
    """Invoke a stage's script as a subprocess, respecting cache guards."""
    outputs = _stage_outputs(stage)
    missing = [p for p in outputs if not os.path.exists(p)]

    if not missing and not args.force:
        print(f'[{stage}] SKIP  {len(outputs)} output(s) exist')
        return
    if missing and len(missing) < len(outputs):
        print(f'[{stage}] RERUN incomplete — missing '
              f'{", ".join(os.path.basename(p) for p in missing)}')

    # For stages without native --force, delete the outputs to force a re-run
    if args.force and stage not in _ACCEPTS_FORCE:
        for p in outputs:
            if os.path.exists(p):
                os.remove(p)

    script = os.path.join(SRC_DIR, _SCRIPT[stage])
    cmd = [sys.executable, script]

    if args.force and stage in _ACCEPTS_FORCE:
        cmd.append('--force')
    if args.verbose and stage in _ACCEPTS_VERBOSE:
        cmd.append('--verbose')
    if extra_args:
        cmd.extend(extra_args)

    print(f'[{stage}] RUN   {" ".join(os.path.basename(c) if c.endswith(".py") else c for c in cmd)}')
    t0 = time.time()
    subprocess.run(cmd, check=True, cwd=SRC_DIR)
    print(f'[{stage}] DONE  ({time.time()-t0:.1f}s)')


# ── build stage: consolidate all per-site-ligand features ────────────────────

def build_site_ligand_features(force: bool) -> None:
    """Merge all per-site-ligand computed features into one CSV.

    Reads the enriched input CSV plus whatever stage outputs are present, and
    writes ``data/site_ligand_features.csv`` — one row per site with all computed
    columns appended. Every join is optional; missing inputs are reported, not fatal.

    Joins:
      smiles          ligand_smiles.csv          on ccd_code  (one CCD → many sites)
      descriptors     ligand_descriptors.csv     on ccd_code  (one CCD → many sites)
      aa_composition  site_aa_composition.csv    on ID
      site            sasa_results.csv + curvature_results.csv  on ID
      interactions    interaction_contacts[_shared].csv         on ID
    """
    out_path = os.path.join(DATA_DIR, 'site_ligand_features.csv')
    if os.path.exists(out_path) and not force:
        print('[build] SKIP  site_ligand_features.csv exists')
        return

    t0 = time.time()
    print('[build] RUN   build_site_ligand_features …')

    # ── Base: enriched input CSV ──────────────────────────────────────────────
    main = pd.read_csv(config.ENRICHED_CSV)
    main = io_utils.assign_site_category_col(main, config.assign_site_category)
    main = io_utils.extract_ccd_code_rsplit(main)

    # ── smiles: canonical SMILES (CCD-level, join on ccd_code) ───────────────
    smiles_path = os.path.join(DATA_DIR, 'ligand_smiles.csv')
    if os.path.exists(smiles_path):
        smiles_df = pd.read_csv(smiles_path, usecols=['ccd_code', 'canonical_smiles'])
        main = main.merge(smiles_df, on='ccd_code', how='left')
        print(f'[build]   smiles canonical SMILES merged ({len(smiles_df):,} CCD rows)')
    else:
        print('[build]   ligand_smiles.csv not found — skipping canonical_smiles')

    excl_path = os.path.join(DATA_DIR, 'ligand_exclusive_categories.csv')
    if os.path.exists(excl_path):
        excl = pd.read_csv(excl_path)
        main['is_exclusive'] = main['ccd_code'].isin(excl['ccd_code'])
    else:
        main['is_exclusive'] = False

    # ── descriptors: ligand descriptors (CCD-level, join on ccd_code) ────────
    lig_path = os.path.join(DATA_DIR, 'ligand_descriptors.csv')
    if os.path.exists(lig_path):
        lig = pd.read_csv(lig_path)
        drop_cols = [c for c in ('modal_category', 'site_category') if c in lig.columns]
        lig = lig.drop(columns=drop_cols)
        # ligand_descriptors.csv holds raw RDKit descriptors only; the per-HA
        # variants are derived here (and identically for the plot CSVs).
        lig = config.add_derived_columns(lig)
        main = main.merge(lig, on='ccd_code', how='left')
        print(f'[build]   descriptors ligand descriptors merged ({len(lig):,} CCD rows)')
    else:
        print('[build]   ligand_descriptors.csv not found — skipping')

    # ── aa_composition: per-site AA fractions (from analyze_aa_composition.py) ─
    aa_path = os.path.join(DATA_DIR, 'site_aa_composition.csv')
    if os.path.exists(aa_path):
        aa_df = pd.read_csv(aa_path)
        aa_df = aa_df.drop(columns=['source_pdb', 'ligand_code', 'site_category', 'is_exclusive'],
                            errors='ignore')
        main  = main.merge(aa_df, on='ID', how='left')
        print(f'[build]   aa_composition AA composition merged ({len(aa_df):,} site rows)')
    else:
        print('[build]   site_aa_composition.csv not found — skipping AA features')

    # ── site: SASA + curvature ────────────────────────────────────────────────
    sasa_path = os.path.join(DATA_DIR, 'sasa_results.csv')
    curv_path = os.path.join(DATA_DIR, 'curvature_results.csv')
    if os.path.exists(sasa_path) or os.path.exists(curv_path):
        frames = []
        join_cols = ['ID', 'source_pdb', 'ligand_code', 'site_category']
        if os.path.exists(sasa_path):
            frames.append(pd.read_csv(sasa_path))
        if os.path.exists(curv_path):
            frames.append(pd.read_csv(curv_path))
        site_m = frames[0]
        for other in frames[1:]:
            site_m = site_m.merge(other, on=join_cols, how='outer')
        site_m = site_m.drop(columns=['source_pdb', 'ligand_code', 'site_category'], errors='ignore')
        main = main.merge(site_m, on='ID', how='left')
        print('[build]   site SASA/curvature merged')
    else:
        print('[build]   sasa_results.csv / curvature_results.csv not found — skipping')

    # ── interactions: contact fractions (optional) ───────────────────────────
    # Merge both the exclusive-category (Analysis A) and shared-ligand (Analysis B)
    # outputs, so every analysed site carries contact fractions here. Only the
    # exclusive set is plotted and tested downstream — that filtering happens in the
    # figure and report code, which read interaction_contacts.csv alone, not by
    # withholding data from this table. Sites whose ligand is RDKit-incompatible are
    # in neither analysis and stay NaN.
    inter_path        = os.path.join(DATA_DIR, 'interaction_contacts.csv')
    inter_shared_path = os.path.join(DATA_DIR, 'interaction_contacts_shared.csv')
    inter_frames = [pd.read_csv(p) for p in (inter_path, inter_shared_path) if os.path.exists(p)]
    if inter_frames:
        # The two sets are disjoint by construction; drop_duplicates is defensive.
        inter = pd.concat(inter_frames, ignore_index=True) \
                  .drop_duplicates(subset=['ID'])
        frac_cols = [c for c in inter.columns if c.startswith('frac_')]
        main = main.merge(inter[['ID'] + frac_cols], on='ID', how='left')
        print(f'[build]   interactions contacts merged ({len(frac_cols)} frac_* columns, exclusive+shared)')

    main = main.drop(columns=[c for c in config.BUILD_EXCLUDE_COLUMNS if c in main.columns])
    left_cols = main.columns[:main.columns.get_loc('is_exclusive') + 1].tolist()
    if 'canonical_smiles' in left_cols:
        left_cols.remove('canonical_smiles')
        left_cols.insert(left_cols.index('ligand_code') + 1, 'canonical_smiles')
    # Every merge above is best-effort, so the ordered column list may name columns
    # no stage produced this run (e.g. frac_* when interactions hasn't run). Keep
    # what exists and say what's missing rather than raising KeyError.
    feature_cols = [c for c in config.SITE_LIGAND_FEATURE_COLUMN_ORDER if c in main.columns]
    absent = [c for c in config.SITE_LIGAND_FEATURE_COLUMN_ORDER if c not in main.columns]
    if absent:
        print(f'[build]   WARNING: {len(absent)} expected columns absent (upstream stage '
              f'not run?): {", ".join(absent)}')
    main = main[left_cols + feature_cols]
    main = main.round(4)
    main.to_csv(out_path, index=False)
    print(f'[build] DONE  {len(main):,} rows × {len(main.columns)} columns → {out_path}'
          f'  ({time.time()-t0:.1f}s)')


# ── Pipeline orchestration ────────────────────────────────────────────────────

def _should_run(stage: str, active_stages: List[str]) -> bool:
    return stage in active_stages


def run_parallel(active_stages: List[str], args: argparse.Namespace,
                 site_extra: List[str], cluster_extra: List[str]) -> None:
    """Run stages with maximum parallelism respecting the dependency graph
    described in the module docstring."""
    # ThreadPoolExecutor: threads are fine here since work is subprocess-bound
    with ThreadPoolExecutor(max_workers=8) as pool:

        # Phase A: exclusive_categories + smiles (both independent)
        def _run_exclusive_categories():
            if _should_run('exclusive_categories', active_stages):
                compute_exclusive_categories(args.force)

        f_step0 = pool.submit(_run_exclusive_categories)

        if _should_run('smiles', active_stages):
            f_smiles: Optional[Future] = pool.submit(run_stage, 'smiles', args)
        else:
            f_smiles = None

        # Phase B: descriptors after smiles; aa_composition after exclusive_categories
        def _run_descriptors():
            if f_smiles is not None:
                f_smiles.result()  # wait for smiles
            if _should_run('descriptors', active_stages):
                run_stage('descriptors', args)

        def _run_aa_composition():
            f_step0.result()  # wait for exclusive_categories (produces ligand_exclusive_categories.csv)
            if _should_run('aa_composition', active_stages):
                run_stage('aa_composition', args)

        f_descriptors    = pool.submit(_run_descriptors)
        f_aa_composition = pool.submit(_run_aa_composition)

        # Phase C: site reads only the enriched CSV and the CIFs, so it has no
        # in-pipeline dependency and starts immediately.
        def _run_site():
            if _should_run('site', active_stages):
                run_stage('site', args, site_extra)

        f_site = pool.submit(_run_site)

        # Phase C: interactions after descriptors (needs ligand_descriptors.csv)
        # and exclusive_categories (needs ligand_exclusive_categories.csv).
        # f_step0 must be awaited explicitly here too — submission order does
        # not imply completion order, and f_descriptors only waits on smiles.
        def _run_interactions():
            f_step0.result()
            f_descriptors.result()
            if _should_run('interactions', active_stages):
                run_stage('interactions', args)

        f_interactions = pool.submit(_run_interactions)

        # Phase C: recap_fragments after smiles (needs ligand_smiles.csv) and
        # exclusive_categories (needs ligand_exclusive_categories.csv). With
        # --only recap_fragments, f_smiles is None (smiles isn't scheduled),
        # so f_step0 is the only barrier and must be awaited explicitly.
        def _run_recap_fragments():
            f_step0.result()
            if f_smiles is not None:
                f_smiles.result()
            if _should_run('recap_fragments', active_stages):
                run_stage('recap_fragments', args)

        f_recap_fragments = pool.submit(_run_recap_fragments)

        futures = {'descriptors': f_descriptors,
                   'aa_composition': f_aa_composition, 'site': f_site,
                   'interactions': f_interactions,
                   'recap_fragments': f_recap_fragments}

        # Collect all results (re-raises any subprocess errors). build() below
        # opportunistically merges interaction_contacts*.csv if present, so it
        # must not start until interactions/recap_fragments are done.
        for name, f in futures.items():
            f.result()

    # build runs after all stages complete (no parallelism benefit)
    if _should_run('build', active_stages):
        build_site_ligand_features(args.force)

    # cluster reads site_ligand_features.csv, so it runs after build
    if _should_run('cluster', active_stages):
        run_stage('cluster', args, cluster_extra)


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Run the analysis pipeline stages.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--force', action='store_true',
                        help='Re-run all stages, ignoring all caches')
    parser.add_argument('--verbose', action='store_true',
                        help='Pass --verbose to all stages')
    parser.add_argument('--only', nargs='+', metavar='STAGE',
                        choices=ALL_STAGES,
                        help='Run only these stages (e.g. --only smiles aa_composition)')
    # site-specific pass-through args (needs CIF_DIR); --cif-dir is also
    # passed to the cluster stage's Foldseek axis
    parser.add_argument('--cif-dir', default=None, metavar='PATH',
                        help='CIF directory (passed to the site and cluster stages\' --cif-dir)')
    parser.add_argument('--freesasa-bin', default=None, metavar='PATH',
                        help='freesasa binary path (passed to the site stage --freesasa-bin)')
    # cluster-specific pass-through args
    parser.add_argument('--foldseek-bin', default=None, metavar='PATH',
                        help='foldseek binary path (passed to the cluster stage --foldseek-bin)')
    parser.add_argument('--threads', type=int, default=None, metavar='N',
                        help='Threads for foldseek (passed to the cluster stage --threads; '
                             'default: unset, foldseek uses all available cores)')
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Determine active stages (--only pulls in each stage's inputs)
    if args.only:
        active = expand_stages(args.only)
        pulled_in = [s for s in active if s not in args.only]
    else:
        active = ALL_STAGES[:]
        pulled_in = []

    # Build site-specific extra args
    site_extra: List[str] = []
    if args.cif_dir:
        site_extra.extend(['--cif-dir', args.cif_dir])
    if args.freesasa_bin:
        site_extra.extend(['--freesasa-bin', args.freesasa_bin])

    # Build cluster-specific extra args
    cluster_extra: List[str] = []
    if args.cif_dir:
        cluster_extra.extend(['--cif-dir', args.cif_dir])
    if args.foldseek_bin:
        cluster_extra.extend(['--foldseek-bin', args.foldseek_bin])
    if args.threads:
        cluster_extra.extend(['--threads', str(args.threads)])

    t_total = time.time()
    print('=' * 60)
    print('analysis pipeline')
    print(f'  stages  : {", ".join(active)}')
    if pulled_in:
        print(f'            (pulled in as dependencies: {", ".join(pulled_in)})')
    print(f'  force   : {args.force}')
    print('=' * 60)

    try:
        run_parallel(active, args, site_extra, cluster_extra)
    except subprocess.CalledProcessError as exc:
        print(f'\nPIPELINE ERROR: {exc}', file=sys.stderr)
        sys.exit(1)

    print('=' * 60)
    print(f'Pipeline complete in {time.time()-t_total:.1f}s')
    print('=' * 60)


if __name__ == '__main__':
    main()
