"""
analyze_recap_fragments.py  –  Fragment Analysis
=======================================================
Fragments exclusive-category ligands using RDKit RECAP (all tree levels) and
reports descriptive per-fragment occurrence counts in orthosteric vs
allosteric ligands. RECAP fragmentation runs in parallel across CPU cores.

Fragment size filters:
  - Skip root node (= input molecule itself)
  - Skip trivial/small:  HAC <= --min-hac  (default 8)
  - Skip too-large:      HAC >  --max-hac  (default 20)

Inputs
------
  data/ligand_exclusive_categories.csv
  data/ligand_smiles.csv
  nonredundant_sites_annotated_enriched.csv

Cache
-----
  data/recap_frag_cache_hac{min}-{max}.pkl   RECAP inverted index (frag_smi → {ccd_codes})
                              + HAC map, keyed by the --min-hac/--max-hac thresholds used
                              to build it. Loaded on subsequent runs with matching thresholds
                              to skip the expensive RECAP step. Rebuilt only when --force is
                              supplied.

Outputs
-------
  data/recap_fragments.csv               per-fragment descriptive counts
  data/recap_fragment_associations.csv   ccd_code ↔ fragment rows

Usage
-----
    python analyze_recap_fragments.py [--verbose] [--force]
                                      [--min-hac N] [--max-hac N] [--workers N]
"""

import argparse
import logging
import multiprocessing as mp
import os
import sys
import threading
from collections import defaultdict

import pandas as pd
from rdkit import Chem, rdBase
from rdkit.Chem.Recap import RecapDecompose

rdBase.DisableLog('rdApp.*')

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import config  # must precede utils/cli_utils/io_utils: injects collection/src onto sys.path
from utils import load_pickle, save_pickle
import cli_utils
import io_utils

logger = logging.getLogger(__name__)

CATEGORIES     = config.CATEGORIES
FRAG_CAP         = 30   # max fragments per molecule
_RECAP_TIMEOUT_S = 15   # per-molecule timeout via daemon thread join

_DUMMY_SMARTS  = None     # lazy-loaded

# Plot data: recap_fragments.csv saved by this module contains fragment SMILES
# and descriptive counts, read directly to generate 2D structure grids.


def _get_dummy_query():
    """Cached query mol matching dummy atoms ([#0]); parsed once on first use."""
    global _DUMMY_SMARTS
    if _DUMMY_SMARTS is None:
        _DUMMY_SMARTS = Chem.MolFromSmarts('[#0]')
    return _DUMMY_SMARTS


# ─────────────────────────────────────────────────────────────────────────────
# Fragmentation
# ─────────────────────────────────────────────────────────────────────────────

_frag_min_hac: int = 8
_frag_max_hac: int = 20


def _strip_dummies(mol):
    """Strip RECAP's [*] attachment-point atoms and sanitize. None if that fails."""
    try:
        result = Chem.DeleteSubstructs(mol, _get_dummy_query())
        if result is None:
            return None
        Chem.SanitizeMol(result)
        return result
    except Exception:
        return None


def _all_nodes(node):
    """Depth-first walk of a RECAP tree, yielding every node — not just leaves."""
    yield node
    for child in node.children.values():
        yield from _all_nodes(child)


def _recap_worker(args):
    """RECAP-fragment one molecule; returns (ccd_code, {(frag_smi, hac): None}).

    Walks the whole tree, strips dummies, filters by HAC, dedups by canonical
    SMILES, and caps at FRAG_CAP fragments. Sentinel keys report outcomes across
    the process boundary without raising: ('__FAILED__', -2) for a parse or
    fragmentation failure, ('__CAPPED__', -1) when the cap was hit.

    A dict rather than a set, so key order is insertion order rather than
    Python's randomized string hashing — run_recap_batch needs that for
    reproducible fragment_index tie-breaking.
    """
    ccd_code, smiles = args
    result: dict = {}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result[('__FAILED__', -2)] = None
        return ccd_code, result

    input_smi = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)

    try:
        tree = RecapDecompose(mol)
    except Exception:
        result[('__FAILED__', -2)] = None
        return ccd_code, result

    capped = False
    for node in _all_nodes(tree):
        frag_mol = node.mol
        if frag_mol is None:
            continue

        cleaned = _strip_dummies(frag_mol)
        if cleaned is None or cleaned.GetNumAtoms() == 0:
            continue

        hac = cleaned.GetNumHeavyAtoms()

        try:
            frag_smi = Chem.MolToSmiles(cleaned, canonical=True, isomericSmiles=False)
        except Exception:
            continue
        if frag_smi == input_smi:
            # Skip the root node (= original molecule before any bond cleavage)
            continue

        if hac <= _frag_min_hac or hac > _frag_max_hac:
            continue

        result[(frag_smi, hac)] = None

        if len(result) >= FRAG_CAP:
            capped = True
            break

    if capped:
        result[('__CAPPED__', -1)] = None

    return ccd_code, result


def _fragment_molecule(item) -> tuple:
    """_recap_worker with a hard timeout, via a daemon thread.

    RDKit's RECAP C extension can't be interrupted once running, so the only way
    to bound it is to abandon the thread. Daemon so a hung one can't keep the
    process alive. Returns the ('__TIMEOUT__', -3) sentinel on expiry.
    """
    result = [None]

    def _target():
        """Run _recap_worker in a daemon thread so it can be killed on timeout."""
        result[0] = _recap_worker(item)

    t = threading.Thread(target=_target, daemon=True)
    t.start()
    t.join(timeout=_RECAP_TIMEOUT_S)
    if t.is_alive():
        return item[0], {('__TIMEOUT__', -3): None}
    return result[0]


def _recap_worker_init(min_hac: int, max_hac: int) -> None:
    """Pool initializer for the HAC thresholds.

    Needed because 'spawn' re-imports this module in each worker, so a global
    set in the parent never reaches them.
    """
    global _frag_min_hac, _frag_max_hac
    _frag_min_hac = min_hac
    _frag_max_hac = max_hac
    rdBase.DisableLog('rdApp.*')


def run_recap_batch(items: list, min_hac: int, max_hac: int,
                    n_workers: int | None = None) -> tuple[dict, dict, dict]:
    """RECAP-fragment `items` in parallel, serially for small inputs.

    Uses ordered imap so results arrive in input order — fragment_index ties are
    broken by insertion order into the inverted index and must stay reproducible.

    Returns (inverted: frag_smi -> {ccd_codes}, frag_hac: frag_smi -> HAC,
    stats: n_capped/n_failed/n_timeout/n_total).
    """
    n_total = len(items)
    n_workers = n_workers or max(1, os.cpu_count() or 1)
    use_parallel = n_total >= 200 and n_workers > 1

    inverted: dict[str, set] = defaultdict(set)
    frag_hac: dict[str, int] = {}
    n_capped = n_failed = n_timeout = 0
    PROGRESS_EVERY = max(1, n_total // 20)

    def _consume(ccd_code, frag_set, n_processed):
        nonlocal n_capped, n_failed, n_timeout
        if n_processed % PROGRESS_EVERY == 0 or n_processed == n_total:
            print(f'\r    {n_processed:,}/{n_total:,} molecules processed ...',
                  end='', flush=True)
        if ('__TIMEOUT__', -3) in frag_set:
            n_timeout += 1
            logger.warning('RECAP timeout (>%ds) for CCD %s — skipped',
                           _RECAP_TIMEOUT_S, ccd_code)
            return
        if ('__FAILED__', -2) in frag_set:
            n_failed += 1
            return
        if ('__CAPPED__', -1) in frag_set:
            frag_set.pop(('__CAPPED__', -1), None)
            n_capped += 1
            logger.warning('Fragment cap hit for CCD %s', ccd_code)
        for frag_smi, hac in frag_set:
            # sys.intern() stores only one copy of each repeated string in memory.
            # Fragment SMILES strings repeat heavily across thousands of molecules;
            # interning them reduces memory usage for the inverted index.
            inverted[sys.intern(frag_smi)].add(ccd_code)
            frag_hac[frag_smi] = hac

    if use_parallel:
        with mp.Pool(n_workers, initializer=_recap_worker_init,
                    initargs=(min_hac, max_hac), maxtasksperchild=200) as pool:
            chunksize = max(1, n_total // (n_workers * 10))
            for i, (ccd_code, frag_set) in enumerate(
                    pool.imap(_fragment_molecule, items, chunksize=chunksize), 1):
                _consume(ccd_code, frag_set, i)
    else:
        global _frag_min_hac, _frag_max_hac
        _frag_min_hac, _frag_max_hac = min_hac, max_hac
        for i, item in enumerate(items, 1):
            ccd_code, frag_set = _fragment_molecule(item)
            _consume(ccd_code, frag_set, i)
    print()

    stats = {'n_capped': n_capped, 'n_failed': n_failed,
             'n_timeout': n_timeout, 'n_total': n_total}
    return dict(inverted), frag_hac, stats


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    """Fragment exclusive-category ligands and write per-fragment counts."""
    parser = argparse.ArgumentParser(description='RECAP fragment analysis')
    cli_utils.add_common_args(parser, force=True)
    parser.add_argument('--min-hac',    type=int, default=8,   dest='min_hac',
                        help='Skip fragments with HAC <= this (default 8)')
    parser.add_argument('--max-hac',    type=int, default=20,  dest='max_hac',
                        help='Skip fragments with HAC > this (default 20)')
    parser.add_argument('--workers',    type=int, default=None, dest='workers',
                        help='Number of worker processes for RECAP fragmentation '
                             '(default: os.cpu_count())')
    args = parser.parse_args()

    out_frags  = os.path.join(config.DATA_DIR, 'recap_fragments.csv')
    out_assoc  = os.path.join(config.DATA_DIR, 'recap_fragment_associations.csv')
    cache_path = os.path.join(config.DATA_DIR,
                              f'recap_frag_cache_hac{args.min_hac}-{args.max_hac}.pkl')
    t0 = cli_utils.init_run('analyze_recap_fragments', args, config.LOG_DIR)
    print(f'  min_hac={args.min_hac}  max_hac={args.max_hac}  frag_cap={FRAG_CAP}')

    # ── [1/4] Load inputs ────────────────────────────────────────────────────
    print('\n[1/4] Loading inputs ...')
    excl_path   = os.path.join(config.DATA_DIR, 'ligand_exclusive_categories.csv')
    smiles_path = os.path.join(config.DATA_DIR, 'ligand_smiles.csv')
    excl_df   = pd.read_csv(excl_path)
    smiles_df = pd.read_csv(smiles_path)
    enriched  = pd.read_csv(config.ENRICHED_CSV, low_memory=False)

    enriched = io_utils.extract_ccd_code_rsplit(enriched)

    # All unique UniProt IDs per CCD (across all rows, handles multiple sites per protein)
    ccd_to_uniprots: dict[str, set] = (
        enriched.dropna(subset=['uniprot_ac'])
        .groupby('ccd_code')['uniprot_ac']
        .apply(set)
        .to_dict()
    )

    merged = excl_df.merge(smiles_df[['ccd_code', 'canonical_smiles']], on='ccd_code', how='inner')

    n_before = len(merged)
    merged = merged.dropna(subset=['canonical_smiles'])
    n_dropped = n_before - len(merged)
    if n_dropped:
        print(f'  Dropped {n_dropped:,} CCD codes with missing canonical_smiles')
        logger.info('Dropped %d CCD codes with missing canonical_smiles', n_dropped)

    cat_map     = merged.set_index('ccd_code')['exclusive_category'].to_dict()
    total_ortho = int((merged['exclusive_category'] == 'orthosteric').sum())
    total_allo  = int((merged['exclusive_category'] == 'allosteric').sum())
    print(f'  {len(merged):,} exclusive CCD codes '
          f'({total_ortho:,} ortho, {total_allo:,} allo)')
    logger.info('Loaded %d exclusive CCD codes (%d ortho, %d allo)',
                len(merged), total_ortho, total_allo)

    # ── [2/4] RECAP fragmentation (or load from cache) ───────────────────────
    if io_utils.output_is_cached(cache_path, args.force):
        print('\n[2/4] Loading fragmentation cache (--force to rebuild) ...')
        print(f'  {cache_path}')
        logger.info('[2/4] Loading fragmentation cache from %s', cache_path)
        inverted, frag_hac = load_pickle(cache_path)
        n_unique = len(inverted)
        print(f'  {n_unique:,} unique fragments loaded from cache')
        logger.info('Cache loaded: %d unique fragments', n_unique)
    else:
        items = list(zip(merged['ccd_code'], merged['canonical_smiles']))

        n_workers_desc = args.workers or os.cpu_count()
        mode = ('parallel, %d workers' % n_workers_desc
                if len(items) >= 200 and n_workers_desc > 1 else 'serial')
        print(f'\n[2/4] RECAP fragmentation ({mode}, timeout={_RECAP_TIMEOUT_S}s/mol) ...')
        logger.info('[2/4] RECAP fragmentation (%s)', mode)

        inverted, frag_hac, rstats = run_recap_batch(
            items, args.min_hac, args.max_hac, n_workers=args.workers)

        n_unique = len(inverted)
        print(f'  {n_unique:,} unique fragments collected '
              f'({rstats["n_capped"]:,} cap hits, {rstats["n_timeout"]:,} timeouts, '
              f'{rstats["n_failed"]:,} failures)')
        logger.info('%d unique fragments; %d cap hits; %d timeouts; %d failures',
                    n_unique, rstats['n_capped'], rstats['n_timeout'], rstats['n_failed'])

        save_pickle((inverted, frag_hac), cache_path)
        print(f'  fragmentation cache saved → {cache_path}')
        logger.info('Fragmentation cache saved to %s', cache_path)

    # ── [3/4] Build per-fragment counts and UniProt stats ─────────────────────
    print(f'\n[3/4] Computing counts and UniProt stats for {n_unique:,} fragments ...')
    logger.info('[3/4] Computing counts/UniProt for %d fragments', n_unique)
    prop_rows = []
    for frag_smi, ccd_set in inverted.items():
        ortho_ccds = {c for c in ccd_set if cat_map.get(c) == 'orthosteric'}
        allo_ccds  = {c for c in ccd_set if cat_map.get(c) == 'allosteric'}
        n_ortho = len(ortho_ccds)
        n_allo  = len(allo_ccds)

        ortho_uniprots = {u for c in ortho_ccds for u in ccd_to_uniprots.get(c, set())}
        allo_uniprots  = {u for c in allo_ccds  for u in ccd_to_uniprots.get(c, set())}

        prop_rows.append({
            'fragment_smiles':  frag_smi,
            'hac':              frag_hac.get(frag_smi, 0),
            'n_ortho':          n_ortho,
            'n_allo':           n_allo,
            'n_total':          n_ortho + n_allo,
            'n_uniprot_ortho':  len(ortho_uniprots),
            'n_uniprot_allo':   len(allo_uniprots),
        })

    frag_df = pd.DataFrame(prop_rows)

    # Assign fragment_index (1-based, sorted by total frequency descending)
    frag_df = frag_df.sort_values('n_total', ascending=False, kind='stable').reset_index(drop=True)
    frag_df['fragment_index'] = frag_df.index + 1

    # Record the thresholds used for this run so downstream consumers
    # can report them without having to re-derive them from the data (the observed
    # min/max of `hac` isn't the same as args.min_hac/max_hac — the filter keeps
    # min_hac < hac <= max_hac, so the true lower bound isn't recoverable from data alone).
    frag_df['min_hac'] = args.min_hac
    frag_df['max_hac'] = args.max_hac

    # ── [4/4] Save outputs ────────────────────────────────────────────────────
    print('\n[4/4] Saving outputs ...')
    logger.info('[4/4] Saving outputs')

    # CSVs
    out_cols = ['fragment_index', 'fragment_smiles', 'hac',
                'n_ortho', 'n_allo', 'n_total',
                'n_uniprot_ortho', 'n_uniprot_allo',
                'min_hac', 'max_hac']
    frag_df[out_cols].to_csv(out_frags, index=False)

    idx_map = frag_df.set_index('fragment_smiles')['fragment_index'].to_dict()
    assoc_rows = []
    for frag_smi, ccd_set in inverted.items():
        fidx = idx_map.get(frag_smi)
        for ccd in ccd_set:
            assoc_rows.append({
                'exclusive_category': cat_map.get(ccd, ''),
                'fragment_index':     fidx,
                'fragment_smiles':    frag_smi,
                'ccd_code':           ccd,
            })
    assoc_df = pd.DataFrame(assoc_rows).sort_values(['fragment_index', 'ccd_code'])
    assoc_df.to_csv(out_assoc, index=False)

    print(f'  recap_fragments.csv             → {out_frags}')
    print(f'  recap_fragment_associations.csv → {out_assoc}  ({len(assoc_df):,} rows)')
    logger.info('Saved recap_fragments.csv (%d fragments)', len(frag_df))
    logger.info('Saved recap_fragment_associations.csv (%d rows)', len(assoc_df))

    cli_utils.finish_run(t0)


if __name__ == '__main__':
    main()
