"""
cluster_protein_similarity.py
==============================
Cluster binding-site rows by protein structural similarity, using Foldseek over
the assembly CIF files in --cif-dir. Two binding-site rows cluster together only if 
they use the exact same set of Foldseek structural clusters across all their chains.

ALL rows and allosteric-only rows are clustered independently, so cluster IDs
are file-local: cluster_struc=3 means different things in clusters_all.csv and
clusters_allosteric.csv.

Foldseek defaults here: --cluster-mode 1 (connected-component) with -c 0.8 /
--cov-mode 0 / -e 0.01 / --tmscore-threshold 0.5. See --help to override.

SCOP2 and CATH annotate rows for the largest_clusters.csv breakdown and never
gate clustering.

Inputs
------
  data/site_ligand_features.csv               (the BUILD step's output)
  {cif_dir}/{PDB_ID}-assembly1.cif            (missing files -> singletons)
  collection/result/mapping/mapped_{scop,cath}.csv

Outputs
-------
  data/clusters_all.csv, data/clusters_allosteric.csv
  data/foldseek_{ALL,allosteric}_cluster.tsv  (easy-cluster's raw adjacency list)
  data/scop2_sf_names.csv, data/cath-names.txt  (annotation caches)
  data/largest_clusters.csv                   (per-cluster SCOP2/CATH breakdown)

Usage
-----
    python cluster_protein_similarity.py [--verbose] [--force]
    python cluster_protein_similarity.py --cif-dir /path/to/cifs --threads 8
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import config  # must precede cli_utils/io_utils: injects collection/src onto sys.path
import cli_utils
import io_utils
from utils import UnionFind, create_session
from site_pipeline_utils import parse_site_res_list

logger = logging.getLogger(__name__)

SCOP2_API_ANCESTRY_URL = 'https://www.ebi.ac.uk/pdbe/scop/api/ancestry/{}'
SCOP2_FETCH_WORKERS = 8

CATH_NAMES_URL = ('https://download.cathdb.info/cath/releases/latest-release/'
                   'cath-classification-data/cath-names.txt')

LARGEST_CLUSTERS_FILE = os.path.join(config.DATA_DIR, 'largest_clusters.csv')


# ─────────────────────────────────────────────────────────────────────────────
# Chain identification (no CIF access needed for this)
# ─────────────────────────────────────────────────────────────────────────────

def _all_chains(binding_site: str) -> list[str]:
    """Every distinct chain contributing to a site, sorted.

    All of them count for clustering — a heterodimer's partner chain matters as
    much as its most-populated one.
    """
    chains = {chain for chain, _resnum in parse_site_res_list(binding_site)}
    return sorted(chains)


def add_chain_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Copy of df with an added all_chains (list[str]) column."""
    df = df.copy()
    df['all_chains'] = df['binding_site'].apply(_all_chains)
    n_multi = (df['all_chains'].apply(len) > 1).sum()
    print(f'  {n_multi:,} / {len(df):,} rows span more than one chain '
          f'({100 * n_multi / max(len(df), 1):.2f}%).')
    logger.info('%d / %d rows are multi-chain sites', n_multi, len(df))
    return df


# ─────────────────────────────────────────────────────────────────────────────
# SCOP2 SF ID resolution + name fetching (side-channel, never gates clustering)
# ─────────────────────────────────────────────────────────────────────────────

def load_scop2_lookup(scop_csv: str) -> dict[tuple[str, str], set[str]]:
    """(PDB, chain) -> {SCOP2B domain ID} from mapped_scop.csv.

    Filtered to xref_db == 'SCOP2B' and domain_name == 'SF'. Despite the 'SF'
    label these accessions are per-structure *domain* IDs, not superfamily IDs —
    reaching the superfamily needs fetch_scop2_superfamily_info. Multi-domain
    chains keep all their IDs. Missing file -> empty dict, never aborts.
    """
    if not os.path.exists(scop_csv):
        print(f'  WARNING: {scop_csv} not found — SCOP2 side-channel will be empty.')
        logger.warning('mapped_scop.csv not found at %s — scop2_sf_ids will be empty', scop_csv)
        return {}

    df = pd.read_csv(scop_csv, dtype=str)
    df = df[(df['xref_db'] == 'SCOP2B') & (df['domain_name'] == 'SF')]

    lookup: dict[tuple[str, str], set[str]] = {}
    for pdb, chain, sf_id in zip(df['pdb_id'].str.upper(), df['asym_id'], df['xref_db_acc']):
        lookup.setdefault((pdb, chain), set()).add(sf_id)
    return lookup


def resolve_scop2_domain_ids(df: pd.DataFrame, lookup: dict[tuple[str, str], set[str]]) -> pd.Series:
    """Per-row frozenset of SCOP2B domain IDs, unioned across all_chains.

    Internal only — these get mapped to superfamily IDs before anything is
    written. Rows with no resolvable chain get an empty set, never dropped.
    """
    pdbs_upper = df['source_pdb'].str.upper()
    domain_sets = [
        frozenset().union(*(lookup.get((pdb, c), frozenset()) for c in chains)) if chains else frozenset()
        for pdb, chains in zip(pdbs_upper, df['all_chains'])
    ]
    n_resolved = sum(1 for s in domain_sets if s)
    print(f'  {n_resolved:,} / {len(df):,} rows resolve to >=1 SCOP2 domain ID.')
    logger.info('%d / %d rows resolve to >=1 SCOP2 domain ID', n_resolved, len(df))
    return pd.Series(domain_sets, index=df.index)


# ─────────────────────────────────────────────────────────────────────────────
# CATH side-channel (best-effort fallback for SCOP2-unresolved rows)
# ─────────────────────────────────────────────────────────────────────────────

def load_cath_lookup(cath_csv: str) -> dict[tuple[str, str], set[str]]:
    """(PDB, chain) -> {full 4-level CATH code} from mapped_cath.csv.

    Every row is already a domain-granularity code (e.g. '1.10.530.40'); they
    are truncated to the topology level later. Missing file -> empty dict.
    """
    if not os.path.exists(cath_csv):
        print(f'  WARNING: {cath_csv} not found — CATH side-channel will be empty.')
        logger.warning('mapped_cath.csv not found at %s — cath_topology_ids will be empty', cath_csv)
        return {}

    df = pd.read_csv(cath_csv, dtype=str)

    lookup: dict[tuple[str, str], set[str]] = {}
    for pdb, chain, code in zip(df['pdb_id'].str.upper(), df['asym_id'], df['xref_db_acc']):
        lookup.setdefault((pdb, chain), set()).add(code)
    return lookup


def resolve_cath_topology_ids(df: pd.DataFrame, lookup: dict[tuple[str, str], set[str]]) -> pd.Series:
    """Per-row frozenset of CATH topology codes, unioned across all_chains.

    Codes are truncated to 3 levels (Class.Architecture.Topology) — the deepest
    level cath-names.txt names. Rows with no resolvable chain get an empty set.
    """
    pdbs_upper = df['source_pdb'].str.upper()
    topology_sets = []
    for pdb, chains in zip(pdbs_upper, df['all_chains']):
        codes = (frozenset().union(*(lookup.get((pdb, c), frozenset()) for c in chains))
                 if chains else frozenset())
        topology_sets.append(frozenset('.'.join(code.split('.')[:3]) for code in codes))
    n_resolved = sum(1 for s in topology_sets if s)
    print(f'  {n_resolved:,} / {len(df):,} rows resolve to >=1 CATH topology ID.')
    logger.info('%d / %d rows resolve to >=1 CATH topology ID', n_resolved, len(df))
    return pd.Series(topology_sets, index=df.index)


def fetch_cath_names(cache_path: str) -> dict[str, str]:
    """{topology code: description} from CATH's cath-names.txt.

    Whole-file cache — downloaded once, then read from cache_path forever
    (delete it to refresh). A download failure degrades to an empty dict.
    """
    if os.path.exists(cache_path):
        print(f'  Using cached CATH names file → {cache_path}')
        with open(cache_path) as f:
            text = f.read()
    else:
        print(f'  Downloading CATH names file from {CATH_NAMES_URL} …')
        try:
            session = create_session()
            r = session.get(CATH_NAMES_URL, timeout=60)
            r.raise_for_status()
            text = r.text
        except Exception as exc:
            print(f'  WARNING: failed to download CATH names file — CATH side-channel will be empty ({exc}).')
            logger.warning('CATH names file download failed: %s', exc)
            return {}
        with open(cache_path, 'w') as f:
            f.write(text)
        print(f'  Saved → {cache_path}')

    names: dict[str, str] = {}
    for line in text.splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.split(None, 2)
        if len(parts) < 3:
            continue
        node_id, _rep_domain, description = parts
        if node_id.count('.') != 2:   # only topology-level (3-part) nodes
            continue
        names[node_id] = description.strip().lstrip(':').strip()
    print(f'  Loaded {len(names):,} CATH topology-level names.')
    return names


def _fetch_one_scop2_superfamily(domain_id: str, session) -> tuple[str, str] | None:
    """Walk one SCOP2B domain ID up to its superfamily (id, name) via the EBI API.

    Domain IDs are per-structure records, so the superfamily node is what gives
    a shared, meaningful grouping. Returns None on a request failure, an
    {"error": ...} body (the API's HTTP-200 way of rejecting an ID), or a
    lineage with no superfamily — unresolved IDs are retried next run, never
    cached as permanent failures.
    """
    try:
        r = session.get(SCOP2_API_ANCESTRY_URL.format(domain_id), timeout=20)
        r.raise_for_status()
        data = r.json()
    except Exception as exc:
        logger.warning('SCOP2 ancestry fetch failed for domain %s: %s', domain_id, exc)
        return None
    if 'error' in data:
        logger.warning('SCOP2 ancestry error for domain %s: %s', domain_id, data['error'])
        return None
    nodes = data.get('lineage', {}).get('nodes', {})
    for node in nodes.values():
        if node.get('type') == 'superfamily':
            name = node.get('name')
            if name:
                return str(node['id']), name
    logger.warning('No superfamily node in SCOP2 ancestry for domain %s', domain_id)
    return None


def fetch_scop2_superfamily_info(domain_ids: set[str], cache_path: str) -> dict[str, tuple[str, str]]:
    """domain ID -> (superfamily ID, description), fetched incrementally and cached.

    Only IDs absent from scop2_sf_names.csv cost a round-trip. Threaded, since
    the API takes one ID per request. IDs that fail to resolve are simply
    absent from the result — never a placeholder.
    """
    cached: dict[str, tuple[str, str]] = {}
    if os.path.exists(cache_path):
        cache_df = pd.read_csv(cache_path, dtype=str, keep_default_na=False)
        cached = {
            domain_id: (sf_id, description)
            for domain_id, sf_id, description in zip(
                cache_df['domain_id'], cache_df['sf_id'], cache_df['description'])
        }

    missing = sorted(domain_ids - cached.keys())
    if not missing:
        print(f'  All {len(domain_ids):,} SCOP2 domain IDs already cached — no fetch needed.')
        return {domain_id: cached[domain_id] for domain_id in domain_ids if domain_id in cached}

    print(f'  Fetching {len(missing):,} new SCOP2 superfamily mappings from the EBI API '
          f'({len(domain_ids) - len(missing):,} already cached)…')
    logger.info('Fetching %d new SCOP2 superfamily mappings (%d already cached)',
                len(missing), len(domain_ids) - len(missing))

    session = create_session()
    fetched: dict[str, tuple[str, str]] = {}
    with ThreadPoolExecutor(max_workers=SCOP2_FETCH_WORKERS) as executor:
        futures = {executor.submit(_fetch_one_scop2_superfamily, domain_id, session): domain_id
                   for domain_id in missing}
        n_done = 0
        for future in as_completed(futures):
            domain_id = futures[future]
            result = future.result()
            if result:
                fetched[domain_id] = result
            n_done += 1
            if n_done % 500 == 0 or n_done == len(missing):
                print(f'    [{n_done}/{len(missing)}] SCOP2 superfamily mappings fetched')

    n_failed = len(missing) - len(fetched)
    if n_failed:
        print(f'  WARNING: {n_failed:,} / {len(missing):,} SCOP2 domain IDs failed to resolve a superfamily.')
        logger.warning('%d / %d new SCOP2 domain IDs failed to resolve', n_failed, len(missing))

    merged = {**cached, **fetched}
    out_df = pd.DataFrame(
        [(domain_id, sf_id, description) for domain_id, (sf_id, description) in sorted(merged.items())],
        columns=['domain_id', 'sf_id', 'description'])
    out_df.to_csv(cache_path, index=False)
    logger.info('Wrote %d SCOP2 superfamily mappings -> %s', len(out_df), cache_path)
    print(f'  Saved → {cache_path}')

    return {domain_id: merged[domain_id] for domain_id in domain_ids if domain_id in merged}


# ─────────────────────────────────────────────────────────────────────────────
# Shared clustering helpers
# ─────────────────────────────────────────────────────────────────────────────

def _union_by_exact_key(uf: UnionFind, keys: list[frozenset | None]) -> None:
    """Union indices sharing an identical key; empty keys stay singletons.

    Exact set matching, not "any shared element" — deliberate.
    """
    buckets: dict[frozenset, list[int]] = {}
    for i, key in enumerate(keys):
        if not key:
            continue
        buckets.setdefault(key, []).append(i)
    for indices in buckets.values():
        if len(indices) < 2:
            continue
        for j in indices[1:]:
            uf.union(indices[0], j)


def _save_intermediates(out_prefix: str, dest_tag: str) -> None:
    """Copy an easy-cluster run's _cluster.tsv into data/ as {dest_tag}_cluster.tsv."""
    src = f'{out_prefix}_cluster.tsv'
    if not os.path.exists(src):
        return
    dest = os.path.join(config.DATA_DIR, f'{dest_tag}_cluster.tsv')
    shutil.copyfile(src, dest)
    logger.info('Copied intermediate %s -> %s', src, dest)


# ─────────────────────────────────────────────────────────────────────────────
# Structure-based clustering (Foldseek)
# ─────────────────────────────────────────────────────────────────────────────

_CIF_SUFFIX = '-assembly1.cif'
_CIF_STEM_SUFFIX = '-assembly1'  # member IDs from foldseek have the extension stripped


def cluster_by_structure(scope_df: pd.DataFrame, foldseek_bin: str, cif_dir: str,
                          coverage: float | None, cov_mode: int | None,
                          tmscore_threshold: float | None, tmscore_threshold_mode: int | None,
                          evalue: float | None, cluster_mode: int, threads: int | None,
                          scope_name: str) -> UnionFind:
    """Union rows whose chains fall in the same Foldseek structural clusters.

    Foldseek splits each multi-chain assembly file into one database entry per
    chain; --chain-name-mode 1 suffixes each entry with _{chain} so member IDs
    map back to (source_pdb, chain) unambiguously.

    Rows are unioned only on an exact match of their chains' cluster-rep ID sets,
    so a heterodimer merges only when both chains match.

    cluster_mode=1 additionally requires --single-step-clustering 1, added
    automatically — easy-cluster's cascaded workflow hard-errors without it.
    The filter arguments are passed through only when not None, so Foldseek's
    own defaults apply otherwise.

    Returns a UnionFind over range(len(scope_df)).
    """
    uf = UnionFind(range(len(scope_df)))

    needed_pdbs = sorted(scope_df['source_pdb'].str.upper().unique())
    available_pdbs = [
        pdb for pdb in needed_pdbs
        if os.path.exists(os.path.join(cif_dir, f'{pdb}{_CIF_SUFFIX}'))
    ]
    n_missing = len(needed_pdbs) - len(available_pdbs)
    if n_missing:
        print(f'  {n_missing:,} / {len(needed_pdbs):,} unique source_pdb have no '
              f'CIF file in {cif_dir} — those rows get a structural singleton.')
        logger.info('%d / %d unique source_pdb missing from cif_dir %s',
                    n_missing, len(needed_pdbs), cif_dir)

    if not available_pdbs:
        print(f'  WARNING: zero CIF files found in cif_dir — all rows singleton in [{scope_name}].')
        logger.warning('Zero CIF files available in cif_dir %s — skipping foldseek subprocess', cif_dir)
        return uf

    with tempfile.TemporaryDirectory(prefix='foldseek_input_') as input_dir, \
         tempfile.TemporaryDirectory(prefix='foldseek_') as tmp:
        for pdb in available_pdbs:
            src = os.path.join(cif_dir, f'{pdb}{_CIF_SUFFIX}')
            os.symlink(os.path.abspath(src), os.path.join(input_dir, f'{pdb}{_CIF_SUFFIX}'))

        out_prefix = os.path.join(tmp, 'clust')
        cmd = [
            foldseek_bin, 'easy-cluster', input_dir, out_prefix, tmp,
            '--cluster-mode', str(cluster_mode),
            '--chain-name-mode', '1',
            '-v', '1',
        ]
        if cluster_mode == 1:
            cmd.extend(['--single-step-clustering', '1'])
        if coverage is not None:
            cmd.extend(['-c', str(coverage)])
        if cov_mode is not None:
            cmd.extend(['--cov-mode', str(cov_mode)])
        if tmscore_threshold is not None:
            cmd.extend(['--tmscore-threshold', str(tmscore_threshold)])
        if tmscore_threshold_mode is not None:
            cmd.extend(['--tmscore-threshold-mode', str(tmscore_threshold_mode)])
        if evalue is not None:
            cmd.extend(['-e', str(evalue)])
        if threads is not None:
            cmd.extend(['--threads', str(threads)])
        logger.info('Running: %s', ' '.join(cmd))
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            logger.error('foldseek easy-cluster failed (rc=%d): %s', exc.returncode, exc.stderr)
            print(f'  ERROR: foldseek easy-cluster failed (rc={exc.returncode}):\n{exc.stderr}')
            raise

        _save_intermediates(out_prefix, f'foldseek_{scope_name}')

        cluster_tsv = f'{out_prefix}_cluster.tsv'
        member_to_rep: dict[str, str] = {}
        with open(cluster_tsv) as f:
            for line in f:
                line = line.rstrip('\n')
                if not line:
                    continue
                rep_id, member_id = line.split('\t')
                member_to_rep[member_id] = rep_id

    # member_id format: '{PDB_ID}-assembly1_{chain}' (chain always present
    # thanks to --chain-name-mode 1).
    pdb_chain_to_rep: dict[tuple[str, str], str] = {}
    for member_id, rep_id in member_to_rep.items():
        stem, _, chain = member_id.rpartition('_')
        pdb = stem[:-len(_CIF_STEM_SUFFIX)] if stem.endswith(_CIF_STEM_SUFFIX) else stem
        pdb_chain_to_rep[(pdb, chain)] = rep_id

    pdbs_upper = scope_df['source_pdb'].str.upper()
    row_keys = [
        frozenset(pdb_chain_to_rep[(pdb, c)] for c in chains if (pdb, c) in pdb_chain_to_rep)
        for pdb, chains in zip(pdbs_upper, scope_df['all_chains'])
    ]
    _union_by_exact_key(uf, row_keys)

    return uf


# ─────────────────────────────────────────────────────────────────────────────
# Deterministic ID assignment
# ─────────────────────────────────────────────────────────────────────────────

def assign_deterministic_ids(uf: UnionFind, scope_df: pd.DataFrame) -> list[int]:
    """Sequential cluster IDs, one per row, assigned deterministically.

    Foldseek rep IDs and Python set ordering aren't stable across runs, so
    groups are sorted by their minimum (source_pdb, ligand_code) first.
    """
    groups: dict[int, list[int]] = {}
    for i in range(len(scope_df)):
        groups.setdefault(uf.find(i), []).append(i)

    def group_key(indices: list[int]) -> tuple:
        return min(
            (scope_df.iloc[i]['source_pdb'], scope_df.iloc[i]['ligand_code'])
            for i in indices
        )

    ordered_roots = sorted(groups, key=lambda root: group_key(groups[root]))
    root_to_id = {root: cid for cid, root in enumerate(ordered_roots)}

    ids = [0] * len(scope_df)
    for root, indices in groups.items():
        for i in indices:
            ids[i] = root_to_id[root]
    return ids


# ─────────────────────────────────────────────────────────────────────────────
# Per-scope clustering
# ─────────────────────────────────────────────────────────────────────────────

_CLUSTER_MODE_LABELS = {0: 'greedy-set-cover', 1: 'connected-component', 2: 'greedy-incremental'}


def cluster_scope(scope_df: pd.DataFrame, scope_name: str, foldseek_bin: str, cif_dir: str,
                   struct_coverage: float | None, struct_cov_mode: int | None,
                   tmscore_threshold: float | None, tmscore_threshold_mode: int | None,
                   evalue: float | None, struct_cluster_mode: int,
                   threads: int | None) -> pd.DataFrame:
    """Cluster one scope (ALL or allosteric-only), one output row per input row."""
    scope_df = scope_df.reset_index(drop=True)
    n = len(scope_df)
    print(f'  [{scope_name}] {n:,} rows.')

    uf = cluster_by_structure(scope_df, foldseek_bin, cif_dir, struct_coverage,
                               struct_cov_mode, tmscore_threshold=tmscore_threshold,
                               tmscore_threshold_mode=tmscore_threshold_mode,
                               evalue=evalue, cluster_mode=struct_cluster_mode,
                               threads=threads, scope_name=scope_name)

    cluster_struc = assign_deterministic_ids(uf, scope_df)

    n_clusters = len(set(cluster_struc))
    mode_label = _CLUSTER_MODE_LABELS.get(struct_cluster_mode, f'mode-{struct_cluster_mode}')
    print(f'  [{scope_name}] clusters — {mode_label}: {n_clusters:,}')

    out = pd.DataFrame({
        'ID':                scope_df['ID'],
        'protein_name_pdb':  scope_df['protein_name_pdb'],
        'source_pdb':        scope_df['source_pdb'],
        'ligand_code':       scope_df['ligand_code'],
        'site_category':     scope_df['site_category'],
        'uniprot_ac':        scope_df['uniprot_ac'],
        'cluster_struc':     cluster_struc,
        'scop2_sf_ids':      scope_df['scop2_sf_ids'].apply(lambda s: ';'.join(sorted(s))),
        'cath_topology_ids': scope_df['cath_topology_ids'].apply(lambda s: ';'.join(sorted(s))),
    })
    assert len(out) == n, 'row count must be preserved exactly'
    return out


# ─────────────────────────────────────────────────────────────────────────────
# largest_clusters.csv: per-cluster SCOP2 ID breakdown (report-facing output)
# ─────────────────────────────────────────────────────────────────────────────

def _scop2_top_groups_rows(out_df: pd.DataFrame, scope_name: str, sf_names: dict[str, str],
                            cath_names: dict[str, str]) -> list[dict]:
    """Rows describing the largest clusters' most frequent SCOP2 ID groups.

    "Largest" is ranked twice (rank_metric 'pdb'/'uniprot'), by unique source_pdb
    and by unique uniprot_ac, since ranking on raw row count over-credits clusters
    dominated by one multi-site structure. Per cluster, the top 5 SCOP2 groups are
    emitted in descending count order; a group is one whole scop2_sf_ids string
    (multi-chain values like '8033861;8036511' stay intact), 'NA' if unresolved.
    `percentage` is against cluster_size_sites. Each group also carries its most
    common CATH topology as a cross-check for where SCOP2 resolved nothing.

    Returns a list of dicts: scope, rank_metric, cluster_id, cluster_size_sites,
    cluster_size_pdb, cluster_size_uniprot, count, percentage, sf_ids, sf_names,
    cath_topology_ids, cath_topology_names.
    """
    cluster_size_sites   = out_df.groupby('cluster_struc').size()
    cluster_size_pdb     = out_df.groupby('cluster_struc')['source_pdb'].nunique()
    cluster_size_uniprot = out_df.groupby('cluster_struc')['uniprot_ac'].nunique()

    rows: list[dict] = []
    for rank_metric, size_series in [('pdb', cluster_size_pdb), ('uniprot', cluster_size_uniprot)]:
        top_clusters = size_series.nlargest(10)

        for cluster_id in top_clusters.index:
            cluster_df = out_df.loc[out_df['cluster_struc'] == cluster_id]
            sf_values = cluster_df['scop2_sf_ids'].fillna('')
            top5 = sf_values.value_counts().head(5)

            for sf_value, count in top5.items():
                if sf_value:
                    sf_ids_field = sf_value
                    sf_names_field = '; '.join(sf_names.get(sf_id, '') for sf_id in sf_value.split(';'))
                else:
                    sf_ids_field = 'NA'
                    sf_names_field = 'NA'

                group_cath_values = cluster_df.loc[sf_values == sf_value, 'cath_topology_ids']
                cath_counts = group_cath_values[group_cath_values != ''].value_counts()
                if len(cath_counts):
                    # min() breaks count ties deterministically (value_counts' own
                    # tie order isn't guaranteed stable across runs/pandas versions).
                    top_count = cath_counts.iloc[0]
                    cath_ids_field = min(cath_counts[cath_counts == top_count].index)
                    cath_names_field = '; '.join(cath_names.get(c, '') for c in cath_ids_field.split(';'))
                else:
                    cath_ids_field = 'NA'
                    cath_names_field = 'NA'

                rows.append({
                    'scope':                scope_name,
                    'rank_metric':          rank_metric,
                    'cluster_id':           cluster_id,
                    'cluster_size_sites':   int(cluster_size_sites[cluster_id]),
                    'cluster_size_pdb':     int(cluster_size_pdb[cluster_id]),
                    'cluster_size_uniprot': int(cluster_size_uniprot[cluster_id]),
                    'count':                count,
                    'percentage':           round(count / cluster_size_sites[cluster_id] * 100, 1),
                    'sf_ids':               sf_ids_field,
                    'sf_names':             sf_names_field,
                    'cath_topology_ids':    cath_ids_field,
                    'cath_topology_names':  cath_names_field,
                })

    return rows


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    """Cluster site_ligand_features.csv by Foldseek structural similarity, two scopes.

    Resolves the SCOP2 and CATH side-channel annotations first (both
    best-effort), then clusters ALL rows and allosteric-only rows independently
    and writes the cluster tables plus largest_clusters.csv.
    """
    parser = argparse.ArgumentParser(
        description='Cluster binding-site rows by protein structural similarity '
                    '(Foldseek only, default --cluster-mode 1 / connected-component).')
    parser.add_argument('--foldseek-bin', default='foldseek',
                        help='Path to foldseek executable (default: foldseek on PATH)')
    parser.add_argument('--cif-dir', default=config.CIF_DIR,
                        help='Directory of {PDB_ID}-assembly1.cif files (default: config.CIF_DIR). '
                             'Rows with no file here become singletons, never dropped.')
    parser.add_argument('--struct-coverage', type=float, default=0.8,
                        help='foldseek -c coverage (default: 0.8). Bounds alignment span, not '
                             'quality — see --tmscore-threshold.')
    parser.add_argument('--struct-cov-mode', type=int, default=0,
                        help='foldseek --cov-mode (default: 0, i.e. coverage of query AND target)')
    parser.add_argument('--struct-cluster-mode', type=int, default=1,
                        help='foldseek --cluster-mode (default: 1, connected-component; 0 for '
                             'greedy set-cover). --single-step-clustering 1 is added when this is 1.')
    parser.add_argument('--evalue', type=float, default=0.01,
                        help='foldseek -e E-value cutoff (default: 0.01, foldseek\'s own is 10.0)')
    parser.add_argument('--tmscore-threshold', type=float, default=0.5,
                        help='foldseek --tmscore-threshold (default: 0.5, foldseek\'s own is 0.0). '
                             'The alignment-quality filter; 0.0 disables it.')
    parser.add_argument('--tmscore-threshold-mode', type=int, default=0,
                        help='foldseek --tmscore-threshold-mode: normalise TM by 0 alignment, '
                             '1 query, or 2 target length (default: 0)')
    parser.add_argument('--threads', type=int, default=None,
                        help='Threads for foldseek (default: unset, foldseek uses all cores). '
                             '--threads 1 is more reproducible across runs.')
    cli_utils.add_common_args(parser, force=True)
    args = parser.parse_args()

    if not args.cif_dir:
        print('  WARNING: no --cif-dir provided (config.CIF_DIR is unset on this machine) — '
              'the Foldseek structural axis has no CIF files to read, so cluster_struc will be '
              'all singletons. Set CIF_DIR in src/config_local.py, the STRUCTURE_CIF_DIR env '
              'var, or pass --cif-dir to get real structural clusters.')
        logger.warning('No --cif-dir provided — cluster columns will be all singletons for every row')

    t0 = cli_utils.init_run('cluster_protein_similarity', args, config.LOG_DIR)

    site_features_csv = os.path.join(config.DATA_DIR, 'site_ligand_features.csv')
    out_all_path  = os.path.join(config.DATA_DIR, 'clusters_all.csv')
    out_allo_path = os.path.join(config.DATA_DIR, 'clusters_allosteric.csv')

    if (io_utils.output_is_cached(out_all_path, args.force)
            and io_utils.output_is_cached(out_allo_path, args.force)
            and io_utils.output_is_cached(LARGEST_CLUSTERS_FILE, args.force)):
        print('[1/1] clusters_all.csv, clusters_allosteric.csv, and largest_clusters.csv '
              'already exist — SKIP (--force to recompute).')
        cli_utils.finish_run(t0)
        return

    # ── [1/4] Load inputs ────────────────────────────────────────────────────
    print('\n[1/4] Loading site_ligand_features.csv …')
    df = pd.read_csv(site_features_csv)
    print(f'  {len(df):,} rows.')

    # ── [2/4] Chain identification + SCOP2 SF ID resolution ──────────────────
    print('\n[2/4] Identifying every chain contributing to each site …')
    df = add_chain_columns(df)

    print('\n  Resolving SCOP2 domain IDs (side-channel) …')
    scop2_lookup = load_scop2_lookup(config.SCOP_MAPPING_CSV)
    df['scop2_domain_ids'] = resolve_scop2_domain_ids(df, scop2_lookup)

    print('\n  Resolving SCOP2 superfamily IDs from domain IDs …')
    all_domain_ids = {domain_id for domain_set in df['scop2_domain_ids'] for domain_id in domain_set}
    scop2_names_path = os.path.join(config.DATA_DIR, 'scop2_sf_names.csv')
    domain_to_sf = fetch_scop2_superfamily_info(all_domain_ids, scop2_names_path)
    df['scop2_sf_ids'] = df['scop2_domain_ids'].apply(
        lambda ids: frozenset(domain_to_sf[i][0] for i in ids if i in domain_to_sf))

    print('\n  Resolving CATH topology-level IDs (side-channel) …')
    cath_lookup = load_cath_lookup(config.CATH_MAPPING_CSV)
    df['cath_topology_ids'] = resolve_cath_topology_ids(df, cath_lookup)

    print('\n  Fetching CATH topology-level names …')
    cath_names_path = os.path.join(config.DATA_DIR, 'cath-names.txt')
    cath_names = fetch_cath_names(cath_names_path)

    # ── [3/4] Cluster both scopes (structure only) ────────────────────────────
    _mode_label = _CLUSTER_MODE_LABELS.get(args.struct_cluster_mode, f'mode-{args.struct_cluster_mode}')
    print(f'\n[3/4] Clustering (Foldseek structure, {_mode_label}) …')
    df_all  = df
    df_allo = df[df['site_category'] == 'allosteric']

    clusters_all = cluster_scope(
        df_all, 'ALL', args.foldseek_bin, args.cif_dir,
        args.struct_coverage, args.struct_cov_mode, args.tmscore_threshold,
        args.tmscore_threshold_mode, args.evalue,
        args.struct_cluster_mode, args.threads)

    clusters_allo = cluster_scope(
        df_allo, 'allosteric', args.foldseek_bin, args.cif_dir,
        args.struct_coverage, args.struct_cov_mode, args.tmscore_threshold,
        args.tmscore_threshold_mode, args.evalue,
        args.struct_cluster_mode, args.threads)

    # ── [4/4] Write outputs ───────────────────────────────────────────────────
    print('\n[4/4] Writing outputs …')
    io_utils.save_csv(clusters_all, out_all_path, logger)
    io_utils.save_csv(clusters_allo, out_allo_path, logger)

    sf_names = {sf_id: description for sf_id, description in domain_to_sf.values() if sf_id and description}

    largest_clusters_rows = []
    for scope_name, out_df in [('ALL', clusters_all), ('allosteric', clusters_allo)]:
        largest_clusters_rows.extend(_scop2_top_groups_rows(out_df, scope_name, sf_names, cath_names))
    largest_clusters_df = pd.DataFrame(largest_clusters_rows)
    io_utils.save_csv(largest_clusters_df, LARGEST_CLUSTERS_FILE, logger)

    cli_utils.finish_run(t0)


if __name__ == '__main__':
    main()
