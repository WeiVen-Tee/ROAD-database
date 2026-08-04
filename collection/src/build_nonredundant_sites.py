"""
Build a non-redundant set of binding sites.

Modes:
  annotated (default) — only sites with a computed functional-site distance
  complete            — all detected binding sites

Each site is represented as an ordered AA one-letter tuple. Sites within the same
(UniProt AC, ligand) group are clustered by normalised Levenshtein distance; one
representative per cluster is kept (most binding residues; ties by PDB ID).

Usage:
    python build_nonredundant_sites.py [--mode annotated|complete] [--cif-dir DIR] [-v]
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

import gemmi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    CIF_DIR, EDIT_DIST_THRESHOLD, MIN_SITE_RESIDUES, LOG_DIR,
    RESULT_MAPPING_DIR, RESULT_LIGANDS_DIR, RESULT_BINDING_DIR,
    RESULT_DISTANCES_DIR, RESULT_NONREDUNDANT_DIR,
)
from utils import UnionFind, cif2seq, load_pickle, save_pickle, setup_logging

log = None  # set in main()


def _parse_ligand_key(key: str) -> tuple:
    """'1abc:1.A:101:ATP' -> ('1abc', '1.A', '101', 'ATP')."""
    parts = key.split(':')
    return parts[0], parts[1], parts[2], parts[3]


def _ligand_code(key: str) -> str:
    """'1abc:1.A:101:ATP' -> 'A.101:ATP' (model prefix dropped)."""
    _, model_chain, resnum, resname = _parse_ligand_key(key)
    chain = model_chain.split('.')[1]
    return f'{chain}.{resnum}:{resname}'


def _parse_res_id(res_id: str) -> tuple:
    """'1.A.100A' -> ('1.A', '100A'); splits at the second dot."""
    dot2 = res_id.index('.', res_id.index('.') + 1)
    return res_id[:dot2], res_id[dot2 + 1:]


def _split_resnum(resnum_str: str) -> tuple | None:
    """'42' -> (42, ''); '100A' -> (100, 'A'); '-5B' -> (-5, 'B').

    Returns None if the token has no leading integer.
    """
    try:
        return int(resnum_str), ''
    except ValueError:
        pass
    for i, ch in enumerate(resnum_str):
        if not ch.isdigit() and ch != '-':
            digits, icode = resnum_str[:i], resnum_str[i:]
            try:
                return int(digits), icode
            except ValueError:
                return None
    return None


def _res_list_chain_dashes(entry: dict) -> int:
    """Sum of '-' characters in the chain portion of all residue IDs in the entry."""
    return sum(res.split('.')[0].count('-') for res in entry['res_list'])


def _build_seqid_to_aa(pdb_chains: dict) -> dict:
    """Invert cif2seq() output to 'MODEL.CHAIN' -> {(seqnum, icode): one_letter}.

    The key is the composite (seqnum, icode) because 100, 100A and 100B all
    share seqid.num == 100.
    """
    result = {}
    for chain_key, chain in pdb_chains.items():
        result[chain_key] = {
            (seqid, icode): aa
            for aa, seqid, icode in zip(chain.seq, chain.auth_nums, chain.icodes)
        }
    return result


def _parse_cif_cache_entry(pdb_id: str, cif_path: str) -> dict:
    """Sequence and metadata for one CIF assembly.

    Returns a dict with keys seqid_to_aa, resolution, n_chains, total_aa,
    organism; any of them may be None if that part could not be parsed.
    """
    entry = {'seqid_to_aa': None, 'resolution': None,
             'n_chains': None, 'total_aa': None, 'organism': None}
    try:
        chains = cif2seq(pdb_id, cif_path)
        # Prefer model 1 chains (the biologically representative model in assemblies);
        # fall back to all chains if no '1.' prefix exists.
        model1 = {k: v for k, v in chains.items() if k.startswith('1.')} or chains
        entry['seqid_to_aa'] = _build_seqid_to_aa(model1)
        entry['n_chains']    = len(model1)
        entry['total_aa']    = sum(len(c.seq) for c in model1.values())
    except Exception as exc:
        log.warning(f'cif2seq failed {pdb_id}: {exc}')

    # Prefer a non-assembly CIF for metadata (smaller file, same header block).
    meta_cif = os.path.join(os.path.dirname(cif_path), f'{pdb_id}.cif')
    if not os.path.exists(meta_cif):
        meta_cif = cif_path
    try:
        block = gemmi.cif.read_file(meta_cif).sole_block()
        # Try multiple CIF tags for resolution — different deposition workflows use different tags.
        for tag in ('_refine.ls_d_res_high', '_refine_hist.d_res_high',
                    '_reflns.d_resolution_high'):
            val = block.find_value(tag)
            if val and val not in ('.', '?', ''):
                try:
                    entry['resolution'] = float(val)
                    break  # Use the first tag that returns a valid numeric value.
                except ValueError:
                    pass
        # Try recombinant expression source first, then natural source.
        for prefix, col in (('_entity_src_gen.', 'pdbx_gene_src_scientific_name'),
                             ('_entity_src_nat.', 'pdbx_organism_scientific')):
            try:
                for row in block.find(prefix, [col]):
                    val = row[0].strip().strip('"\'')
                    if val and val not in ('.', '?'):
                        entry['organism'] = val
                        break
            except Exception:
                pass
            if entry['organism']:
                break
    except Exception as exc:
        log.warning(f'CIF metadata parse failed {pdb_id}: {exc}')

    return entry



def _site_to_aa_tuple(res_list: list, seqid_to_aa: dict) -> tuple:
    """Residue IDs -> one-letter AA tuple, the canonical form used for clustering.

    Sorted by (seqnum, icode, chain) so a site always yields the same tuple.
    Residues absent from seqid_to_aa are skipped.
    """
    items = []
    for res_id in res_list:
        chain_key, resnum_str = _parse_res_id(res_id)
        parsed = _split_resnum(resnum_str)
        if parsed is None:
            continue  # Malformed residue ID — skip.
        seqid, icode = parsed
        aa = seqid_to_aa.get(chain_key, {}).get((seqid, icode))
        if aa:
            items.append((seqid, icode, chain_key, aa))
    # Sort by numeric position, then insertion code, then chain name (e.g. homodimer).
    items.sort(key=lambda x: (x[0], x[1], x[2]))
    return tuple(aa for _, _, _, aa in items)


def _dominant_chain_unp(res_list: list, pdb_id: str, pdb_chain_to_unp: dict,
                         pdb_chain_to_unp_all: dict, fallback: set) -> str:
    """UniProt accession for the chain contributing the most site residues.

    Resolved against SIFTS best mappings first, then all mappings, then the
    PDB-level set (pdb2unp), so mapping quality outranks the lexicographic
    min() used to break ties within a chain.
    """
    chain_counts: dict = defaultdict(int)
    for res_id in res_list:
        chain_key, _ = _parse_res_id(res_id)
        chain_counts[chain_key.split('.')[1]] += 1
    # Most-contributing chain first, ties broken by name.
    ordered = sorted(chain_counts, key=lambda c: (-chain_counts[c], c))
    for chain_map in (pdb_chain_to_unp, pdb_chain_to_unp_all):
        for chain in ordered:
            unp_set = chain_map.get((pdb_id, chain))
            if unp_set:
                return min(unp_set)
    return min(fallback)


def _cluster(entries: list) -> list:
    """Group entries whose aa_tuples are within EDIT_DIST_THRESHOLD of each other.

    Clustering is transitive (Union-Find): if A~B and B~C then A, B and C share a
    cluster even when A and C are further apart than the threshold.
    """
    from rapidfuzz.distance import Levenshtein
    n = len(entries)
    uf = UnionFind(range(n))
    strs = [''.join(e['aa_tuple']) for e in entries]
    for i in range(n):
        for j in range(i + 1, n):
            if Levenshtein.normalized_distance(strs[i], strs[j]) <= EDIT_DIST_THRESHOLD:
                uf.union(i, j)
    clusters: dict = defaultdict(list)
    for i in range(n):
        clusters[uf.find(i)].append(i)
    return [[entries[i] for i in idx_list] for idx_list in clusters.values()]


def main(mode: str = 'annotated', cif_dir: str = CIF_DIR, verbose: bool = False,
         force: bool = False) -> None:
    """Build the non-redundant binding-site set and write it as a pickle and a CSV.

    Sites are converted to AA tuples, grouped by (UniProt AC, ligand), clustered
    within each group, and one representative per cluster is kept.

    Args:
        mode:  'annotated' — sites with a measured functional-site distance;
               'complete'  — every detected binding site.
        force: Re-parse every CIF instead of reusing pdb_cif_cache.p.
    """
    global log
    os.makedirs(RESULT_NONREDUNDANT_DIR, exist_ok=True)
    log = setup_logging(verbose, os.path.join(LOG_DIR, f'build_nonredundant_sites_{mode}.log'))
    log.info(f'Mode: {mode}')

    distance_results = load_pickle(os.path.join(RESULT_DISTANCES_DIR, 'distance_results.p'))
    site_res_db      = load_pickle(os.path.join(RESULT_BINDING_DIR, 'site_res_db.p'))
    ions             = load_pickle(os.path.join(RESULT_LIGANDS_DIR, 'ions.p'))
    pdb2unp          = load_pickle(os.path.join(RESULT_MAPPING_DIR, 'pdb2unp.p'))
    pdb_chain_to_unp = load_pickle(os.path.join(RESULT_MAPPING_DIR, 'pdb_chain_to_unp.p'))
    # Fallback map covering chains SIFTS gave no best mapping for (see _dominant_chain_unp).
    chain_all_path       = os.path.join(RESULT_MAPPING_DIR, 'pdb_chain_to_unp_all.p')
    pdb_chain_to_unp_all = load_pickle(chain_all_path) if os.path.exists(chain_all_path) else {}
    log.info(f'distance_results={len(distance_results)}  site_res_db={len(site_res_db)}'
             f'  pdb2unp={len(pdb2unp)}  chain_map={len(pdb_chain_to_unp)}')

    # Build / load CIF cache — parsing all CIF files takes ~10 min on first run,
    # so the result is persisted and reused on subsequent runs.
    cif_cache_path = os.path.join(RESULT_NONREDUNDANT_DIR, 'pdb_cif_cache.p')
    if os.path.exists(cif_cache_path) and not force:
        cif_cache: dict = load_pickle(cif_cache_path)
        log.info(f'Loaded CIF cache: {len(cif_cache)} PDBs')
    else:
        if force and os.path.exists(cif_cache_path):
            log.info('Force mode — rebuilding CIF cache from scratch')
        # Only parse PDBs that appear in both site_res_db and have a UniProt mapping.
        needed = {k.split(':')[0] for k in site_res_db} & set(pdb2unp)
        log.info(f'Building CIF cache for {len(needed)} PDBs...')
        print(f"Building CIF cache for {len(needed):,} PDBs (first run)...")
        cif_cache = {}
        for i, pdb_id in enumerate(sorted(needed), 1):
            cif_path = os.path.join(cif_dir, f'{pdb_id}-assembly1.cif')
            if os.path.exists(cif_path):
                cif_cache[pdb_id] = _parse_cif_cache_entry(pdb_id, cif_path)
            else:
                # Placeholder entry keeps the PDB ID in the cache to avoid re-checking later.
                cif_cache[pdb_id] = {'seqid_to_aa': None, 'resolution': None,
                                     'n_chains': None, 'total_aa': None, 'organism': None}
            if i % 5000 == 0:
                print(f"  {i}/{len(needed)}...")
                save_pickle(cif_cache, cif_cache_path)  # Incremental saves guard against crashes.
        save_pickle(cif_cache, cif_cache_path)

    # Determine which ligand keys to include based on the chosen mode.
    if mode == 'annotated':
        # 'annotated' mode: only include sites that have a *computed* distance to a
        # functional site. distance_results stores (None, None, None) for ligands whose
        # binding residues had no extractable coords — exclude those by requiring a
        # non-None min_dist so the set matches the "with a distance" definition.
        lig_key_source = sorted({k for d in distance_results.values()
                                 for k, v in d.items() if v[0] is not None})
    else:
        # 'complete' mode: use every detected binding site, even without a distance.
        lig_key_source = sorted(site_res_db)
    log.info(f'lig_key_source: {len(lig_key_source)}')
    print(f"Building AA-type tuples for {len(lig_key_source):,} ligand keys...")

    # Group valid entries by (UniProt AC, ligand name) ready for clustering.
    groups: dict = defaultdict(list)
    skipped_no_unp = skipped_no_cif = skipped_short = 0

    for lig_key in lig_key_source:
        res_list = site_res_db.get(lig_key)
        if res_list is None:
            continue
        pdb_id, _, _, resname = _parse_ligand_key(lig_key)
        if resname in ions:
            continue
        unp_set = pdb2unp.get(pdb_id, set())
        if not unp_set:
            skipped_no_unp += 1
            continue
        seqid_to_aa = (cif_cache.get(pdb_id) or {}).get('seqid_to_aa')
        if seqid_to_aa is None:
            skipped_no_cif += 1
            continue
        aa_tuple = _site_to_aa_tuple(res_list, seqid_to_aa)
        if len(aa_tuple) < MIN_SITE_RESIDUES:
            skipped_short += 1
            continue
        uniprot_ac = _dominant_chain_unp(res_list, pdb_id, pdb_chain_to_unp,
                                         pdb_chain_to_unp_all, unp_set)
        groups[(uniprot_ac, resname)].append({
            # Strip the model prefix ('1.') from residue IDs — downstream steps don't need it.
            'res_list':   [r.split('.', 1)[1] for r in res_list],
            'aa_tuple':   aa_tuple,
            'source_key': lig_key,
            'n_residues': len(aa_tuple),
        })

    log.info(f'Skipped: no_unp={skipped_no_unp} no_cif={skipped_no_cif} short={skipped_short}')

    # Cluster entries within each (UniProt, ligand) group and pick one representative per cluster.
    print(f"Clustering {sum(len(v) for v in groups.values()):,} entries...")
    nonredundant = []
    for (uniprot_ac, _resname), entries in groups.items():
        for cluster in _cluster(entries):
            best_n = max(e['n_residues'] for e in cluster)
            # Tiebreak order:
            #   1. most residues
            #   2. fewest '-' in binding-residue chain IDs (prefer canonical chains over assembly chains)
            #   3. lexicographically earliest PDB ID (deterministic)
            rep = min(
                (e for e in cluster if e['n_residues'] == best_n),
                key=lambda e: (_res_list_chain_dashes(e), e['source_key'].split(':')[0]),
            )
            nonredundant.append({
                'uniprot_ac':    uniprot_ac,
                'source_pdb':    rep['source_key'].split(':')[0],
                'ligand_code':   _ligand_code(rep['source_key']),
                'site_res_list': rep['res_list'],
            })

    out_p   = os.path.join(RESULT_NONREDUNDANT_DIR, f'_nonredundant_sites_{mode}.p')
    out_csv = os.path.join(RESULT_NONREDUNDANT_DIR, f'_nonredundant_sites_{mode}.csv')
    save_pickle(nonredundant, out_p)
    with open(out_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['uniprot_ac', 'source_pdb', 'ligand_code', 'site_res_list'])
        writer.writeheader()
        for row in nonredundant:
            # site_res_list is a Python list; join with ';' for a single CSV cell.
            writer.writerow({**row, 'site_res_list': ';'.join(row['site_res_list'])})

    log.info(f'Saved {len(nonredundant)} non-redundant records → {out_p}')
    print(f"Done. {len(nonredundant):,} non-redundant records → {RESULT_NONREDUNDANT_DIR}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['annotated', 'complete'], default='annotated')
    parser.add_argument('--cif-dir', default=CIF_DIR)
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-f', '--force', action='store_true',
                        help='Rebuild pdb_cif_cache.p instead of reusing it')
    args = parser.parse_args()
    main(args.mode, args.cif_dir, args.verbose, args.force)
