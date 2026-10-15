"""
Enrich nonredundant binding sites with distance, structure, and metadata.

Usage:
    python enrich_nonredundant_sites.py --input PATH [-v]

Adds columns: site_category, min_dist, median_dist, avg_dist, pdb_resolution,
n_chains, total_aa, organism, ec_number, ec_source, protein_name_uniprot, protein_name_pdb,
chains_in_site, nearby_ions, cath_superfamily, functional_site_res.
"""

import argparse
import csv as _csv
import os
import sys

import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    RESULT_MAPPING_DIR, RESULT_BINDING_DIR, RESULT_FUNCTIONAL_DIR,
    RESULT_DISTANCES_DIR, RESULT_NONREDUNDANT_DIR, RESULT_PROTEIN_NAMES_DIR, LOG_DIR,
    RCSB_META_WORKERS, RCSB_META_TIMEOUT, ORTHOSTERIC_MAX_DIST,
)
from utils import canonical_accession, create_session, load_pickle, save_pickle, setup_logging

log = None  # set in main()

RCSB_ENTRY_URL  = 'https://data.rcsb.org/rest/v1/core/entry/{}'
RCSB_ENTITY_URL = 'https://data.rcsb.org/rest/v1/core/polymer_entity/{}/1'

FIELDNAMES = [
    'ID', 'source_pdb', 'ligand_code', 'binding_site', 'functional_site_res',
    'site_category', 'min_dist', 'median_dist', 'avg_dist', 'nearby_ions', 'uniprot_ac',
    'protein_name_uniprot', 'protein_name_pdb', 'chains_in_site', 'organism',
    'pdb_resolution', 'n_chains', 'total_aa', 'cath_superfamily', 'ec_number',
    'ec_source',
]


def _build_id(source_pdb: str, ligand_code: str) -> str:
    return f"{source_pdb}_{ligand_code.replace('.', '').replace(':', '')}"


def _assign_site_category(min_dist: float | None) -> str | None:
    """Classify a site by min_dist. Returns None when min_dist is unknown —
    a bare `min_dist <= threshold` would evaluate None/NaN as 'allosteric'.
    """
    if min_dist is None:
        return None
    return 'orthosteric' if min_dist <= ORTHOSTERIC_MAX_DIST else 'allosteric'


# ── RCSB metadata fetch ───────────────────────────────────────────────────────

def _fetch_pdb_meta(pdb_id: str, session) -> tuple:
    """(resolution, organism) from the RCSB REST API; either may be None.

    Used only for entries where the local CIF gave no value. Failures are
    ignored so one bad response does not block the rest.
    """
    resolution, organism = None, None
    try:
        r = session.get(RCSB_ENTRY_URL.format(pdb_id), timeout=RCSB_META_TIMEOUT)
        if r.ok:
            # resolution_combined is a list; take the first (best) value.
            vals = r.json().get('rcsb_entry_info', {}).get('resolution_combined') or []
            if vals:
                resolution = float(vals[0])
    except Exception:
        pass
    try:
        r = session.get(RCSB_ENTITY_URL.format(pdb_id), timeout=RCSB_META_TIMEOUT)
        if r.ok:
            orgs = r.json().get('rcsb_entity_source_organism') or []
            if orgs:
                organism = orgs[0].get('scientific_name')
    except Exception:
        pass
    return resolution, organism


# ── Distance lookup ───────────────────────────────────────────────────────────

def _build_dist_lookup(dist_results: dict) -> dict:
    """Flatten distance_results.p to {(pdb_id, chain, resnum, resname): (min, median, avg)}.

    The model prefix is dropped so the keys match the 'A.101:ATP' ligand codes
    used in the nonredundant records.
    """
    lookup = {}
    for pdb_id, lig_dict in dist_results.items():
        if not isinstance(lig_dict, dict):
            continue
        for lig_key, dists in lig_dict.items():
            parts = lig_key.split(':')
            if len(parts) != 4:
                continue
            # parts[1] is 'MODEL.CHAIN' (e.g. '1.A'); strip the model prefix.
            chain = parts[1].split('.')[1] if '.' in parts[1] else parts[1]
            lookup[(pdb_id, chain, parts[2], parts[3])] = dists
    return lookup


def _lookup_dists(lookup: dict, pdb_id: str, ligand_code: str) -> tuple:
    """Return (min_dist, median_dist, avg_dist) or (None, None, None)."""
    try:
        chain_resnum, resname = ligand_code.split(':')
        chain, resnum = chain_resnum.split('.', 1)
    except ValueError:
        return None, None, None
    return lookup.get((pdb_id, chain, resnum, resname), (None, None, None))


# ── EC lookup (from structured UniProt ec field) ──────────────────────────────

def _build_chain_range_lookup(uniprot_csv: str) -> dict:
    """(pdb_id, asym_id, unp_acc) -> (span_start, span_end, [(start, end), ...]).

    The span is the union of the chain's SIFTS segments; the segments themselves
    are kept so _narrowed_ec can tell an internal gap from covered sequence.
    Only best_mapping == 'y' rows are used.
    """
    segs: dict = {}
    with open(uniprot_csv) as f:
        for row in _csv.DictReader(f):
            if row.get('best_mapping') != 'y':
                continue
            key = (row['pdb_id'].upper(), row['asym_id'], row['unp_acc'])
            segs.setdefault(key, []).append((int(row['unp_start']), int(row['unp_end'])))
    return {
        k: (min(s for s, _ in v), max(e for _, e in v), v)
        for k, v in segs.items()
    }


_EC_SOURCES = ['name', 'catalytic_activity', 'brenda']


def _narrowed_ec(unp: str, chain: str, pdb_id: str, ec_lookup: dict,
                  ec_segments_lookup: dict, chain_range_lookup: dict) -> tuple[list[str], str]:
    """(EC numbers, ec_source) for one accession's contribution to a site.

    For a cleaved polyprotein the pooled ec_lookup value unions every cleavage
    product's EC, which is wrong for any single structure. Where per-product
    segments and the chain's UniProt range are both known, this narrows to the
    segment containing the chain midpoint; overlapping segments resolve to the
    narrowest, so a site on a fully-cleaved domain gets that domain's own EC
    rather than its uncleaved precursor's.

    Falls back to the pooled value whenever narrowing is not possible, so it
    never returns less than ec_lookup alone. Isoforms use their canonical entry's EC.
    """
    segments = ec_segments_lookup.get(unp)
    rng = chain_range_lookup.get((pdb_id, chain, unp))
    if segments and rng:
        span_start, span_end, chain_segs = rng
        midpoint = (span_start + span_end) / 2
        # Only narrow when the midpoint falls in sequence the chain really covers.
        if any(s <= midpoint <= e for s, e in chain_segs):
            matches = [(start, end, ecs) for start, end, ecs in segments
                       if start <= midpoint <= end]
            if matches:
                # Overlapping segments: the narrowest (most specific) product wins.
                _, _, ecs = min(matches, key=lambda seg: seg[1] - seg[0])
                return list(ecs), 'name'
    return ec_lookup.get(canonical_accession(unp), ([], ''))


def _ec_for_site(site_res: list, pdb_id: str, pdb_chain_to_unp: dict, ec_lookup: dict,
                  ec_segments_lookup: dict, chain_range_lookup: dict,
                  fallback_uniprot_ac: str | None = None):
    """Semicolon-separated EC numbers for every accession contributing to a site.

    A site can span several chains and each chain map to several accessions.
    pdb_chain_to_unp is keyed by SIFTS asym_id, which lacks the NCS suffixes
    ('A-2') that site_res_list carries, so each chain is looked up under its
    exact label and then its canonical form. If every chain misses, falls back
    to fallback_uniprot_ac. Returns None when nothing is found.

    Args:
        site_res:         List of residue IDs in 'CHAIN.RESNUM' format (model prefix stripped).
        pdb_chain_to_unp: Dict mapping (pdb_id, chain) → set of UniProt ACs.
        ec_lookup:        Dict mapping canonical UniProt AC → ([pooled EC, ...], ec_source).
        ec_segments_lookup: Dict mapping UniProt AC → list of (unp_start, unp_end, [ec, ...])
            cleavage-product segments, for polyprotein accessions only (see _narrowed_ec).
        chain_range_lookup: Dict mapping (pdb_id, asym_id, unp_acc) →
            (unp_start, unp_end, [(seg_start, seg_end), ...]) from mapped_uniprot.csv —
            the chain's union span plus its individual SIFTS segments.
        fallback_uniprot_ac: The site's already-resolved dominant UniProt AC, used when
            no chain in site_res maps to anything in pdb_chain_to_unp.

    Returns:
        (EC numbers, ec_source), each ';'-separated, or (None, None) if no EC found.
    """
    # Extract unique chain letters from the residue list.
    chains = {r.split('.')[0] for r in site_res}
    seen: set = set()
    ec_numbers: list = []
    sources: list = []
    for chain in sorted(chains):
        # Try the exact chain label, then the canonical suffix-stripped one.
        canonical = chain.split('-')[0]
        mapped = (pdb_chain_to_unp.get((pdb_id, chain))
                  or pdb_chain_to_unp.get((pdb_id, canonical))
                  or [])
        for unp in sorted(mapped):
            if unp in seen:
                continue
            seen.add(unp)
            # Range lookups are keyed by the canonical asym_id for the same reason.
            ecs, source = _narrowed_ec(unp, canonical, pdb_id, ec_lookup,
                                       ec_segments_lookup, chain_range_lookup)
            for e in ecs:
                if e not in ec_numbers:
                    ec_numbers.append(e)
            if ecs and source not in sources:
                sources.append(source)
    if not ec_numbers and fallback_uniprot_ac:
        ec_numbers, source = ec_lookup.get(canonical_accession(fallback_uniprot_ac), ([], ''))
        sources = [source] if ec_numbers else []
    if not ec_numbers:
        return None, None
    return ';'.join(ec_numbers), ';'.join(sorted(sources, key=_EC_SOURCES.index))


# ── Ion proximity lookup ──────────────────────────────────────────────────────

def _lookup_ions(site_ion_db: dict, pdb_id: str, ligand_code: str) -> list[str]:
    """Return list of nearby ion strings for this ligand, e.g. ['A.100:ZN']."""
    try:
        chain_resnum, resname = ligand_code.split(':')
        chain, resnum = chain_resnum.split('.', 1)
    except ValueError:
        return []
    key = f"{pdb_id}:1.{chain}:{resnum}:{resname}"
    return site_ion_db.get(key, [])


# ── CATH lookup ───────────────────────────────────────────────────────────────

def _build_cath_lookup(cath_csv: str) -> dict:
    """Build (pdb_id, asym_id) → [superfamily, ...] from mapped_cath.csv."""
    lookup: dict = {}
    if not os.path.exists(cath_csv):
        return lookup
    with open(cath_csv) as f:
        for row in _csv.DictReader(f):
            key = (row['pdb_id'].upper(), row['asym_id'])
            sf = row['xref_db_acc'].strip()
            if sf and sf not in lookup.get(key, []):
                lookup.setdefault(key, []).append(sf)
    return lookup


def _build_protein_names_lookup(names_csv: str) -> dict:
    """(pdb_id, ligand_code) -> (protein_name_pdb, chains_in_site); {} if absent."""
    lookup: dict = {}
    if not os.path.exists(names_csv):
        return lookup
    with open(names_csv) as f:
        for row in _csv.DictReader(f):
            key = (row['source_pdb'].upper(), row['ligand_code'])
            lookup[key] = (row.get('protein_name_pdb') or None, row.get('chains_in_site') or None)
    return lookup


def _cath_for_site(site_res_list: list, pdb_id: str, cath_lookup: dict) -> str | None:
    """CATH superfamilies for the chain contributing the most residues, or None."""
    chain_counts: dict = {}
    for res in site_res_list:
        chain = res.split('.')[0]
        chain_counts[chain] = chain_counts.get(chain, 0) + 1
    if not chain_counts:
        return None
    # Pick the chain with the highest residue count.
    dominant = max(chain_counts, key=chain_counts.get)
    # cath_lookup is keyed by SIFTS asym_id, which never carries the NCS/symmetry
    # suffixes (e.g. 'A-2') that appear in site_res_list for assembly-expanded chains.
    # Try the exact label first, then fall back to the canonical (suffix-stripped) chain.
    sfs = cath_lookup.get((pdb_id, dominant)) or cath_lookup.get((pdb_id, dominant.split('-')[0]))
    return ';'.join(sfs) if sfs else None


# ── Main ──────────────────────────────────────────────────────────────────────

def main(input_path: str, verbose: bool = False) -> None:
    """Join a nonredundant site pickle with every auxiliary source and write the result.

    Columns added, and where each comes from:
      site_category, min_dist, median_dist, avg_dist   distance_results.p
      pdb_resolution, organism, n_chains, total_aa     pdb_cif_cache.p
                                                       (resolution/organism gaps
                                                       filled from the RCSB API here)
      protein_name_uniprot, ec_number, ec_source       uniprot_names.p, narrowed
                                                       per cleavage product where
                                                       segments exist
      protein_name_pdb, chains_in_site                 protein_names_pdb.csv
      nearby_ions                                      site_ion_db.p
      cath_superfamily                                 mapped_cath.csv
      functional_site_res                              func_site_dict.p
    """
    global log
    os.makedirs(RESULT_NONREDUNDANT_DIR, exist_ok=True)
    mode    = 'annotated' if 'annotated' in os.path.basename(input_path) else 'complete'
    out_p   = os.path.join(RESULT_NONREDUNDANT_DIR, f'nonredundant_sites_{mode}_enriched.p')
    out_csv = os.path.join(RESULT_NONREDUNDANT_DIR, f'nonredundant_sites_{mode}_enriched.csv')
    log = setup_logging(verbose, os.path.join(LOG_DIR, 'enrich_nonredundant_sites.log'))

    nonredundant     = load_pickle(input_path)
    dist_results     = load_pickle(os.path.join(RESULT_DISTANCES_DIR, 'distance_results.p'))
    uniprot_names    = load_pickle(os.path.join(RESULT_MAPPING_DIR, 'uniprot_names.p'))
    pdb_chain_to_unp = load_pickle(os.path.join(RESULT_MAPPING_DIR, 'pdb_chain_to_unp.p'))
    ec_segments_path = os.path.join(RESULT_MAPPING_DIR, 'uniprot_ec_segments.p')
    ec_segments_lookup  = load_pickle(ec_segments_path) if os.path.exists(ec_segments_path) else {}
    chain_range_lookup = _build_chain_range_lookup(os.path.join(RESULT_MAPPING_DIR, 'mapped_uniprot.csv'))
    func_site_dict   = load_pickle(os.path.join(RESULT_FUNCTIONAL_DIR, 'func_site_dict.p'))
    cif_cache_path   = os.path.join(RESULT_NONREDUNDANT_DIR, 'pdb_cif_cache.p')
    cif_cache        = load_pickle(cif_cache_path)
    site_ion_db      = load_pickle(os.path.join(RESULT_BINDING_DIR, 'site_ion_db.p'))
    log.info(f'nonredundant={len(nonredundant)} dist_results={len(dist_results)}'
             f' uniprot_names={len(uniprot_names)} cif_cache={len(cif_cache)}'
             f' site_ion_db={len(site_ion_db)}')

    # Back-fill resolution/organism for representative PDBs missing either field.
    # Only representative PDB IDs are needed — enrich never reads non-representative entries.
    rep_pdbs = {rec['source_pdb'] for rec in nonredundant}
    need_meta = [
        p for p in rep_pdbs
        if cif_cache.get(p) and cif_cache[p].get('seqid_to_aa') is not None
        and (cif_cache[p].get('resolution') is None or cif_cache[p].get('organism') is None)
    ]
    if need_meta:
        log.info(f'Fetching RCSB metadata for {len(need_meta)} representative PDBs...')
        print(f"Fetching RCSB metadata for {len(need_meta):,} representative PDBs...")
        session = create_session()
        done = 0
        with ThreadPoolExecutor(max_workers=RCSB_META_WORKERS) as executor:
            futures = {executor.submit(_fetch_pdb_meta, p, session): p for p in need_meta}
            for future in as_completed(futures):
                p = futures[future]
                try:
                    res, org = future.result()
                    # Only fill gaps; never overwrite a CIF-parsed value.
                    if res is not None:
                        cif_cache[p]['resolution'] = res
                    if org is not None:
                        cif_cache[p]['organism'] = org
                except Exception:
                    pass
                done += 1
                if done % 1000 == 0:
                    print(f"  {done}/{len(need_meta)}...")
        save_pickle(cif_cache, cif_cache_path)
        log.info(f'Updated cif_cache with RCSB metadata → {cif_cache_path}')

    # Pre-build all lookup structures so the main loop below is a simple dict lookup.
    dist_lookup  = _build_dist_lookup(dist_results)
    cath_lookup  = _build_cath_lookup(os.path.join(RESULT_MAPPING_DIR, 'mapped_cath.csv'))
    names_lookup = _build_protein_names_lookup(os.path.join(RESULT_PROTEIN_NAMES_DIR, 'protein_names_pdb.csv'))

    # uniprot_names is a pandas DataFrame; convert to plain dicts for fast per-row access.
    name_lookup: dict = {}
    ec_lookup:   dict = {}
    for _, r in uniprot_names.iterrows():
        acc = r['unp_acc']
        name_lookup[acc] = r.get('protein_name')
        val = r.get('ec_number')
        ecs = [e.strip() for e in val.split(';') if e.strip()] if pd.notna(val) else []
        ec_lookup[acc] = (ecs, r.get('ec_source'))

    # For each PDB, collect the set of functional residue IDs (strip model prefix '1.').
    # e[1] is the residue ID; split on '.' once and take the suffix (CHAIN.RESNUM).
    func_lookup = {pdb: sorted({e[1].split('.', 1)[1] for e in entries if len(e) >= 2})
                   for pdb, entries in func_site_dict.items()}

    print(f"Enriching {len(nonredundant):,} records...")
    enriched = []
    for rec in nonredundant:
        pdb_id     = rec['source_pdb']
        lig_code   = rec['ligand_code']
        uniprot_ac = rec['uniprot_ac']
        meta       = cif_cache.get(pdb_id, {})
        min_d, med_d, avg_d = _lookup_dists(dist_lookup, pdb_id, lig_code)
        func_res   = func_lookup.get(pdb_id)
        ions       = _lookup_ions(site_ion_db, pdb_id, lig_code)
        name_pdb, chains_in_site = names_lookup.get((pdb_id, lig_code), (None, None))
        ec_number, ec_source = _ec_for_site(rec['site_res_list'], pdb_id, pdb_chain_to_unp,
                                            ec_lookup, ec_segments_lookup, chain_range_lookup,
                                            fallback_uniprot_ac=uniprot_ac)
        enriched.append({
            'ID':                  _build_id(pdb_id, lig_code),
            'source_pdb':          pdb_id,
            'ligand_code':         lig_code,
            'binding_site':        rec['site_res_list'],
            'functional_site_res': func_res,
            'site_category':       _assign_site_category(min_d),
            'min_dist':            min_d,
            'median_dist':         med_d,
            'avg_dist':            avg_d,
            'nearby_ions':         ions,
            'uniprot_ac':          uniprot_ac,
            'protein_name_uniprot': name_lookup.get(canonical_accession(uniprot_ac)) if uniprot_ac else None,
            'protein_name_pdb':    name_pdb,
            'chains_in_site':      chains_in_site,
            # Upper-case only the first character so lowercase CIF-sourced names become
            # 'Homo sapiens' without corrupting already-cased RCSB names like
            # 'Influenza A virus' or 'Escherichia coli (strain K12)'.
            'organism':            (meta['organism'][:1].upper() + meta['organism'][1:]) if meta.get('organism') else None,
            'pdb_resolution':      meta.get('resolution'),
            'n_chains':            meta.get('n_chains'),
            'total_aa':            meta.get('total_aa'),
            'cath_superfamily':    _cath_for_site(rec['site_res_list'], pdb_id, cath_lookup),
            'ec_number':           ec_number,
            'ec_source':           ec_source,
        })

    save_pickle(enriched, out_p)
    with open(out_csv, 'w', newline='') as f:
        writer = _csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        for row in enriched:
            out = dict(row)
            # List fields must be serialised as semicolon-joined strings for CSV output.
            out['binding_site']        = ';'.join(row['binding_site'])
            out['functional_site_res'] = ';'.join(row['functional_site_res']) if row['functional_site_res'] else ''
            out['nearby_ions']         = ';'.join(row['nearby_ions']) if row['nearby_ions'] else ''
            out['cath_superfamily']    = row['cath_superfamily'] or ''
            writer.writerow(out)

    n_dist = sum(1 for r in enriched if r['min_dist'] is not None)
    log.info(f'Enriched {len(enriched)} records; {n_dist} with distance data → {out_p}')
    print(f"Done. {len(enriched):,} records ({n_dist} with distance) → {out_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()
    main(args.input, args.verbose)
