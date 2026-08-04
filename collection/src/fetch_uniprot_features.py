"""
Fetch UniProt functional site annotations (Active site, Binding site, DNA binding)
for all proteins linked to PDB entries with ligands.

Automatically retries failed accessions with a longer timeout if any errors occur.

Usage:
    python fetch_uniprot_features.py [-v]

Inputs:  result/binding_sites/pdb_with_ligand.p, result/mapping/pdb2unp.p
Outputs: result/functional/unp2feature.p, unp2site.p, seq_dict.p
"""

import argparse
import os
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    RESULT_MAPPING_DIR, RESULT_BINDING_DIR, RESULT_FUNCTIONAL_DIR,
    LOG_FUNCTIONAL_DIR, UNIPROT_FETCH_TIMEOUT, UNIPROT_RETRY_TIMEOUT,
)
from utils import create_session, load_pickle, save_pickle, setup_logging

log = None  # set in main()


def _fetch_json(protein_id: str, session, timeout: int = 5) -> dict:
    """JSON record for one UniProt accession. Raises HTTPError on a non-2xx status."""
    r = session.get(f"https://rest.uniprot.org/uniprotkb/{protein_id}.json", timeout=timeout)
    r.raise_for_status()
    return r.json()


def _extract_sites(data: dict, uniprot_id: str) -> tuple[list[dict], str]:
    """Parse a UniProt JSON record to extract functional site annotations and the sequence.

    Extracts three feature types: 'Active site', 'Binding site', 'DNA binding'.

    A feature yields one row *per residue it spans*, not one per feature:
    'DNA binding' is a domain range. 'Active site' is single-residue.

    A protein with no relevant features still gets one all-None row so it
    appears in the DataFrame.

    Returns (rows, sequence).
    """
    rec_name = (data.get('proteinDescription', {})
                    .get('recommendedName', {})
                    .get('fullName', {})
                    .get('value', ''))
    organism = data.get('organism', {}).get('scientificName', '')
    sequence = data.get('sequence', {}).get('value', '')

    target_types = {'Active site', 'Binding site', 'DNA binding'}
    rows = []
    for feat in data.get('features', []):
        ftype = feat.get('type')
        if ftype not in target_types:
            continue
        loc = feat.get('location', {})
        start_pos = loc.get('start', {}).get('value')
        end_pos   = loc.get('end', {}).get('value')
        if start_pos is None:
            continue
        # A feature with no end, or an end before its start, is a single residue.
        if end_pos is None or end_pos < start_pos:
            end_pos = start_pos
        # Clamp to the sequence so a stale annotation can't emit phantom residues.
        last = len(sequence) if sequence else end_pos
        for pos in range(start_pos, min(end_pos, last) + 1):
            # Look up the one-letter residue at this position (1-based → 0-based index).
            residue = sequence[pos - 1] if sequence and 0 < pos <= len(sequence) else ''
            rows.append({
                'uniprot_id':      uniprot_id,
                'protein_name':    rec_name,
                'organism':        organism,
                'sequence_length': len(sequence) or None,
                'feature_type':    ftype,
                'position':        pos,
                # The parent feature's end, repeated on every residue of the span.
                'end_position':    end_pos,
                'description':     feat.get('description', ''),
                'residue':         residue,
            })

    if not rows:
        # Placeholder row ensures the protein is present in the DataFrame even without sites.
        rows.append({
            'uniprot_id': uniprot_id, 'protein_name': rec_name, 'organism': organism,
            'sequence_length': len(sequence) or None,
            'feature_type': None, 'position': None, 'end_position': None,
            'description': None, 'residue': None,
        })
    return rows, sequence


def _fetch_batch(protein_ids: list[str], max_workers: int = 15,
                 timeout: int = UNIPROT_FETCH_TIMEOUT) -> tuple[pd.DataFrame, dict[str, str], list[str]]:
    """Fetch many accessions in parallel.

    Returns (feature-row DataFrame, {accession: sequence}, list of failed IDs).
    """
    all_rows, seq_dict, errors = [], {}, []
    total = len(protein_ids)
    completed = 0

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_fetch_one, pid, timeout): pid for pid in protein_ids}
        for future in as_completed(futures):
            pid = futures[future]
            completed += 1
            result = future.result()
            if result['success']:
                all_rows.extend(result['rows'])
                seq_dict[pid] = result['sequence']
            else:
                errors.append(pid)
                log.debug(f"Failed: {pid}: {result['error']}")
            if completed % 100 == 0 or completed == total:
                print(f"  [{completed}/{total}] fetched")

    return pd.DataFrame(all_rows), seq_dict, errors


def _fetch_one(protein_id: str, timeout: int = UNIPROT_FETCH_TIMEOUT) -> dict:
    """Fetch and parse one accession.

    Returns {'success': True, 'rows': [...], 'sequence': '...'} or
    {'success': False, 'error': '...'} — never raises, so one failure does not
    cancel the batch. Each call makes its own Session (Session is not
    thread-safe).
    """
    try:
        session = create_session()
        data = _fetch_json(protein_id, session, timeout)
        rows, sequence = _extract_sites(data, protein_id)
        return {'success': True, 'protein_id': protein_id, 'rows': rows, 'sequence': sequence}
    except Exception as e:
        return {'success': False, 'protein_id': protein_id, 'error': str(e)}


def main(verbose: bool = False, force: bool = False) -> None:
    """Fetch UniProt features for every protein linked to a PDB entry with a ligand.

    Accessions already in seq_dict.p are skipped; failures are retried once with
    a longer timeout. unp2site.p is the subset of rows that are 'Active site' or
    'DNA binding' — the categories the distances step measures against.

    force re-fetches everything. It is required whenever _extract_sites()
    changes, since unp2site.p is otherwise rebuilt from the cached
    unp2feature.p.
    """
    global log
    os.makedirs(RESULT_FUNCTIONAL_DIR, exist_ok=True)
    os.makedirs(LOG_FUNCTIONAL_DIR, exist_ok=True)
    log = setup_logging(verbose, os.path.join(LOG_FUNCTIONAL_DIR, 'fetch_uniprot_features.log'))

    pdb_with_ligand = load_pickle(os.path.join(RESULT_BINDING_DIR, 'pdb_with_ligand.p'))
    pdb2unp         = load_pickle(os.path.join(RESULT_MAPPING_DIR, 'pdb2unp.p'))
    # Flatten: collect every UniProt AC that appears in any PDB with a ligand.
    all_unp = {unp for pdb in pdb_with_ligand for unp in pdb2unp.get(pdb, [])}
    log.info(f"UniProt IDs needed: {len(all_unp):,}")

    # Load existing cache — seq_dict is the authoritative cache key because it is
    # populated only when a full successful fetch has been completed.
    cache_feature = os.path.join(RESULT_FUNCTIONAL_DIR, 'unp2feature.p')
    cache_seq     = os.path.join(RESULT_FUNCTIONAL_DIR, 'seq_dict.p')
    if not force and os.path.exists(cache_feature) and os.path.exists(cache_seq):
        unp2feature = load_pickle(cache_feature)
        seq_dict    = load_pickle(cache_seq)
        to_fetch    = list(all_unp - set(seq_dict))  # Only fetch IDs not yet in cache.
        log.info(f"Cache: {len(seq_dict):,} entries. Fetching {len(to_fetch):,} new.")
    else:
        unp2feature, seq_dict = pd.DataFrame(), {}
        to_fetch = list(all_unp)
        if force:
            log.info(f"Force mode — ignoring cache, re-fetching {len(to_fetch):,} entries.")
            print(f"Force mode — re-fetching all {len(to_fetch):,} UniProt records "
                  f"(the cached ones were parsed by an older _extract_sites).")
        else:
            log.info(f"No cache. Fetching {len(to_fetch):,} entries.")

    if to_fetch:
        print(f"Fetching UniProt features for {len(to_fetch):,} proteins...")
        new_df, new_seq, errors = _fetch_batch(to_fetch)
        unp2feature = pd.concat([unp2feature, new_df], ignore_index=True)
        seq_dict.update(new_seq)

        if errors:
            error_path = os.path.join(LOG_FUNCTIONAL_DIR, 'error.txt')
            with open(error_path, 'w') as f:
                f.writelines(e + '\n' for e in errors)
            log.warning(f"{len(errors)} failed accessions written to {error_path}")

            # Auto-retry with longer timeout — first-attempt failures are often transient.
            log.info(f"Retrying {len(errors)} failed accessions (timeout={UNIPROT_RETRY_TIMEOUT}s)...")
            print(f"Retrying {len(errors)} failed accessions...")
            retry_df, retry_seq, still_failed = _fetch_batch(errors, timeout=UNIPROT_RETRY_TIMEOUT)
            if not retry_df.empty:
                unp2feature = pd.concat([unp2feature, retry_df], ignore_index=True)
                seq_dict.update(retry_seq)
                log.info(f"Retry recovered {len(retry_seq)}/{len(errors)} accessions")
            # Always rewrite error.txt so a fully successful retry truncates it —
            # otherwise the file keeps listing accessions that were since recovered.
            with open(error_path, 'w') as f:
                f.writelines(e + '\n' for e in still_failed)
            if still_failed:
                log.warning(f"{len(still_failed)} still failed after retry")
            else:
                log.info("Retry recovered every failed accession")

    # unp2site is the subset of features used by map_functional_sites.py:
    # only active sites and DNA binding sites with a known residue position.
    unp2site = pd.concat([
        unp2feature[(unp2feature.feature_type == 'Active site') & unp2feature.position.notna()],
        unp2feature[(unp2feature.feature_type == 'DNA binding') & unp2feature.position.notna()],
    ], ignore_index=True)

    save_pickle(unp2feature, os.path.join(RESULT_FUNCTIONAL_DIR, 'unp2feature.p'))
    save_pickle(unp2site,    os.path.join(RESULT_FUNCTIONAL_DIR, 'unp2site.p'))
    save_pickle(seq_dict,    os.path.join(RESULT_FUNCTIONAL_DIR, 'seq_dict.p'))
    log.info(f"Saved: unp2feature({len(unp2feature)}), unp2site({len(unp2site)}), seq_dict({len(seq_dict)})")
    print(f"Done. unp2site: {len(unp2site):,} rows  |  seq_dict: {len(seq_dict):,} proteins")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-f", "--force", action="store_true",
                        help="Ignore unp2feature.p/seq_dict.p and re-fetch everything")
    args = parser.parse_args()
    main(args.verbose, args.force)
