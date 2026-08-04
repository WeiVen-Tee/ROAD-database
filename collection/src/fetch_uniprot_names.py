"""
Fetch UniProt protein names for all unique UniProt accessions mapped to PDB entries.

Reads pdb2unp.p to derive the set of canonical accessions, then fetches protein
names in parallel batches via the UniProt search API.

For accessions with more than one pooled EC number (polyproteins/multi-catalytic-
domain proteins, e.g. viral replicase polyproteins), the pooled `ec` field from the
search API conflates every cleavage product's EC number with no positional
information. A second pass fetches the full JSON record for just these accessions
and extracts per-cleavage-product residue segments from `proteinDescription.contains`
(physically cleaved mature chains) joined to `features` of type 'Chain' by name.
`proteinDescription.includes` (bifunctional single-chain enzymes, e.g. PFKFB3, where
every activity is present in every structure) is intentionally left untouched — those
accessions have no matching Chain feature to narrow against, so their EC numbers stay
pooled, which is correct for that case.

Usage:
    python fetch_uniprot_names.py [-v]

Outputs (result/mapping/):
    uniprot_names.p        — DataFrame: unp_acc, protein_name, ec_number
    uniprot_names.csv      — same as CSV
    uniprot_ec_segments.p  — dict[unp_acc] -> list[(unp_start, unp_end, [ec_number, ...])],
                              populated only for cleaved multi-EC accessions where every
                              EC-bearing cleavage product could be matched to a residue range
    uniprot_ec_segments.csv — same as CSV (unp_acc, unp_start, unp_end, ec_numbers)
"""

import argparse
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import RESULT_MAPPING_DIR, UNIPROT_SEARCH_URL, UNIPROT_CHUNK_SIZE, UNIPROT_WORKERS, LOG_DIR
from utils import create_session, load_pickle, save_pickle, setup_logging

log = None  # set in main()

_ISOFORM_RE = re.compile(r"-\d+$")


def _fetch_chunk(accessions: list[str], session) -> list[dict]:
    """Names and EC numbers for a batch of accessions in one search-API call.

    Returns [{unp_acc, protein_name, ec_number}, ...] parsed from the TSV response.
    """
    query = "accession:(" + " OR ".join(accessions) + ")"
    r = session.get(UNIPROT_SEARCH_URL, params={
        "query": query,
        "fields": "accession,protein_name,ec",
        "format": "tsv",
        "size": len(accessions),
    }, timeout=60)
    r.raise_for_status()
    rows = []
    for line in r.text.splitlines()[1:]:  # Skip TSV header row.
        parts = line.split("\t")
        if len(parts) >= 2:
            rows.append({
                "unp_acc":      parts[0],
                "protein_name": parts[1],
                "ec_number":    parts[2].strip() if len(parts) >= 3 else "",
            })
    return rows


# ── EC segment disambiguation (for cleaved polyproteins) ──────────────────────

def _extract_ec_segments(data: dict) -> list[tuple[int, int, list[str]]] | None:
    """Per-cleavage-product EC segments from a UniProt JSON entry.

    Uses only `proteinDescription.contains` — chains the polyprotein is
    physically cleaved into, so one structure matches at most one. `includes`
    (bifunctional single-chain enzymes) is skipped: those have no residue range
    of their own, and pooling is correct for them.

    Returns [(unp_start, unp_end, [ec, ...]), ...], or None if any EC-bearing
    entry cannot be matched to a Chain feature — a partial set would silently
    drop EC numbers instead of falling back to the pooled list.
    """
    chain_ranges: dict[str, tuple[int, int]] = {}
    for feat in data.get('features', []):
        if feat.get('type') != 'Chain':
            continue
        desc = feat.get('description')
        loc = feat.get('location', {})
        start = loc.get('start', {}).get('value')
        end = loc.get('end', {}).get('value')
        if desc and start is not None and end is not None:
            chain_ranges[desc] = (start, end)

    segments = []
    for entry in data.get('proteinDescription', {}).get('contains', []):
        rec = entry.get('recommendedName', {})
        ec_numbers = [e['value'] for e in rec.get('ecNumbers', [])]
        if not ec_numbers:
            continue
        name = rec.get('fullName', {}).get('value')
        rng = chain_ranges.get(name)
        if rng is None:
            return None
        segments.append((rng[0], rng[1], ec_numbers))
    return segments or None


def _fetch_ec_segments_one(acc: str) -> list | None:
    """EC segments for one accession, or None.

    Called only for accessions whose pooled ec_number holds more than one value.
    """
    session = create_session()
    r = session.get(f"https://rest.uniprot.org/uniprotkb/{acc}.json", timeout=30)
    r.raise_for_status()
    return _extract_ec_segments(r.json())


def _build_ec_segments(df: 'pd.DataFrame') -> dict:
    """{unp_acc: [(unp_start, unp_end, [ec, ...]), ...]} for every multi-EC accession."""
    ambiguous = sorted(df.loc[df['ec_number'].str.contains(';', na=False), 'unp_acc'])
    if not ambiguous:
        return {}
    log.info(f"Resolving EC segments for {len(ambiguous):,} multi-EC accessions")
    print(f"Resolving EC segments for {len(ambiguous):,} multi-EC accessions...")

    ec_segments: dict = {}
    with ThreadPoolExecutor(max_workers=UNIPROT_WORKERS) as executor:
        futures = {executor.submit(_fetch_ec_segments_one, acc): acc for acc in ambiguous}
        for future in as_completed(futures):
            acc = futures[future]
            try:
                segments = future.result()
            except Exception as e:
                log.warning(f"EC segment fetch failed for {acc}: {e}")
                continue
            if segments:
                ec_segments[acc] = segments

    log.info(f"Resolved EC segments for {len(ec_segments):,}/{len(ambiguous):,} accessions")
    return ec_segments


def main(verbose: bool = False) -> None:
    """Fetch protein names and EC numbers for every UniProt AC linked to a PDB entry.

    Isoform suffixes are stripped before querying, since isoforms share the
    canonical entry's name and EC. Results are saved as both a pickle and a CSV.
    """
    global log
    log = setup_logging(verbose, os.path.join(LOG_DIR, 'mapping', 'fetch_uniprot_names.log'))
    os.makedirs(RESULT_MAPPING_DIR, exist_ok=True)

    pdb2unp = load_pickle(os.path.join(RESULT_MAPPING_DIR, 'pdb2unp.p'))
    # Flatten all UniProt ACs from the PDB→UniProt map.
    raw_accs = {unp for unps in pdb2unp.values() for unp in unps}
    # Strip isoform suffixes (e.g. 'P00533-2' → 'P00533') and deduplicate.
    unique_accs = sorted({_ISOFORM_RE.sub("", a) for a in raw_accs})
    total = len(unique_accs)
    log.info(f"Fetching names for {total:,} canonical UniProt accessions")
    print(f"Fetching UniProt names for {total:,} accessions...")

    # Split into fixed-size chunks for the batch API call.
    chunks = [unique_accs[i:i + UNIPROT_CHUNK_SIZE] for i in range(0, total, UNIPROT_CHUNK_SIZE)]
    all_rows: list[dict] = []
    errors: list[str] = []
    completed = 0

    with ThreadPoolExecutor(max_workers=UNIPROT_WORKERS) as executor:
        # Each chunk gets its own Session so threads don't share connection state.
        futures = {executor.submit(_fetch_chunk, chunk, create_session()): i
                   for i, chunk in enumerate(chunks)}
        for future in as_completed(futures):
            completed += 1
            try:
                all_rows.extend(future.result())
            except Exception as e:
                errors.append(str(e))
                log.warning(f"Chunk failed: {e}")
            if completed % 10 == 0 or completed == len(chunks):
                print(f"  [{completed}/{len(chunks)}] chunks done  ({len(all_rows):,} names so far)")

    if errors:
        log.warning(f"{len(errors)} chunk(s) failed")

    if not all_rows:
        log.error("No rows collected — all chunks failed")
        return

    # Drop duplicates that may arise when isoform and canonical were both requested.
    df = pd.DataFrame(all_rows).drop_duplicates(subset="unp_acc").reset_index(drop=True)
    resolved = (df["protein_name"].notna() & (df["protein_name"] != "")).sum()
    log.info(f"Resolved: {resolved:,}/{total:,} ({100 * resolved / total:.1f}%)")

    csv_path = os.path.join(RESULT_MAPPING_DIR, "uniprot_names.csv")
    pkl_path = os.path.join(RESULT_MAPPING_DIR, "uniprot_names.p")
    df.to_csv(csv_path, index=False)
    save_pickle(df, pkl_path)
    log.info(f"Saved {len(df):,} rows → {pkl_path}")
    print(f"Saved {len(df):,} UniProt names → {RESULT_MAPPING_DIR}")

    # Second pass: disambiguate pooled multi-EC accessions into per-cleavage-product
    # residue segments (see module docstring).
    ec_segments = _build_ec_segments(df)
    ec_seg_pkl = os.path.join(RESULT_MAPPING_DIR, "uniprot_ec_segments.p")
    ec_seg_csv = os.path.join(RESULT_MAPPING_DIR, "uniprot_ec_segments.csv")
    save_pickle(ec_segments, ec_seg_pkl)
    seg_rows = [
        {"unp_acc": acc, "unp_start": start, "unp_end": end, "ec_numbers": ";".join(ecs)}
        for acc, segments in ec_segments.items()
        for start, end, ecs in segments
    ]
    pd.DataFrame(seg_rows).to_csv(ec_seg_csv, index=False)
    log.info(f"Saved EC segments for {len(ec_segments):,} accessions → {ec_seg_pkl}")
    print(f"Saved EC segments for {len(ec_segments):,} accessions → {RESULT_MAPPING_DIR}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(args.verbose)
