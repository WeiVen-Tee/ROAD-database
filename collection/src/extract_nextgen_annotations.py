"""
Extract NextGen-exclusive SIFTS annotations from PDB NextGen Archive.

Downloads xyz-enrich.cif.gz per PDB entry and extracts 2 categories as flat CSVs:
  pdbx_sifts_unp_segments.csv        — UniProt segment mappings
  pdbx_sifts_xref_db_segments.csv    — External DB segment mappings

Usage:
    python extract_nextgen_annotations.py [--cif-dir DIR] [--workers N] [--no-del] [-v]

The CIF directory defines the set of PDB entries to process; it defaults to
CIF_DIR from config.py and can be overridden with --cif-dir.
"""

import argparse
import csv
import functools
import gzip
import os
import re
import shlex
import socket
import sys
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import (
    CIF_DIR, NEXTGEN_OUT_DIR, NEXTGEN_BASE_URL,
    NEXTGEN_TARGET_CATEGORIES, NEXTGEN_CATEGORY_CSV, LOG_DIR,
)
from utils import pdb_ids_from_cif_dir, setup_logging

_NEXTGEN_DOWNLOAD_DIR = os.path.join(os.path.dirname(NEXTGEN_OUT_DIR), 'mmcif_nextgen')

_locks: dict[str, threading.Lock] = {}
_lock_guard = threading.Lock()
log = None  # set in main()


def _get_lock(path: str) -> threading.Lock:
    """Per-file lock, created on first use, so writes to different CSVs don't serialise."""
    with _lock_guard:
        if path not in _locks:
            _locks[path] = threading.Lock()
    return _locks[path]


def _nextgen_url(pdb_id: str) -> str:
    """'1ABC' -> '.../ab/pdb_00001abc/pdb_00001abc_xyz-enrich.cif.gz'."""
    low = pdb_id.lower()
    accession = f"pdb_0000{low}"
    return f"{NEXTGEN_BASE_URL}/{low[1:3]}/{accession}/{accession}_xyz-enrich.cif.gz"


def _download_enrich(pdb_id: str, retries: int = 3, backoff: float = 2.0) -> str | None:
    """Download one xyz-enrich.cif.gz with exponential backoff.

    Returns the local path, an already-downloaded file if present, or None if
    every attempt fails.
    """
    low = pdb_id.lower()
    dest = os.path.join(_NEXTGEN_DOWNLOAD_DIR, f"pdb_0000{low}_xyz-enrich.cif.gz")
    if os.path.exists(dest):
        return dest
    for attempt in range(1, retries + 1):
        try:
            urllib.request.urlretrieve(_nextgen_url(pdb_id), dest)
            return dest
        except Exception as e:
            try:
                os.remove(dest)
            except FileNotFoundError:
                pass
            if attempt < retries:
                wait = backoff * (2 ** (attempt - 1))  # Exponential backoff.
                log.warning(f"{pdb_id}: download failed (attempt {attempt}/{retries}), retrying in {wait:.0f}s — {e}")
                time.sleep(wait)
            else:
                log.warning(f"{pdb_id}: download failed after {retries} attempts — {e}")
    return None


def _parse_row(line: str) -> list[str]:
    """Tokenise a CIF data row, honouring quotes; falls back to split() on bad quoting."""
    try:
        return shlex.split(line)
    except ValueError:
        return line.split()


def _parse_enrich(filepath: str) -> dict:
    """Extract the NEXTGEN_TARGET_CATEGORIES rows from an xyz-enrich.cif(.gz) file.

    Handles both mmCIF record forms: inline key-value tags and loop_ tables.
    Returns {category: (fields, rows)}.
    """
    results = {cat: ([], []) for cat in NEXTGEN_TARGET_CATEGORIES}
    opener = gzip.open if filepath.endswith(".gz") else open

    with opener(filepath, "rt", encoding="utf-8", errors="replace") as fh:
        in_loop = False
        active_cats: set[str] = set()
        loop_fields: dict[str, list[str]] = {}

        for raw in fh:
            stripped = raw.strip()

            if stripped == "loop_":
                in_loop = True
                active_cats = set()
                loop_fields = {}
                continue

            if not in_loop:
                # Handle non-loop inline key-value items (single-entry records).
                token = stripped.split()[0] if stripped else ""
                # Match longest category prefix first to avoid ambiguous matches.
                for cat in sorted(NEXTGEN_TARGET_CATEGORIES, key=len, reverse=True):
                    if token.startswith(cat + "."):
                        field = token[len(cat) + 1:]
                        value = stripped[len(token):].strip()
                        fields, rows = results[cat]
                        if not fields:
                            fields.append(field)
                            rows.append([value])
                        else:
                            fields.append(field)
                            rows[-1].append(value)
                        break
                continue

            if stripped.startswith("_"):
                # Column-name line inside a loop block.
                token = stripped.split()[0]
                for cat in sorted(NEXTGEN_TARGET_CATEGORIES, key=len, reverse=True):
                    if token.startswith(cat + "."):
                        active_cats.add(cat)
                        loop_fields.setdefault(cat, []).append(token[len(cat) + 1:])
                        break
                continue

            if not stripped or stripped in ("#",) or stripped.startswith(("data_", "save_")):
                # Any of these signals the end of the current loop block.
                in_loop = False
                active_cats = set()
                continue

            if active_cats:
                values = _parse_row(stripped)
                for cat in active_cats:
                    fields, rows = results[cat]
                    if not fields:
                        # First data row — copy the collected column names into fields.
                        fields.extend(loop_fields.get(cat, []))
                    rows.append(values)
            else:
                in_loop = False

    return results


def _append_csv(out_path: str, pdb_id: str, fields: list[str], rows: list[list[str]]) -> None:
    """Append rows to a shared CSV under its per-file lock, writing the header if new.

    pdb_id is prepended as the first column of every row.
    """
    lock = _get_lock(out_path)
    with lock:
        is_new = not os.path.exists(out_path) or os.path.getsize(out_path) == 0
        with open(out_path, "a", newline="") as fh:
            writer = csv.writer(fh)
            if is_new:
                writer.writerow(["pdb_id"] + fields)
            for row in rows:
                writer.writerow([pdb_id] + row)


def _load_skip_set() -> frozenset[str]:
    """PDB IDs already present in the output CSVs, for --resume."""
    skip: set[str] = set()
    for cat_csv in NEXTGEN_CATEGORY_CSV.values():
        path = os.path.join(NEXTGEN_OUT_DIR, cat_csv)
        try:
            with open(path, newline="") as fh:
                reader = csv.reader(fh)
                next(reader)  # skip header
                for row in reader:
                    if row:
                        skip.add(row[0].upper())
        except FileNotFoundError:
            pass
    return frozenset(skip)


def _process_one(pdb_id: str, delete_after: bool, skip_set: frozenset[str] = frozenset()) -> tuple[str, int]:
    """Download, parse and append one PDB entry's NextGen SIFTS data.

    Returns (pdb_id, rows_written); 0 on download failure, -1 when skipped via
    skip_set. delete_after removes the .cif.gz once parsed.
    """
    if pdb_id in skip_set:
        return pdb_id, -1

    filepath = _download_enrich(pdb_id)
    if filepath is None:
        return pdb_id, 0

    parsed = _parse_enrich(filepath)
    total_rows = 0

    for cat, (fields, rows) in parsed.items():
        if not fields or not rows:
            continue
        _append_csv(os.path.join(NEXTGEN_OUT_DIR, NEXTGEN_CATEGORY_CSV[cat]), pdb_id, fields, rows)
        total_rows += len(rows)

    if delete_after:
        try:
            os.remove(filepath)
        except OSError:
            pass

    return pdb_id, total_rows


def _parse_failed_from_log(log_path: str) -> list[str]:
    """PDB IDs that failed every download attempt in a prior run, read back from the log."""
    pattern = re.compile(r"(\w{4}): download failed after \d+ attempts")
    failed: set[str] = set()
    try:
        with open(log_path) as fh:
            for line in fh:
                m = pattern.search(line)
                if m:
                    failed.add(m.group(1).upper())
    except FileNotFoundError:
        pass
    return list(failed)


def main(cif_dir: str = CIF_DIR, workers: int = 32, delete: bool = True,
         verbose: bool = False, retry_failed: bool = False, resume: bool = False,
         pdb_ids_arg: list[str] | None = None, force: bool = False) -> None:
    """Download and extract NextGen SIFTS annotations for every PDB entry in cif_dir.

    Output CSVs are appended to as each entry completes, so an interrupted run
    keeps its partial results.

    Args:
        cif_dir:      Directory of {PDB_ID}-assembly1.cif files defining the entry set.
        retry_failed: Re-process only the IDs that failed in the previous run.
        resume:       Skip IDs already present in the output CSVs.
        pdb_ids_arg:  Process only these IDs instead of all in cif_dir.
        force:        Truncate the output CSVs first, for a clean full rebuild —
                      required because the writes are appends. Ignored for the
                      three targeted modes above, which add to what exists.
    """
    global log
    log_path = os.path.join(LOG_DIR, 'mapping', 'extract_nextgen_annotations.log')
    log = setup_logging(verbose, log_path)

    socket.setdefaulttimeout(30)

    os.makedirs(_NEXTGEN_DOWNLOAD_DIR, exist_ok=True)
    os.makedirs(NEXTGEN_OUT_DIR, exist_ok=True)

    # Full rebuild: clear the append-mode outputs first (see the `force` arg docs).
    # The targeted modes below deliberately add to existing CSVs, so skip them.
    if force and not (resume or retry_failed or pdb_ids_arg):
        for cat_csv in NEXTGEN_CATEGORY_CSV.values():
            path = os.path.join(NEXTGEN_OUT_DIR, cat_csv)
            if os.path.exists(path):
                os.remove(path)
                log.info(f"Force mode — removed existing {cat_csv}")

    if pdb_ids_arg:
        pdb_ids = [p.upper() for p in pdb_ids_arg]
        log.info(f"Processing {len(pdb_ids)} specified PDB IDs...")
        print(f"Processing {len(pdb_ids)} specified PDB IDs ({workers} workers)...")
    elif retry_failed:
        pdb_ids = _parse_failed_from_log(log_path)
        if not pdb_ids:
            log.info("No failed downloads found in log.")
            print("No failed downloads found in log.")
            return
        log.info(f"Retrying {len(pdb_ids)} failed PDB IDs...")
        print(f"Retrying {len(pdb_ids)} failed downloads ({workers} workers)...")
    else:
        if not os.path.isdir(cif_dir):
            log.warning(f"CIF directory not found: {cif_dir!r} — set STRUCTURE_CIF_DIR, "
                        f"create src/config_local.py, or pass --cif-dir")
            return
        pdb_ids = pdb_ids_from_cif_dir(cif_dir)
        if not pdb_ids:
            log.warning(f"No *-assembly1.cif files found in: {cif_dir}")
            return
        log.info(f"Found {len(pdb_ids)} PDB IDs in {cif_dir}")
        print(f"Extracting NextGen annotations for {len(pdb_ids)} PDB IDs ({workers} workers)...")

    skip_set: frozenset[str] = frozenset()
    if resume:
        skip_set = _load_skip_set()
        log.info(f"Resume mode: {len(skip_set)} PDB IDs already in CSVs will be skipped")
        print(f"  Resume: skipping {len(skip_set)} already-processed IDs")

    total = len(pdb_ids)
    completed = 0
    worker_fn = functools.partial(_process_one, delete_after=delete, skip_set=skip_set)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(worker_fn, pid): pid for pid in pdb_ids}
        for future in as_completed(futures):
            pdb_id, n = future.result()
            if n >= 0:
                log.debug(f"{pdb_id}: {n} rows written")
            completed += 1
            if completed % 1000 == 0 or completed == total:
                print(f"  [{completed}/{total}] done")

    log.info(f"Done. CSV files written to: {NEXTGEN_OUT_DIR}")
    print(f"Done. CSV files → {NEXTGEN_OUT_DIR}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cif-dir", default=CIF_DIR)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--no-del", dest="delete", action="store_false",
                        help="Keep downloaded xyz-enrich.cif.gz files after parsing")
    parser.set_defaults(delete=True)
    parser.add_argument("--retry-failed", action="store_true",
                        help="Retry only PDB IDs that failed in a previous run (parsed from log)")
    parser.add_argument("--resume", action="store_true",
                        help="Skip PDB IDs already present in output CSVs; only processes missing entries")
    parser.add_argument("--pdb-ids", nargs="+", metavar="PDB_ID",
                        help="Process only these specific PDB IDs instead of all IDs in the CIF directory")
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-f", "--force", action="store_true",
                        help="Truncate output CSVs first for a clean full rebuild")
    args = parser.parse_args()
    main(cif_dir=args.cif_dir, workers=args.workers, delete=args.delete, verbose=args.verbose,
         retry_failed=args.retry_failed, resume=args.resume, pdb_ids_arg=args.pdb_ids,
         force=args.force)
