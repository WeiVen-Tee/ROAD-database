"""
Bundle the assembly CIF files behind the non-redundant annotated site table.

Writes a tar archive: result/nonredundant/PDB.tar.gz, holding one plain `{PDB_ID}-assembly1.cif`
per unique `source_pdb` in the enriched annotated table, plus
PDB_allosteric.tar.gz for the subset of structures with at least one allosteric site.

Usage:
    python copy_nonredundant_pdb.py [--cif-dir DIR] [-v]

Inputs:  result/nonredundant/nonredundant_sites_annotated_enriched.csv, CIF files
Outputs: result/nonredundant/PDB.tar.gz             (members under PDB/)
         result/nonredundant/PDB_manifest.txt       (sorted member names)
         result/nonredundant/PDB_missing.txt        (IDs with no CIF in --cif-dir)
         result/nonredundant/PDB_allosteric.tar.gz  (members under PDB_allosteric/)
         result/nonredundant/PDB_allosteric_manifest.txt
         result/nonredundant/PDB_allosteric_missing.txt
"""

import argparse
import csv
import os
import sys
import tarfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import RESULT_NONREDUNDANT_DIR, CIF_DIR, LOG_DIR
from utils import setup_logging

log = None  # set in main()

GZIP_LEVEL = 6


def _read_source_pdbs(csv_path: str) -> tuple:
    """Return (all_ids, allosteric_ids), sorted uppercase source_pdb values.

    A structure counts as allosteric if any of its rows does, so the second list is
    a subset of the first.
    """
    all_ids, allosteric_ids = set(), set()
    with open(csv_path, newline='') as fh:
        for row in csv.DictReader(fh):
            pdb_id = (row.get('source_pdb') or '').strip().upper()
            if not pdb_id:
                continue
            all_ids.add(pdb_id)
            if (row.get('site_category') or '').strip() == 'allosteric':
                allosteric_ids.add(pdb_id)
    return sorted(all_ids), sorted(allosteric_ids)


def _normalise(ti: tarfile.TarInfo) -> tarfile.TarInfo:
    """Strip the building machine's ownership and clock precision from a member.

    Whole-second mtimes matter beyond tidiness: a fractional one does not fit a plain
    ustar record, so tarfile emits a ~512-byte pax extended header for every member
    that has one.
    """
    ti.uid = ti.gid = 0
    ti.uname = ti.gname = ''
    ti.mode = 0o644
    ti.mtime = int(ti.mtime)
    return ti


def _write_bundle(pdb_ids: list, dest_tgz: str, arc_prefix: str, manifest_path: str,
                  missing_path: str, cif_dir: str) -> str:
    """Archive one plain CIF per ID into dest_tgz; return a summary line.

    Members are stored uncompressed under arc_prefix/ and the whole tar is gzipped. 
    IDs with no CIF in cif_dir are written to missing_path rather than raising.
    """
    written, missing, raw_bytes = [], [], 0
    print(f'{os.path.basename(dest_tgz)}: archiving {len(pdb_ids):,} structures', flush=True)
    try:
        with tarfile.open(dest_tgz, 'w:gz', compresslevel=GZIP_LEVEL) as tar:
            for i, pdb_id in enumerate(pdb_ids, 1):
                name = f'{pdb_id}-assembly1.cif'
                src = os.path.join(cif_dir, name)
                if not os.path.exists(src):
                    missing.append(pdb_id)
                    continue
                tar.add(src, arcname=f'{arc_prefix}/{name}', filter=_normalise)
                raw_bytes += os.path.getsize(src)
                written.append(name)

                if i % 2000 == 0:
                    print(f'  {i:,}/{len(pdb_ids):,}', flush=True)
    except BaseException:
        # Leave no half-written archive behind: the pipeline's cache-skip check
        # tests only for the file's existence and would accept a truncated one.
        if os.path.exists(dest_tgz):
            os.remove(dest_tgz)
        raise

    with open(manifest_path, 'w') as fh:
        for name in written:
            fh.write(f'{name}\n')

    # Written unconditionally (empty when nothing is missing) so a stale list from a
    # previous run can never be mistaken for the current one.
    with open(missing_path, 'w') as fh:
        for pdb_id in missing:
            fh.write(f'{pdb_id}\n')

    size = os.path.getsize(dest_tgz)
    summary = (f'{os.path.basename(dest_tgz)}: requested: {len(pdb_ids):,}  '
               f'archived: {len(written):,}  missing: {len(missing):,}  '
               f'({size / 1e9:.2f} GB from {raw_bytes / 1e9:.2f} GB of CIFs)')
    log.info(summary)
    if missing:
        log.warning(f'{len(missing):,} structures had no CIF in {cif_dir!r} '
                    f'— see {missing_path}')
    return summary


def main(cif_dir: str = CIF_DIR, verbose: bool = False) -> None:
    global log
    os.makedirs(LOG_DIR, exist_ok=True)
    log = setup_logging(verbose, os.path.join(LOG_DIR, 'copy_nonredundant_pdb.log'))

    in_csv = os.path.join(RESULT_NONREDUNDANT_DIR,
                          'nonredundant_sites_annotated_enriched.csv')
    if not os.path.exists(in_csv):
        log.error(f'Input table not found: {in_csv}')
        sys.exit(1)

    pdb_ids, allosteric_ids = _read_source_pdbs(in_csv)
    log.info(f'{len(pdb_ids):,} unique source_pdb entries in {os.path.basename(in_csv)} '
             f'({len(allosteric_ids):,} with at least one allosteric site)')

    summaries = [_write_bundle(
        pdb_ids, os.path.join(RESULT_NONREDUNDANT_DIR, 'PDB.tar.gz'), 'PDB',
        os.path.join(RESULT_NONREDUNDANT_DIR, 'PDB_manifest.txt'),
        os.path.join(RESULT_NONREDUNDANT_DIR, 'PDB_missing.txt'),
        cif_dir,
    ), _write_bundle(
        allosteric_ids, os.path.join(RESULT_NONREDUNDANT_DIR, 'PDB_allosteric.tar.gz'),
        'PDB_allosteric',
        os.path.join(RESULT_NONREDUNDANT_DIR, 'PDB_allosteric_manifest.txt'),
        os.path.join(RESULT_NONREDUNDANT_DIR, 'PDB_allosteric_missing.txt'),
        cif_dir,
    )]
    print('\n'.join(summaries))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cif-dir", default=CIF_DIR)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()
    main(args.cif_dir, args.verbose)
