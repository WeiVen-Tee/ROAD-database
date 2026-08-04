"""Shared CIF-driven helpers for the per-structure analysis stages."""

import logging
import os
import socket
import urllib.request

# Seconds; applied before each RCSB fetch in resolve_cif_path (see note there).
_DOWNLOAD_TIMEOUT = 30


def iter_site_residues(s: str):
    """Tokenize a binding_site string into (chain, seqnum, icode) triples.

    Input format is 'A.7;A.11;B.100A;...' — each token is 'CHAIN.RESNUM', where
    RESNUM may carry a trailing insertion code ('100A' -> 100, 'A') and may be
    negative. Tokens whose numeric part won't parse are skipped.

    icode is returned separately because callers differ: most key residues on
    (seqnum, icode), but analyze_interactions matches arpeggio's auth_seq_id and
    ignores it.
    """
    for token in str(s).split(';'):
        token = token.strip()
        if not token or '.' not in token:
            continue
        chain, resnum = token.split('.', 1)
        resnum = resnum.strip()
        icode = ''
        if resnum and resnum[-1].isalpha():
            resnum, icode = resnum[:-1], resnum[-1]
        try:
            seqnum = int(resnum)
        except ValueError:
            continue
        yield chain.strip(), seqnum, icode


def parse_site_res_list(s: str) -> list[tuple[str, str]]:
    """Residue list as (chain, resnum_str) pairs, insertion code kept inline."""
    return [(chain, f'{seqnum}{icode}') for chain, seqnum, icode in iter_site_residues(s)]


def _build_pdb_sites(enrich_df, all_pdbs: list[str]) -> dict[str, list[dict]]:
    """Build a per-structure index of binding sites from the enriched DataFrame.

    Each site dict has keys: ID, ligand_code, site_category, residues
    (list of (chain, seqnum, icode) triples from iter_site_residues).
    """
    pdb_set = set(all_pdbs)
    pdb_sites: dict[str, list[dict]] = {p: [] for p in all_pdbs}
    for _, row in enrich_df.iterrows():
        pdb_id = str(row['source_pdb']).upper()
        if pdb_id not in pdb_set:
            continue
        residues = list(iter_site_residues(row['binding_site']))
        pdb_sites[pdb_id].append({
            'ID':            row['ID'],
            'ligand_code':   row['ligand_code'],
            'site_category': row['site_category'],
            'residues':      residues,
        })
    return pdb_sites


def resolve_cif_path(pdb_id: str, cif_dir: str, download_dir: str | None,
                      logger: logging.Logger, rcsb_url_template: str) -> str | None:
    """Locate a structure's assembly CIF, fetching it from RCSB if necessary.

    Looks in cif_dir first. If the file isn't there and download_dir is given,
    the CIF is fetched into download_dir and cached there for later runs.

    Returns None (after logging a warning) if the CIF cannot be found or fetched.
    """
    cif_file = os.path.join(cif_dir, f'{pdb_id.upper()}-assembly1.cif')
    if not os.path.exists(cif_file):
        if download_dir is not None:
            cif_file = os.path.join(download_dir, f'{pdb_id}-assembly1.cif')
            if not os.path.exists(cif_file):
                try:
                    url = rcsb_url_template.format(pdb_id=pdb_id.upper())
                    logger.info('Downloading %s → %s', url, cif_file)
                    # urlretrieve honours no timeout of its own, so a stalled
                    # connection would hang the stage indefinitely.
                    socket.setdefaulttimeout(_DOWNLOAD_TIMEOUT)
                    urllib.request.urlretrieve(url, cif_file)
                except Exception as exc:
                    logger.warning('%s: download failed: %s', pdb_id, exc)
                    return None
        else:
            logger.warning('%s: CIF not found', pdb_id)
            return None
    return cif_file
