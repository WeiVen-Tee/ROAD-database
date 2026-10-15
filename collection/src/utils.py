"""
Helpers shared across every step of the collection pipeline.

Grouped by concern: pickle I/O (load_pickle/save_pickle), logging setup
(setup_logging), HTTP sessions with retry (create_session), CIF parsing
(cif2seq, pdb_ids_from_cif_dir), UniProt accessions (canonical_accession),
and the UnionFind used by the non-redundant clustering step.

`../analysis` imports setup_logging, create_session, load_pickle, save_pickle,
cif2seq and UnionFind from here directly, so this module is part of the
cross-pipeline contract — see the re-export note in analysis/src/config.py
before renaming anything it imports from config.
"""

import logging
import os
import pickle
import re
from collections import OrderedDict
from functools import lru_cache
from typing import NamedTuple

import gemmi
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


class ChainSeq(NamedTuple):
    """Per-chain sequence/numbering data returned by cif2seq(), one entry per residue."""
    seq: str                      # one-letter sequence, always upper-case; unknown residues are 'X'
    auth_nums: list                # residue.seqid.num — PDB residue number used to pick atoms
    label_seqs: list               # residue.label_seq — entity-sequence position (incl. unmodeled residues)
    entity_ids: list               # entity.name resolved via residue.subchain; None if unresolved
    icodes: list                   # residue.seqid.icode, stripped; '' when the residue has none


class UnionFind:
    """Disjoint-set structure with path splitting, keyed by any hashable."""

    def __init__(self, keys):
        self.parent = {k: k for k in keys}

    def find(self, x):
        """Return the representative of the set containing x."""
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        """Merge the sets containing x and y."""
        rx, ry = self.find(x), self.find(y)
        if rx != ry:
            self.parent[ry] = rx


def save_pickle(obj, filepath: str) -> None:
    """Save object to pickle file."""
    try:
        with open(filepath, 'wb') as f:
            pickle.dump(obj, f)
    except Exception as e:
        raise IOError(f"Failed to save pickle to {filepath}: {e}")


def load_pickle(filepath: str):
    """Load object from pickle file."""
    try:
        with open(filepath, 'rb') as f:
            return pickle.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Pickle file not found: {filepath}")
    except Exception as e:
        raise IOError(f"Failed to load pickle from {filepath}: {e}")


def create_session() -> requests.Session:
    """Create a requests Session with connection pooling and retry strategy."""
    session = requests.Session()
    adapter = HTTPAdapter(
        pool_connections=20,
        pool_maxsize=20,
        max_retries=Retry(
            total=3,
            backoff_factor=0.3,
            status_forcelist=[429, 500, 502, 503, 504],
        ),
    )
    session.mount('http://',  adapter)
    session.mount('https://', adapter)
    return session


_ISOFORM_RE = re.compile(r'-\d+$')


def canonical_accession(acc: str) -> str:
    """UniProt accession without its isoform suffix ('P00533-2' -> 'P00533')."""
    return _ISOFORM_RE.sub('', acc)


_CIF_ID_RE = re.compile(r'^([A-Za-z0-9]{4})-assembly1\.cif$', re.IGNORECASE)


def pdb_ids_from_cif_dir(cif_dir: str) -> list[str]:
    """Return sorted unique PDB IDs extracted from *-assembly1.cif filenames."""
    seen: set[str] = set()
    ids: list[str] = []
    for fname in sorted(os.listdir(cif_dir)):
        m = _CIF_ID_RE.match(fname)
        if m:
            pid = m.group(1).upper()
            if pid not in seen:
                seen.add(pid)
                ids.append(pid)
    return ids


def setup_logging(verbose: bool, log_path: str) -> logging.Logger:
    """Configure root logger.

    File handler: INFO and above (DEBUG suppressed).
    Console handler: INFO when verbose, WARNING when quiet.
    Progress updates should use print() directly — always visible.
    """
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()

    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO if verbose else logging.WARNING)
    ch.setFormatter(logging.Formatter('%(levelname)s %(message)s'))

    logger.addHandler(fh)
    logger.addHandler(ch)
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    return logger


@lru_cache(maxsize=1000)
def cif2seq(cif_id: str, filepath: str):
    """Per-chain sequence and residue numbering from an mmCIF file.

    Returns an OrderedDict mapping 'MODEL.CHAIN' → ChainSeq. Only residues
    covalently part of a polymer chain (label_seq is not None) are included, so
    modified residues such as SEP/TPO/MSE are kept but free-standing amino acids
    bound as ligands are not. Alternative conformations are collapsed to the
    first conformer, and chains shorter than MIN_CHAIN_LENGTH are dropped.

    cif_id is part of the lru_cache key; filepath is the file actually read.
    """
    from config import MIN_CHAIN_LENGTH
    structure = gemmi.read_structure(filepath)
    structure.remove_alternative_conformations()
    # label_asym_id (residue.subchain) → entity_id: symmetry copies of one
    # molecule get separate subchains but share an entity.
    subchain_to_entity = {sc: entity.name for entity in structure.entities for sc in entity.subchains}
    chains = OrderedDict()
    for model_idx, model in enumerate(structure):
        model_id = str(model_idx + 1)
        for chain in model:
            seq = []
            auth_nums = []
            label_seqs = []
            entity_ids = []
            icodes = []
            for residue in chain:
                if residue.label_seq is None:   # non-polymer: ligand, water, ion
                    continue
                res_info = gemmi.find_tabulated_residue(residue.name)
                atoms = [a for a in residue if a.element.name != 'H']
                if not atoms:
                    continue
                one_letter = res_info.one_letter_code if res_info else '\x00'
                if not one_letter.strip() or one_letter == '\x00':
                    fc = res_info.fasta_code() if res_info else 'X'
                    one_letter = fc if fc and fc.strip() else 'X'
                # gemmi lower-cases modified residues; callers compare sequences
                # with case-sensitive edit distance.
                seq.append(one_letter.upper())
                auth_nums.append(residue.seqid.num)
                label_seqs.append(residue.label_seq)
                entity_ids.append(subchain_to_entity.get(residue.subchain))
                icodes.append(residue.seqid.icode.strip())
            if len(seq) >= MIN_CHAIN_LENGTH:
                chains[f'{model_id}.{chain.name}'] = ChainSeq(
                    ''.join(seq), auth_nums, label_seqs, entity_ids, icodes)
    return chains
