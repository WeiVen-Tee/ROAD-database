# collection

Data collection pipeline: builds the non-redundant, annotated ligand-binding-site
table that `../analysis` consumes. Maps PDB structures to UniProt via SIFTS,
extracts ligand-binding sites from structure files, clusters near-duplicate
sites across redundant PDB entries, and enriches the result with distance,
structural, and functional-annotation metadata.

## Prerequisites

Before running the pipeline, these inputs must be in place:

### 1. PDB assembly CIF files

One `{PDB_ID}-assembly1.cif` per entry. See **Setup** below for how to point the
pipeline at your local copy.

### 2. PDB Chemical Component Dictionary (CCD)

1. Go to <https://www.wwpdb.org/data/ccd>
2. Download **Chemical Component Dictionary: mmCIF (gz)**
3. Decompress and rename the file to `Components-rel-alt.cif`, then place it at
   `collection/data/Components-rel-alt.cif`

### 3. BioLiP artifact list

Already included at `collection/data/possible_artifact_molecule_list.txt` — no action
needed.

### 4. Downloaded PDB entry list

Already included at `collection/data/downloaded_PDB_{1,2}.txt` — no action needed.
The two files together record the 128,141 PDB entries that were present in the CIF
directory for the published run, split only for transfer convenience. They document
the input inventory for reproducibility; no pipeline step reads them.

---

## Setup

`src/config.py` needs a directory of PDB assembly CIF files
(`{PDB_ID}-assembly1.cif`) to run the CIF-dependent steps (`nextgen_annotations`,
`uniprot_mapping`, `site_extraction`, `functional_site_mapping`, `distances`,
`nonredundant_build_annotated`, `nonredundant_build_complete`, `protein_names`,
`nonredundant_pdb_copy`).
This path is machine-specific and is **not** hardcoded — set it via:

```bash
export STRUCTURE_CIF_DIR=/path/to/your/PDB/assembly/cifs/
```

or copy `src/config_local.py.example` to a gitignored `src/config_local.py`:

```bash
cp src/config_local.py.example src/config_local.py
# then edit CIF_DIR in src/config_local.py
```

## Running

```bash
python run_collection_pipeline.py                                    # full run, all steps
python run_collection_pipeline.py --steps site_extraction site_filter # run specific steps only
python run_collection_pipeline.py -v                                  # verbose console logging
python run_collection_pipeline.py -f                                  # force re-run, ignore cached outputs
```

A step is skipped when all of its outputs already exist; pass `-f`/`--force` to
recompute. `--force` also clears the caches individual steps keep for themselves
(the distances checkpoint, `pdb_cif_cache.p`, the UniProt feature cache) and
truncates `nextgen_annotations`' output CSVs, which are otherwise appended to.

## Steps

| Step                          | Script                                  | Output                                                              |
|---------------------------------|--------------------------------------------|------------------------------------------------------------------------|
| `nextgen_annotations`           | extract_nextgen_annotations                 | `result/mapping/nextgen_annotations/`                                  |
| `uniprot_mapping`                | map_pdb_nextgen_annotations                  | `result/mapping/mapped_uniprot.csv`, `mapped_cath.csv`, `mapped_scop.csv` |
| `chain_unp_map`                  | build_chain_unp_map                          | `result/mapping/pdb2unp.p`, `pdb_chain_to_unp{,_all}.p`                 |
| `uniprot_names`                   | fetch_uniprot_names                           | `result/mapping/uniprot_names.p`                                         |
| `ligand_blacklist`                | build_badlig                                  | `result/ligands/excluded_lig.p`                                          |
| `site_extraction`                 | find_resi_cif                                 | `result/binding_sites/site_res.txt`                                      |
| `site_filter`                      | filter_binding_sites                          | `result/binding_sites/site_res_db.p`                                      |
| `uniprot_features`                 | fetch_uniprot_features                        | `result/functional/unp2site.p`                                            |
| `functional_site_mapping`          | map_functional_sites                          | `result/functional/func_site_dict.p`                                       |
| `distances`                         | calculate_distances                           | `result/distances/distance_results.p`                                       |
| `nonredundant_build_annotated`      | build_nonredundant_sites (annotated)          | `result/nonredundant/_nonredundant_sites_annotated.{p,csv}`, `pdb_cif_cache.p` |
| `nonredundant_build_complete`       | build_nonredundant_sites (complete)           | `result/nonredundant/_nonredundant_sites_complete.{p,csv}`                    |
| `protein_names`                      | better_protein_names (Protein Identity)       | `result/protein_names/protein_names_pdb.csv`                                  |
| `nonredundant_enrich_annotated`     | enrich_nonredundant_sites (annotated)         | `result/nonredundant/nonredundant_sites_annotated_enriched.{p,csv}`           |
| `nonredundant_enrich_complete`      | enrich_nonredundant_sites (complete)          | `result/nonredundant/nonredundant_sites_complete_enriched.{p,csv}`            |
| `nonredundant_pdb_copy`              | copy_nonredundant_pdb                          | `result/nonredundant/PDB{,_allosteric}.tar.gz`, `PDB{,_allosteric}_manifest.txt` |

`nonredundant_pdb_copy` bundles the assembly CIF behind every `source_pdb` in
`nonredundant_sites_annotated_enriched.csv` (~23.7k structures, ~4 GB) into
`result/nonredundant/PDB.tar.gz`. 
Members are plain, uncompressed `.cif` files under a `PDB/` prefix.

The same step then writes `result/nonredundant/PDB_allosteric.tar.gz` (~5.2k
structures), holding the structures with at least one row of `site_category == 'allosteric'`.

`nonredundant_sites_annotated_enriched.csv` (from `nonredundant_enrich_annotated`)
is the primary input consumed by `../analysis`.
