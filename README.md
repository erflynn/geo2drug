
### geo2drug
### E Flynn
### Last updated: 9/25/2019

The code in this repository is used for mapping GEO datasets to drugs (using DrugBank IDs or other).


#### GEO
GEO metadata is downloaded using `get_list_all_gse.R`
The GEO data itself is then downloaded server-side using the scripts in the `geo_download` directory.

#### Cell line labeling
  * `parse_cellosaurus.py`
  * `process_cell_df.R`
  * `cell_line_labeling.R`  labels cell line in all GSE metadata

#### Drug labeling
  * `drugbank_synonyms.py`
  * `process_drugbank.R`
  * `drug_gse_labeling.R` : search the GSE study-level metadata for drug names





#### Data files:
 * `data/db_data/`: contains output of database downloads and parsing
 * `data/ref_data/`: contains reference data from other sources
 

#### Old:
* MetaSRA:
[* Pubmed XML: `parse_pubmed_xml.py`]
[* MeSH: `mesh_to_ids.py`, `condense_mesh.R`: Downloads information from MeSH data]
All of the extracted data then goes into the `db_data` directory.


* `geo2pubmed_id.R`: Takes geo data, maps to mesh
* `mesh2drugbank.R`: Maps MeSH IDs to drugbank

* `extract_sex_labels.R`: Extract sex label data from drug mesh (run server-side after GSE download)
* `drug_gse_annot.Rmd`: Takes the results of the annotations and visualizes the results.
* `examine_drug_freq.Rmd`: Looks at the frequency of drugs, visualizes the results
