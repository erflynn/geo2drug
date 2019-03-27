
### geo2drug
### E Flynn
### Last updated: 3/5/2019

The code in this repository is used for mapping GEO datasets to drugs (using DrugBank IDs or other).

#### Database extraction:
This code is used to extract relevant data from DrugBank, MeSH, and Cellosaurus
* Pubmed XML: `parse_pubmed_xml.py`
* Drugbank: `drugbank_synonyms.py`, `processDrugbank.R`: Extract drugbank data and synonyms
* MeSH: `mesh_to_ids.py`, `condenseMeSH.R`: Downloads information from MeSH data
* Cellosaurus: `parse_cellosaurus.py`: Parse cell line data from cellosaurus

#### Data mapping:
* `geo2pubmed_id.R`: Takes geo data, maps to mesh
* `mesh2drugbank.R`: Maps MeSH IDs to drugbank
* `extract_sex_labels.R`: Extract sex label data from drug mesh (run server-side after GSE download)

#### Data analysis:
* `drug_gse_annot.Rmd`: Takes the results of the annotations and visualizes the results.
* `examine_drug_freq.Rmd`: Looks at the frequency of drugs, visualizes the results

#### Data files:
 * `data/db_data/`: contains output of database downloads and parsing
 * `data/ref_data/`: contains reference data from other sources
 