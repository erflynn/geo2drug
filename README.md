
### geo2drug
### E Flynn
### Last updated: 3/5/2019

The code in this repository is used for mapping GEO datasets to drugs (using DrugBank IDs or other).

#### Database extraction:
This code is used to extract relevant data from DrugBank, MeSH, and Cellosaurus
* `drugbank_synonyms.py`: Extract drugbank data and synonyms
* `mesh_to_ids.py`: Downloads information from MeSH data
* `parse_pubmed_xml.py`: Parse pubmed XML
* `parse_cellosaurus.py`: Parse cell line data from cellosaurus

#### Data mapping:
* `geo2pubmed_id.R`: Takes geo data, maps to pubmed
* `extract_sex_labels.R`: Extract sex label data from drug mesh (run server-side after GSE download)

#### Data analysis:
* `drug_gse_annot.R`: Takes the results of the annotations and visualizes the results.
* `examine_drug_freq.R`: Looks at the frequency of drugs, visualizes the results

#### Data files:
 * `data/db_data/`: contains output of database downloads and parsing
 * `data/ref_data/`: contains reference data from other sources
 