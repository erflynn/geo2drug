
### geo2drug
### E Flynn
### Last updated: 9/25/2019

The code in this repository is used for mapping GEO datasets to drugs (using DrugBank IDs or other).


#### GEO Metadata and download
GEO metadata is downloaded using `get_list_all_gse.R`
The GEO data itself is then downloaded server-side using the scripts in the `geo_download` directory.

#### Cell line labeling
  * `parse_cellosaurus.py`
  * `process_cell_df.R`
  * `cell_line_labeling.R` 
  
#### Drug labeling
  * `drugbank_synonyms.py`
  * `process_drugbank.R`
  * `drug_gse_labeling.R` 

#### Data files:
 * `data/db_data/`: contains output of database downloads and parsing
 * `data/ref_data/`: contains reference data from other sources
 


