### geo2drug
### E Flynn


This repository annotates GEO datasets using DrugBank and Cellosaurus labels.
Then, we download and sex label all data.


The code directory is set up as follows:

`00_annot/`
  1. GEO metadata is downloaded using `00_get_list_all_gse.R`
   - output: (in `data/01_sample_lists/`)
    * `<>_gse_gsm.csv`
    * `gse_<>.csv`
    * `gse_gsm_all_geo_dedup.csv`
    * `gse_all_geo_info.csv`
  2. `01_cell_labeling/`
   * `00_parse_cellosaurus.py`
     - input: cellosaurus XML, located in `data/00_db_data/cellosaurus.xml`
     - output: `data/00_db_data/cellosaurus.json`
   * `01_process_cell_df.R`
     - output: `data/00_db_data/cellosaurus_df.txt`
   * `02_cell_line_labeling.R` 
     - uses `data/01_sample_lists/gse_all_geo_info.csv`
     - output: `data/02_labeled_data/cell_line_mapped_gse.txt`
  3. `02_drug_labeling/`
   * `00_drugbank_synonyms.py`
     - input: `drugbank.xml`
     - output: `data/00_db_data/drugbank_info.json`
   * `01_process_drugbank.R`
     - output: `data/00_db_data/drugbank_parsed.txt`
   * `02_drug_gse_labeling.R` 
     - uses `data/01_sample_lists/gse_all_geo_info.csv`
     - output: `data/02_labeled_data/drugbank_mapped_gse.txt` 

`01_download/`
  This uses the NCBI aspera client and relies on wrenlab software packages and is slightly complicated to install.
  `sbatch 00_download_gse_wrapper.sh ${ID} ${GSE_LIST}`

This runs `download_geo_chunk.sh` which runs `downloadGEO.py` for each individual GSE.
We do this for ID=["mouse", "rat", and "human"], and use the files "data/sample_lists/gse_${ID}.csv". 


`02_process/`
