README for download pipeline
- all of this is run server-side


This code serves to download a large set of gses and convert it to an expression matrix.
It also contains methods for sex-labeling / filtering (this is used for the silver standard).

Input:
 - list of GSEs (formatted as a file with one per line)
 - list of GSE/GSMs in deduplicated data
 - list of genes
 
1. To download:
   `download_gse_wrapper.sh`
    contains code that runs `download_geo_chunk.sh`
    which runs `downloadGEO.py`
    The end result is that all of these files are downloaded as series matrices to a specified directory.
    Then to get the failed GSEs:
    `grep "ERROR" *.out > error_gses.txt`


2. To check and convert to R objects:
    `submit_convert_to_obj.sh`
    which runs `seriesToObj.R`
    This checks the files for missing annotations/expression data, and then loads them
    To get the output counts:
    ```
      grep "unknown" logs/check_chunk* | grep -o "GSE*\w*" | sort | uniq > unknown_err.txt
      grep "expression data is missing" logs/check_chunk* | grep -o "GSE*\w*" | sort | uniq > exp_missing.txt
      grep "keys" logs/check_chunk* | grep -o "GSE*\w*" | sort | uniq > keys_missing.txt
      grep "30 percent" logs/check_chunk* | grep -o "GSE*\w*" | sort | uniq > warning_nas.txt
      grep "10k genes" logs/check_chunk* | grep -o "GSE*\w*" | sort | uniq > warning_genes.txt
    ```

3. To convert to an expression matrix based on a list of genes
   `Rscript create_exp_matrix.R <chunk_num>`
  Then combine these together using:
  `Rscript combine_exp_matrix.R`
  This also filters the expression data according to the de-duplicated list of GSEs/GSMs -AND- the list of datasets with warnings 


optional: To run sex labeling and filtering 
    `Rscript study_sex_label.R <chunk_num>`