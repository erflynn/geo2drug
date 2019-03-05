
geo2drug
E Flynn
Last updated: 3/5/2019

The code in this repository is used for mapping GEO datasets to drugs (using DrugBank IDs or other).

Data mapping and prep:
* geo2pubmed_id.R: Takes geo data, maps to pubmed
* extract_sex_labels.R: Extract sex label data from drug mesh (run server-side)
* drugbank_synonyms.py: Extract drugbank data and synonyms
* mesh_to_ids.py: Downloads information from MeSH data
* xml_commands.sh: Commands for downloading XML for PMIDs
* parse_pubmed_xml.py: Parse pubmed XML

Data analysis:
* drug_gse_annot.R: Takes the results of the annotations and visualizes the results.
* examine_drug_freq.R: Looks at the frequency of drugs, visualizes the results
