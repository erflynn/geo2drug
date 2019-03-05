# drugbank_synonyms.py
# E Flynn
# 2/27/2019
#
# Get drugbank ID, name, synonym, and ATCCC code from the drugbank XML.
#
# TODO:
# - may want to parse other OTC, FDA approved fields


import xml.etree.ElementTree as ET
import json

drugBankFile = "../../../spring2018/drug-target/DrugBank.xml"
tree = ET.parse(drugBankFile) # slowest step
root = tree.getroot()
treeDrugs = root.findall("{http://www.drugbank.ca}drug")

drug_data = []

for drug in treeDrugs:
	dbID = ""
	drugBankIDs = drug.findall("{http://www.drugbank.ca}drugbank-id")
	for drugBankID in drugBankIDs:
		if drugBankID.get("primary")=="true":
			dbID = drugBankID.text
	nameField = drug.find("{http://www.drugbank.ca}name")
	if nameField is not None:
		drugName = nameField.text
	else:
		drugName=""
	cas = drug.find("{http://www.drugbank.ca}cas-number")
	if cas is not None:
		drugCas = cas.text
	else:
		drugCas = ""
	unii = drug.find("{http://www.drugbank.ca}unii")
	if unii is not None:
		drugUnii = unii.text
	else:
		drugUnii=""
	synonyms = drug.findall("{http://www.drugbank.ca}synonyms")[0].findall("{http://www.drugbank.ca}synonym")
	syn_list = []
	for synonym in synonyms:
		if synonym is not None:
			syn_list.append(synonym.text)
	chebi_id = ""
	external_ids = drug.findall("{http://www.drugbank.ca}external-identifiers")[0].findall("{http://www.drugbank.ca}external-identifier")
	for external_id in external_ids:
		resource = external_id.find("{http://www.drugbank.ca}resource")
		if resource is not None:
			if resource.text=="ChEBI":
				chebi_id = external_id.find("{http://www.drugbank.ca}identifier").text
	categories = drug.findall("{http://www.drugbank.ca}categories")[0].findall("{http://www.drugbank.ca}category")
	category_list = []
	for category in categories:
		# grab the category, meshid 
		cat=category.find("{http://www.drugbank.ca}category")
		if cat is not None:
			cat_label = cat.text
		else:
			cat_label = ""
		mesh = category.find("{http://www.drugbank.ca}mesh-id")
		if mesh is not None:
			meshid = mesh.text
		else:
			meshid=""
		category_list.append((cat_label, meshid))
	atc_code = drug.findall("{http://www.drugbank.ca}atc-codes")[0].find("{http://www.drugbank.ca}atc-code")
	if atc_code is not None:
		atc_id = atc_code.get("code")
	else:
		atc_id=""
	drug_info = {"dbID":dbID, "name":drugName, "chebi": chebi_id, "cas":drugCas, "unii":drugUnii, "synonyms" :syn_list, "categories":category_list, "ATC":atc_id}
	drug_data.append(drug_info)

# json dump the whole thing
with open("data/drugbank_info.json", 'w') as f:
	drug_str = json.dumps(drug_data)
	f.write(drug_str)
