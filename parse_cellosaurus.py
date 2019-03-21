# parse_cellosaurus.py
# E Flynn
# 3/7/2019
#
# Parse the cellosaurus XML
#
# GOAL:
# - extract

import xml.etree.ElementTree as ET

cellFile = "data/cellosaurus.xml"
tree = ET.parse(cellFile)
root = tree.getroot() 

cl_list = root.find("cell-line-list")
cls = cl_list.findall("cell-line") # 110,948
cl_dict = {}
for cl in cls:
	cl_sex = ""
	cl_category = ""
	cl_age=""
	cl_sex = cl.get("sex")
	cl_category = cl.get("category")
	cl_age = cl.get("age")
	accs = cl.find("accession-list").findall("accession")
	primary_acc = []
	secondary_acc = []
	for acc in accs:
		acc_num = acc.text
		acc_type = acc.get("type") # primary or secondary
		if acc_type=="primary":
			primary_acc.append(acc_num)
		else:
			secondary_acc.append(acc_num)
	acc_info = {"primary" : primary_acc, "secondary" : secondary_acc}
	names = cl.find("name-list").findall("name")
	identifier = ""
	synonyms = []
	for name in names:
		name_text = name.text
		name_type = name.get("type") # identifier or synonym
		if name_type == "identifier":
			if identifier != "":
				print "Error multiple identifiers"
				print identifier, name_text
			else:
				identifier = name_text
		else:
			synonyms.append(name_text)
	if identifier == "":
		print "Error no identifier for cell line"
		continue
	# extract the disease of origin
	dz_list = cl.find("disease-list")
	dzs=[]
	if dz_list is not None:
		dzs =[dz.text for dz in dz_list.findall("cv-term")]
	# extract database mappings - Gene Expression databases
	xref_l = cl.find("xref-list")
	gex_list = []
	mesh_acc = ""
	if xref_l is not None:
		xrefs = xref_l.findall("xref")
		for xref in xrefs:
			cat = xref.get("category")
			if cat == "Gene expression databases":
				gex_acc = xref.get("accession")
				if gex_acc is not None:
					gex_list.append(gex_acc)
			if cat == "Ontologies":
				if xref.get("database")=="MeSH":
					mesh_acc = xref.get("accession")
	cl_dict[identifier] = {"dz":dzs, "sex": cl_sex, "category": cl_category, "age" : cl_age, \
	 "synonyms": synonyms, "accession":acc_info, "mesh":mesh_acc, "gex": gex_list}


# write out the json
with open("data/cellosaurus.json", 'w') as f:
	cl_str = json.dumps(cl_dict)
	f.write(cl_str)

# names/synonyms <-- what is the complete list?



# what is the sex breakdown of the data???