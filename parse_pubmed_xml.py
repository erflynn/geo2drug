import xml.etree.ElementTree as ET
import json

pmid_to_mesh = dict()
for i in range(1, 14):
	f_string = "pubmed_%s.xml" %i
	tree = ET.ElementTree(file=f_string)
	root = tree.getroot()
	for pubmed_article in root:
		pmid = str(pubmed_article[0].findall("PMID")[0].text)
		mesh_headings = pubmed_article[0].findall("MeshHeadingList")
		if len(mesh_headings) == 0:
			continue
		for heading_list in mesh_headings[0]:
			mesh = []
			for mesh_heading in heading_list:
				if "UI" in mesh_heading.attrib:
					mesh_id = mesh_heading.attrib['UI']
					mesh.append(mesh_id)
		pmid_to_mesh[pmid] = mesh

with open("data/pmid_to_mesh.json", "w") as outfile:
	json.dump(pmid_to_mesh, outfile)
