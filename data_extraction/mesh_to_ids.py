# mesh_to_ids.py
# E Flynn
# 2/26/2019
# 
# Grabs mesh ID info using webscraper
#   we want CAS, synonym, UNII
#
# NOTE - this is SLOW and run locally because we are running headless chrome. do not re-run!
#  approx time 30s per ID
# TODO: speed up


from selenium import webdriver
from bs4 import BeautifulSoup
import time
import json

# general purpose function for stripping text from fields
def extractText(text_field):
	return text_field.strip().split("\n")[0]


# setup selenium with a headless chrome browser
browser = webdriver.Chrome()
options = webdriver.ChromeOptions()
options.add_argument('headless')
mesh_id_to_info = {}

# get list of mesh to download - NOTE: this is from geo2pubmed_id.R
with open("data/mesh_to_download_v3.txt", 'r') as f:
#with open("data/list_mesh.txt", 'r') as f:
	lines = f.readlines()
list_mesh_ids = [line.strip() for line in lines]

for mesh_id in list_mesh_ids:
	print mesh_id
	cas_number, registry_number, name = "", "", ""
	entry_list = []

  # grab the data
	url="https://meshb.nlm.nih.gov/record/ui?ui=%s" %(mesh_id)
	browser = webdriver.Chrome(chrome_options=options)
	browser.get(url)
	innerHTML = browser.execute_script("return document.body.innerHTML")
	time.sleep(20) # wait 20s for page to load
	innerHTML = browser.execute_script("return document.body.innerHTML")

	# parse the data
	soup = BeautifulSoup(innerHTML, "lxml")
	my_tab = soup.find(class_="tab-content").div.div # grab the table
	name = extractText( my_tab.dl.dd.text ) # grab the name
	
	# try to parse, but print error message if not present
	list_span = my_tab.findAll("span")
	try:
		for span in list_span:
			if "ng-if" in span.attrs:
				if span['ng-if']=="descriptor._generated.relatedRegistryNumberList[0]":
					cas_number = span.find("dd").text.strip()
			if "ng-show" in span.attrs:
				if span["ng-show"] == "descriptor.ConceptList.Concept[0].RegistryNumber.t":
					registry_number = span.find("dd").text.strip()
				if span["ng-show"] == "descriptor._generated.originalEntryTerms.length > 1":
					entries = span.findAll("dd")
					entry_list = [extractText(entry.text) for entry in entries]
	except:
		print "Info not present for %s" %(mesh_id)
	mesh_info = {"name" : name, "CAS" : cas_number, "registryNum" : registry_number, "synonyms" : entry_list}
	mesh_id_to_info[mesh_id]=mesh_info

# write out the info
with open("data/db_data/mesh_info4.json", 'w') as f:
	mesh_str = json.dumps(mesh_id_to_info)
	f.write(mesh_str)

# 57 of these have actual values
sum([value['CAS']=="" for key,value in mesh_id_to_info.iteritems()]) # 176
len(mesh_id_to_info) # 233

