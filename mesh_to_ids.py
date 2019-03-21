# grab mesh ID info
#
# we want CAS, synonym, UNII



from selenium import webdriver
from bs4 import BeautifulSoup
import time
import json

def extractText(text_field):
	return text_field.strip().split("\n")[0]

browser = webdriver.Chrome()
options = webdriver.ChromeOptions()
options.add_argument('headless')
mesh_id_to_info = {}

#with open("data/list_mesh.txt", 'r') as f:
with open("data/mesh_to_download.txt", 'r') as f:
	lines = f.readlines()

list_mesh_ids = [line.strip() for line in lines]

for mesh_id in list_mesh_ids:
	print mesh_id
	cas_number = "" 
	registry_number = ""
	entry_list = []
	name = ""
	url="https://meshb.nlm.nih.gov/record/ui?ui=%s" %(mesh_id)
	#html = lxml.html.fromstring(meshHTML)
	browser = webdriver.Chrome(chrome_options=options)
	browser.get(url)
	innerHTML = browser.execute_script("return document.body.innerHTML")
	time.sleep(20)
	innerHTML = browser.execute_script("return document.body.innerHTML")
	soup = BeautifulSoup(innerHTML, "lxml")
	#print soup.prettify()
	# grab the table
	my_tab = soup.find(class_="tab-content").div.div
	name = extractText( my_tab.dl.dd.text )
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
		# write this out

with open("data/mesh_info2.json", 'w') as f:
	mesh_str = json.dumps(mesh_id_to_info)
	f.write(mesh_str)

# 57 of these have actual values
sum([value['CAS']=="" for key,value in mesh_id_to_info.iteritems()]) # 176
len(mesh_id_to_info) # 233

### ---- old code, do not use ---- ###

import requests

def patient_request(url):
    i = 1
    result = requests.get(url, 
	headers={"User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.13; rv:62.0) Gecko/20100101 Firefox/62.0"})
    if result.status_code == 404:
        return(result)
    while(result.status_code != 200 and result.status_code != 404 and result.status_code != 403):
        time.sleep(10*i)
        print(result.status_code)
        result = requests.get(url)
        print("being patient for " + str(10*i) + " seconds...")
        i += 1
    return(result)
patient_request(url)

from requests_html import HTMLSession
session = HTMLSession()
r = session(url)
res = r.html.render()


# from SPARQLWrapper import SPARQLWrapper, JSON

# sparql = SPARQLWrapper("http://id.nlm.nih.gov/mesh/sparql")

# #https://meshb.nlm.nih.gov/record/ui?ui=D004074

# # looks like I should really run a SPARQL query for this :(

# sparql.setQuery(
# """
# PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
# PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
# PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
# PREFIX owl: <http://www.w3.org/2002/07/owl#>
# PREFIX meshv: <http://id.nlm.nih.gov/mesh/vocab#>
# PREFIX mesh: <http://id.nlm.nih.gov/mesh/>
# PREFIX mesh2015: <http://id.nlm.nih.gov/mesh/2015/>
# PREFIX mesh2016: <http://id.nlm.nih.gov/mesh/2016/>
# PREFIX mesh2017: <http://id.nlm.nih.gov/mesh/2017/>

# SELECT *
# FROM <http://id.nlm.nih.gov/mesh>

# WHERE {
# ?s ?p meshv:D004074 . }



# """)
# sparql.setReturnFormat(JSON)
# results = sparql.query().convert()
# print results