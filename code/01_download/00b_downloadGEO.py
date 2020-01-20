
import wrenlab.util.aspera

import sys
#from joblib import Parallel, delayed

gse_file = sys.argv[1]
dir_id = sys.argv[2]
OUTDIR="gses_%s/" %dir_id

def group_prefix(accession): # from wrenlab
    if len(accession) <= 6:
        group = "GPLnnn"
    else:
        group = "%snnn" %(accession[:(len(accession)-3)])
#    elif len(accession) == 7:
#        group = "{}nnn".format(accession[:4])
#    else:
#        group = "{}nnn".format(accession[:5])
    return group



def grab_gse(gse, c):
    try:
        group = group_prefix(gse)
        uri = "/geo/series/%s/%s/matrix/" %(group, gse)
        c.download(uri, OUTDIR)
        print("Downloaded %s" %gse)
    except Exception as e:
        print("ERROR with %s" %gse)
        print(e)

with open(gse_file, 'r') as f:
    gse_list = [gse.strip() for gse in f.readlines()]

c = wrenlab.util.aspera.Client(host="ftp.ncbi.nlm.nih.gov")
for gse in gse_list:
    grab_gse(gse, c)
