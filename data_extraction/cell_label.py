

def cell_label(cell_dict, text):


import re
import pandas as pd


gse_text = pd.read_csv("data/gse_text.txt", sep="\t")
gse_dict = gse_text.set_index('gse').T.to_dict()

texts = gse_text['str'].tolist()

cell_df = pd.read_csv("data/multiword_cell.txt", sep="\t").drop_duplicates()
cell_df['cl'] = cell_df['cl'].astype(str)
cell_dict = cell_df.T.to_dict().values()

reverse = {d['cl']:d['accession'] for d in sorted(cell_dict, key=lambda x: x['cl'])}
re_phrases = re.compile('({})'.format('|'.join(d['cl'] for d in cell_dict)))

list_matches = []
for row in gse_text.itertuples():
    text = row[2]
    matches = [(reverse[re_whitespace.sub(' ', match.group(1))]) for match in re_phrases.finditer(text)]
    matches_set = set(matches)
    if (len(matches_set) > 0):
        list_matches.append((row[1], matches_set))

# save this list of matches