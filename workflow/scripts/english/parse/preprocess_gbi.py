"""
Preprocess GBI dataset to make various mappings.
"""
import collections
import pickle
import json
from pathlib import Path
from gbi_functions import id2ref

gbi_niv = json.loads(Path(snakemake.input.niv).read_text())
gbi_esv = json.loads(Path(snakemake.input.esv).read_text())

sources = (('niv', gbi_niv), ('esv', gbi_esv))
word_data = collections.defaultdict(lambda: collections.defaultdict(list))
verse2words = collections.defaultdict(lambda: collections.defaultdict(list)) 
linkbyid = collections.defaultdict(list) # list of 2-tuples, each containing word IDs
id2link = collections.defaultdict(dict) # select a link based on a single ID 

# follow the layout of the GBI data as it's been provided:
# i.e. dictionaries of verse data which are stored in lists;
# along the way harvest all of the data we need in a more straightforward format
for name, source in sources:
    for verse in source:
        
        # unpack words for processing
        trans_words =  verse['translation']['words']
        manu_words = verse['manuscript']['words']
        
        # map translation word data
        for w in trans_words:
            ref_tuple = id2ref(w['id'], 'translation')
            verse2words[name][ref_tuple].append(w['id'])
            word_data[name][w['id']] = w
        
        # map WLC word data
        # arbitrarily use the copy stored under NIV
        if name == 'niv':
            for w in manu_words:
                ref_tuple = id2ref(w['id'])
                verse2words['wlc'][ref_tuple].append(w['id'])
                word_data['wlc'][w['id']] = w
                
        # map links to word ids
        # the alignment data just contains indices pointing
        # to the various lists, so these have to be used to 
        # identify the specific word in question
        for wlc_indices, trans_indices in verse['links']:
            wlc_ids = tuple(manu_words[i]['id'] for i in wlc_indices)
            trans_ids = tuple(sorted(trans_words[i]['id'] for i in trans_indices))
            linkbyid[name].append((wlc_ids, trans_ids))
            for wid in wlc_ids:
                id2link[name][wid] = trans_ids

# dump the files
data2path = [
    (word_data, snakemake.output.word_data),
    (verse2words, snakemake.output.verse2words),
    (linkbyid, snakemake.output.linkbyid),
    (id2link, snakemake.output.id2link),
]

for data, path in data2path:
    data = dict(data)
    with open(path, 'wb') as outfile:
        pickle.dump(data, outfile)
