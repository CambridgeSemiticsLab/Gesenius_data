"""
Here we match up 1 bhsa verb node to 1 gbi verb ID using TF and GBI word features.
"""

import json
import pickle
from pathlib import Path
from tf.fabric import Fabric

# -- SET UP --

# fire up Text-Fabric with BHSA data
TF = Fabric(snakemake.input['tf_mods'], silent='deep')
bhsa = TF.load('pdp', silent='deep')
F, E, T, L = bhsa.F, bhsa.E, bhsa.T, bhsa.L

# load the bhsa2gbi mapping
bhsa2gbi_file = Path(snakemake.input.bhsa2gbi)
bhsa2gbi = json.loads(bhsa2gbi_file.read_text())

# load preprocessed GBI data
with open(snakemake.input.word_data, 'rb') as infile:
    word_data = pickle.load(infile)

# -- MAKE THE ALIGNMENTS --

# track where BHSA and WLC disagree on the classification of a verb
no_match = []
verb_bhsa2gbi = {}

# find 1-to-1 matches of BHSA and WLC verbs
for bhsa_nodes, gbi_ids in bhsa2gbi:
    
    # filter out non-verbs from the links
    bhsa_verbs = [w for w in bhsa_nodes if F.pdp.v(w) == 'verb'] 
    wlc_verbs = [w for w in gbi_ids if word_data['wlc'][w]['pos'] == 'verb']
    data = (T.text(bhsa_nodes), T.sectionFromNode(bhsa_nodes[0]), bhsa_nodes, gbi_ids) # track null matches
    
    # one case, Jer 51:3, has a double verb mapping caused by 
    # ידרך ידרך, which BHSA maps to a single word node, and gbi 
    # keeps as 2 words; we disambig that here and keep only 
    # first gbi word
    if bhsa_verbs and bhsa_verbs[0] == 262780:
         wlc_verbs = wlc_verbs[:1]
    
    # skip non-verbal contexts
    if not bhsa_verbs + wlc_verbs:
        continue
    
    # track disagreements between 2 sources
    elif (bhsa_verbs and not wlc_verbs) or (wlc_verbs and not bhsa_verbs):
        no_match.append(data)
    
    # store alignment
    elif len(bhsa_verbs) == 1 and len(wlc_verbs) == 1:
        bhsa_verb, wlc_verb = bhsa_verbs[0], wlc_verbs[0]
        verb_bhsa2gbi[bhsa_verb] = word_data['wlc'][wlc_verb]
    
    # or there's a problem...
    else:
        raise Exception(f'Misalignment at {data}')
        
print(len(no_match), 'verbs are unmatched')
print(len(verb_bhsa2gbi), 'verbs aligned')

# dump the files
with open(snakemake.output.matches, 'w') as outfile:
    json.dump(verb_bhsa2gbi, outfile, ensure_ascii=False, indent=2)

no_match_data = '\n'.join(str(nm) for nm in no_match)
Path(snakemake.output.no_match).write_text(no_match_data)
