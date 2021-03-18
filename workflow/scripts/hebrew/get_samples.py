import json
import pickle
import collections
from tf.fabric import Fabric
from pathlib import Path
from verb_form import get_verbform

# load basic BHSA data with Text-Fabric
TF = Fabric(snakemake.input, silent='deep')
bhsa = TF.load('pdp lex vt language', silent='deep')
F, L = bhsa.F, bhsa.L

# load GBI data for verb_form creation
with open(snakemake.input.bhsa2gbi, 'rb') as infile:
    bhsa2gbi = pickle.load(infile)

# loop through all verbs stored in the BHSA
# and select those forms specified by the wildcard
samples = []
for node in F.pdp.s('verb'):

    # skip non-hebrew words
    if F.language.v(node) != 'Hebrew':
        continue 

    verb_form = get_verbform(node, bhsa, bhsa2gbi)
    get_form = snakemake.wildcards.verb
    
    # handle cohortatives / jussives
    if get_form == 'yqtl' and verb_form in {'jussM', 'cohoM'}:
        samples.append(node)

    # handle all other matching verbs
    elif verb_form == get_form:
        samples.append(node)

# export the samples as {verb}.json
with open(snakemake.output.file, 'w') as outfile:
    json.dump(samples, outfile)
