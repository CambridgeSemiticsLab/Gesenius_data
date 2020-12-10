import json
import collections
from tf.fabric import Fabric
from pathlib import Path
from verb_form import get_verbform

# load basic BHSA data with Text-Fabric
TF = Fabric(snakemake.input, silent='deep')
bhsa = TF.load('pdp lex vt language', silent='deep')
F, L = bhsa.F, bhsa.L

samples = collections.defaultdict(list)

# loop through all verbs stored in the BHSA
# and store those that are specified for collection in config
for verb in F.pdp.s('verb'):

    # skip non-hebrew words
    if F.language.v(verb) != 'Hebrew':
        continue 

    verb_form = get_verbform(verb, bhsa)
    if verb_form in snakemake.config['verb_forms']:
        samples[verb_form].append(verb)

# export the samples as {verb}.json
for verb_form in samples:
    filepath = Path(snakemake.output[0]).joinpath(f'{verb_form}.json')
    with open(filepath, 'w') as outfile:
        json.dump(samples[verb_form], outfile)
