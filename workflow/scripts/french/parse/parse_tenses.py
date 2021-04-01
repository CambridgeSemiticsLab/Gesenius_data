"""
Use Spacy Matcher (2.x) to match tenses of 
interest in French.
"""

import re
import json
import pickle
import spacy
from spacy.matcher import Matcher
from tense_rules import rules

def match_tenses(paths):
    """Match tenses from Spacy Doc objects.

    French verses have already been parsed into 
    Spacy Doc objects. We can apply a Matcher 
    along with some custom-written rules to match
    verb tenses of interest. We do that here.
    """

    # load the parsed Spacy Docs, linked by verse ref
    with open(paths['parses'], 'rb') as infile:
        ref2parse = pickle.load(infile)

    # load the bhsa2french aligned data
    with open(paths['bhsa2french'], 'r') as infile:
        bhsa2french = json.load(infile)

    # initialize the Spacy model and the Matcher with the custom rules
    nlp = spacy.load('fr_core_news_lg')
    matcher = Matcher(nlp.vocab)
    for tag, ruleset in rules:
        matcher.add(tag, None, ruleset)

    # for all verses, get match items and store by 
    # the verse reference; match items are stored 
    # in a list
    ref2matches = {}
    for ref, spacydoc in ref2parse.items():
        ref2matches[ref] = matcher(spacydoc)

    # link the bhsa2french data to the accepted matches;
    # we do this by accepting the first matching span in the 
    # matches list and then popping it from the list;
    # this ensures that any duplicate matches will not interfere
    for bhsa_node, data in bhsa2french.items():
        french = data['french']
        ref = data['ref']
        
        # skip null matches
        matches = ref2matches[ref]
        if not matches:
            continue

        # add any matches to the bhsa2french data
        # NB: iterate over generator so we 
        # can pop items off from the original list
        for match in (m for m in matches):
            str_id, start, end = match
            span = str(ref2parse[ref][start:end])
            span_re = fr'\b{span}\b' 
            french_re = fr'\b{french}\b'
            if re.search(span_re, french) or re.search(french_re, span):
                data['french_tense'] = nlp.vocab.strings[str_id]
                data['tense_span'] = span
                matches.remove(match) 
                break

    # re-export the modified bhsa2french data
    with open(paths['new_bhsa2french'], 'w') as outfile:
        json.dump(bhsa2french, outfile, indent=2, ensure_ascii=False)

# -- run with Snakemake --

paths = {
    'parses': snakemake.input.parses,
    'bhsa2french': snakemake.input.bhsa2french,
    'new_bhsa2french': snakemake.output.bhsa2french,
}

match_tenses(paths)
