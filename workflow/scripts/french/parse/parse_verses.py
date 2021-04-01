"""
Use Spacy 2.x to parse French verses
into a bunch of Doc objects that can
later be matched with tense rules.
"""

import json
import pickle
import spacy

def parse_verses(paths):
    """Parse French verses."""

    # open French verse data 
    with open(paths['verses'], 'r') as infile:
        frenchverses = json.load(infile)

    # load Spacy model
    nlp = spacy.load('fr_core_news_lg')

    # parse the verses
    print(f'Parsing {len(frenchverses)} verses...')
    ref2parse = {}
    for ref, verse in frenchverses.items():
        parse = nlp(verse)
        ref2parse[ref] = parse
    print(f'done! All verses parsed.')

    # dump the parsed objects
    with open(paths['parses'], 'wb') as outfile:
        pickle.dump(ref2parse, outfile)

# -- Run with Snakemake --
paths = {
    'verses': snakemake.input.verses,
    'parses': snakemake.output.parses
}

parse_verses(paths)
