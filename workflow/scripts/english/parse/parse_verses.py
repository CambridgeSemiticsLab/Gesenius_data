# parse a dict of verses with Spacy and output the results

import collections
import pickle
import json
import spacy
from tf.core.timestamp import Timestamp

ts = Timestamp() # for reporting

# set up the Spacy processor as well as some customized attributes
nlp = spacy.load('en_core_web_sm')

def parse_verse(verse_tuple, translation, verse2words, word_data):
    """Parse translation verse with Spacy."""
    word_ids = verse2words[translation][verse_tuple]
    words = [word_data[translation][w] for w in word_ids]
    text = ' '.join(w['text'] for w in words)
    text = text.replace(';', '.') # remove ambiguity for parser
    parsed_doc = nlp(text) # magic happens here
    return parsed_doc

def parse_verses(transdict, verse2words, word_data):
    """Iterate through all verses and parse them.
    
    Args:
        versedict: dict with structure of e.g. {'niv': set(('Genesis', 1, 1)...}}
    Returns:
        dict w/ structure of e.g. {'niv': {('Genesis', 1, 1): Spacy.Doc}}
    """
    
    parsed_verses = collections.defaultdict(dict)
    
    ts.indent(reset=True)
    ts.info(f'Parsing translations...')
    
    for translation, verse_set in transdict.items():
    
        # it takes a long time so we time it
        ts.indent(1, reset=True)
        ts.info(f'Beginning {translation}...')
        ts.indent(2, reset=True)
    
        # parse the verse and put Spacy.Doc in a dict
        for i, ref_tuple in enumerate(verse_set):
            
            if i % 5000 == 0 and i != 0:
                ts.info(f'done with verse {i}')
                
            parsed_verses[translation][ref_tuple] = parse_verse(ref_tuple, translation, verse2words, word_data)
            
        ts.indent(1)
        ts.info('done!')
    
    return parsed_verses

# RUN THE PARSINGS:
def open_pickle(path):
    with open(path, 'rb') as infile:
        return pickle.load(infile)

versestoparse = open_pickle(snakemake.input.versestoparse) 
verse2words = open_pickle(snakemake.input.verse2words)
word_data = open_pickle(snakemake.input.word_data)

# export parse data
parsed_verses = parse_verses(versestoparse, verse2words, word_data)
with open(snakemake.output.parsedverses, 'wb') as outfile:
    pickle.dump(parsed_verses, outfile)

# export plain-text verse data
verse_text_path = snakemake.output.verse2text
verse_text_str = collections.defaultdict(dict)
for version, verses in parsed_verses.items():
    for ref_tuple, text in verses.items():
        verse_text_str[version][str(ref_tuple)] = str(text)        

with open(verse_text_path, 'w') as outfile:
    json.dump(verse_text_str, outfile, ensure_ascii=False)
