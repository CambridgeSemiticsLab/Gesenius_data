"""
Ingest the French dataset and output 
JSON files that can be used for subsequent
parsing and processing.
"""

import re
import csv
import json
from tf.fabric import Fabric
from Levenshtein import distance as levdist

eng_order = [
    'Genesis',
    'Exodus',
    'Leviticus',
    'Numbers',
    'Deuteronomy',
    'Joshua',
    'Judges',
    'Ruth',
    '1_Samuel',
    '2_Samuel',
    '1_Kings',
    '2_Kings',
    '1_Chronicles',
    '2_Chronicles',
    'Ezra',
    'Nehemiah',
    'Esther',
    'Job',
    'Psalms',
    'Proverbs',
    'Ecclesiastes',
    'Song_of_songs',
    'Isaiah',
    'Jeremiah',
    'Lamentations',
    'Ezekiel',
    'Daniel',
    'Hosea',
    'Joel',
    'Amos',
    'Obadiah',
    'Jonah',
    'Micah',
    'Nahum',
    'Habakkuk',
    'Zephaniah',
    'Haggai',
    'Zechariah',
    'Malachi',
]
int2book = {i+1: book for i, book in enumerate(eng_order)}

ref_re = re.compile(r'(\d\d\d)(\d\d\d)(\d\d\d)(\d\d)(\d\d\d)')

def parse_refstring(string):
    """Parse a refstring from UBS.
    
    String consists of:
        BBBCCCVVVSSWWW
    where:
        B = Book, C = chapter, V = Verse, 
        S = segment (can be ignored), W = Word
    """
    data = ref_re.match(string).groups()
    return [int(s) for s in data]

class BhsaWord:
    def __init__(self, node, dist):
        """Store BHSA word node and distance from target word."""
        self.node = node
        self.dist = dist

def ingest_french(paths):
    """Match the French data to our dataset."""

    # load the French dataset
    with open(paths['source'], 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        french_data = list(reader)
    
    # load the BHSA Hebrew data for matching the Hebrew text
    TF = Fabric(locations=paths['bhsa'])
    API = TF.load('g_word_utf8')
    F, T, L = API.F, API.T, API.L
    
    # match the Hebrew verbs in the French data with the 
    # Hebrew verbs in BHSA
    # we treat the ref strings as unique ID's
    # we use 2 dicts; one to hold ID 2 BHSA node mappings
    # another to hold the IDs 2 french data
    french2bhsa = {}
    french2data = {}
    frenchverses = {}
    
    for row in french_data:
        
        # parse French data
        wid = row[0]
        hb_txt, hb_lex, hb_tag, hb_prev = row[1:5]
        fr_words, fr_verse = row[5:7]
        bk, ch, vs, sg, wnum = parse_refstring(wid)
        french2data[wid] = {
            'wid': wid,
            'hebrew': hb_txt,
            'hebrew_parse': hb_tag,
            'french': fr_words,
        }

        # look up BHSA data and get the verse node
        tf_book = int2book[bk]
        vrs_node = T.nodeFromSection((tf_book, ch, vs))
        if vrs_node is None:
            raise Exception((tf_book, ch, vs), wid, hb_txt)

        # save the French verse text
        ref_string = str((tf_book, ch, vs))
        frenchverses[ref_string] = fr_verse        
        french2data[wid]['ref'] = ref_string
       
        # get the closest matching word from the verse;
        # NB we iterate over the verse words in reversed order
        # so that if there are 2+ words with equivalent distances,
        # we always end on the one that is first in the verse;
        # the match is then added to a set so that it is not 
        # available for subsequent matches
        french2bhsa[wid] = BhsaWord(0, float('inf')) # initialize with dummy 
        matched = set() 
        for word_node in reversed(L.d(vrs_node, 'word')):
            if word_node in matched:
                continue
            bhsa_txt = T.text(word_node)
            dist = levdist(bhsa_txt, hb_txt)
            if french2bhsa[wid].dist > dist:
                french2bhsa[wid] = BhsaWord(word_node, dist)
        matched.add(french2bhsa[wid].node) 
                
    # iterate over both french dicts and assemble
    # into one BHSA dict
    bhsa2french = {}
    for wid, bhsa_word in french2bhsa.items():
        bhsa_node = bhsa_word.node 
        if bhsa_node != 0:
            bhsa2french[bhsa_node] = french2data[wid]

    # the linking is complete
    with open(paths['out'], 'w') as outfile:
        json.dump(bhsa2french, outfile, indent=2, ensure_ascii=False)

    with open(paths['out_verses'], 'w') as outfile:
        json.dump(frenchverses, outfile, indent=2, ensure_ascii=False)

# -- Execution with Snakemake -- 

paths = {
    'source': snakemake.input.source,
    'bhsa': snakemake.input.bhsa,
    'out': snakemake.output.bhsa2french,
    'out_verses': snakemake.output.verses, 
}

ingest_french(paths)
