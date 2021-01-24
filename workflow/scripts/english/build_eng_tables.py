"""
Build data rows for tables using English
translation data.
"""

import re
import sys
import json
import pickle
import itertools
from pathlib import Path
from modify_tense_tag import process_TAM

# NB that working directory when script is executed is 
# /workflow; because we have some utilities that we want
# to run from above directory, we need to append it to path
sys.path.append('scripts')
from build_tables import build_sample_tables

# load up the translation data
bhsa2esv = json.loads(Path(snakemake.input.esv).read_text())
bhsa2niv = json.loads(Path(snakemake.input.niv).read_text())
with open(snakemake.input.parsedverses, 'rb') as infile:
    trans2parse = pickle.load(infile)

#trans2text = json.loads(Path(snakemake.input.txt).read_text())

# load up corrections data
TAM_corrections = {
    'niv': json.loads(Path(snakemake.input.niv_corr).read_text()),
    'esv': json.loads(Path(snakemake.input.esv_corr).read_text()),
}

def get_word(node, default={}):
    """Retrieve word data for both translations"""
    str_node = str(node)
    esv_word = bhsa2esv.get(str_node, {})
    niv_word = bhsa2niv.get(str_node, {})
    return [('esv', esv_word), ('niv', niv_word)]

def compare_transs(tdata):
    """Compare translation words and return data."""

    data = {
        'eng_fullparse': 0,
        'eng_TAM': '',
        'eng_TAMsimp': '',
        'eng_agree': 0,
        'eng_simp_agree': 0,
    }

    # add features based on both translations
    niv_tam, esv_tam = tdata['niv_TAM'], tdata['esv_TAM']
    niv_tsimp, esv_tsimp = tdata['niv_TAMsimp'], tdata['esv_TAMsimp']

    if niv_tam and esv_tam:

        data['eng_fullparse'] = 1
        
        # record interchanges of TAM categories between ESV and NIV where the disagree
        if niv_tam != esv_tam:
            eng_tam = sorted([niv_tam, esv_tam])
            data['eng_TAM'] = ' ~ '.join(eng_tam)
        else:
            data['eng_TAM'] = niv_tam
            data['eng_agree'] = 1

        # repeat for simplified TAM tags
        if niv_tsimp != esv_tsimp:
            eng_tam = sorted([niv_tsimp, esv_tsimp])
            data['eng_TAMsimp'] = ' ~ '.join(eng_tam)
        else:
            data['eng_TAMsimp'] = niv_tsimp
            data['eng_simp_agree'] = 1

    return data

find_is = re.compile(r'\b[iI]s\b')

def build_data_row(node):
    """Build a row for the translations tables"""
    
    trans2word_data = get_word(node)
    row_data = {'bhsa_node': node}

    # add TAM data from translations
    for transl, word_data in trans2word_data:
        ref_tuple = tuple(word_data.get('eng_ref', ''))
        row_data[f'{transl}_tags'] = word_data.get('tags', '')
        row_data[f'{transl}_VBtags'] = word_data.get('vb_tags', '')

        # process TAM / tense tags
        raw_tam = word_data.get('tense', '')
        str_node = str(node)
        corr_tam = TAM_corrections[transl].get(str_node, {}).get(raw_tam, raw_tam) # apply corrections
        tam_data = process_TAM(corr_tam)
        for key, value in tam_data.items():
            key = f'{transl}_{key}'
            row_data[key] = value
        
    # add features based on comparisons between translations
    row_data.update(compare_transs(row_data))

    return row_data       

def build_text_row(node):
    """Build a row for translation tables with text."""
    trans2word_data = get_word(node)
    row_data = {'bhsa_node': node}
    for transl, word_data in trans2word_data:
        ref_tuple = tuple(word_data.get('eng_ref', ''))
        row_data[f'{transl}'] = word_data.get('words', '')
        row_data[f'{transl}_TAMspan'] = word_data.get('tense_span', '')

        # process verse text
        parsed_verse = trans2parse[transl].get(ref_tuple, '')
        verse_text = str(parsed_verse)
        row_data[f'{transl}_verse'] = verse_text

        # process sentence text
        sent_i = word_data.get('sentence_i', None)
        if sent_i != None:
            verse_sentences = list(parsed_verse.sents)
            sent_parse = verse_sentences[sent_i]
            sent_text = str(sent_parse)
        else:
            sent_text = ''
        row_data[f'{transl}_sent'] = sent_text

    # detect 'is' as a proxy for stative translations
    row_data['esv_is'] = 1 if find_is.search(row_data.get('esv', '')) else 0
    row_data['niv_is'] = 1 if find_is.search(row_data.get('niv', '')) else 0

    return row_data

rowmakers = [build_data_row, build_text_row]

# make the tables
build_sample_tables(
    rowmakers, 
    snakemake.input.sample, 
    snakemake.output
)
