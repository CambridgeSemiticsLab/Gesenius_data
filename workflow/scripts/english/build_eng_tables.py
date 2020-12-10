"""
Build data rows for tables using English
translation data.
"""

import json
import itertools
from pathlib import Path
from modify_tense_tag import split_TAM
from build_tables import build_sample_tables

# load up the translation data
bhsa2esv = json.loads(Path(snakemake.input.esv).read_text())
bhsa2niv = json.loads(Path(snakemake.input.niv).read_text())
trans2text = json.loads(Path(snakemake.input.txt).read_text())

def get_word(node, default={}):
    """Retrieve word data for both translations"""
    str_node = str(node)
    esv_word = bhsa2esv.get(str_node, {})
    niv_word = bhsa2niv.get(str_node, {})
    return [('esv', esv_word), ('niv', niv_word)]

def build_data_row(node):
    """Build a row for the translations tables"""
    
    transs = get_word(node)
    row_data = {'bhsa_node': node}

    # add TAM data from translations
    for trans, tdata in transs:
        ref_tuple = tuple(tdata.get('eng_ref', ''))
        row_data[f'{trans}_tags'] = tdata.get('tags', '')
        row_data[f'{trans}_VBtags'] = tdata.get('vb_tags', '')
        for tam_key, tam_data in split_TAM(tdata.get('TAM_cx', '')).items():
            tam_key = f'{trans}_{tam_key}'
            row_data[tam_key] = tam_data
        
    return row_data       

def build_text_row(node):
    """Build a row for translation tables with text."""
    transs = get_word(node)
    row_data = {'bhsa_node': node}
    for trans, tdata in transs:
        ref_tuple = tuple(tdata.get('eng_ref', ''))
        row_data[f'{trans}'] = tdata.get('words', '')
        row_data[f'{trans}_verse'] = trans2text[trans].get(str(ref_tuple), '')
        row_data[f'{trans}_TAMspan'] = tdata.get('TAM_span', '')

    return row_data

rowmakers = [build_data_row, build_text_row]

# make the tables
build_sample_tables(
    rowmakers, 
    snakemake.input.samples, 
    snakemake.output
)
