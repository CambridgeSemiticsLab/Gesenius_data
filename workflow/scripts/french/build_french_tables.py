"""
Build french data.
"""
import sys
import json
from pathlib import Path

# NB that working directory when script is executed is 
# /workflow; because we have some utilities that we want
# to run from above directory, we need to append it to path
sys.path.append('scripts')
from build_tables import build_sample_tables

with open(snakemake.input.bhsa2french, 'r') as infile:
    bhsa2french = json.load(infile)
with open(snakemake.input.verses, 'r') as infile:
    frenchverses = json.load(infile)

def build_row(node):
    """Build a row for the french data"""
    french_word = bhsa2french.get(str(node), {})
    row_data = {'bhsa_node': node}
    row_data.update({
        'french_tense': french_word.get('french_tense', ''),
    })
    
    return row_data

def build_text_row(node):
    """Build a row for the french data"""
    french_word = bhsa2french.get(str(node), {})
    row_data = {'bhsa_node': node}
    ref = french_word.get('ref')
    row_data.update({
        'french': french_word.get('french', ''),
        'french_verse': frenchverses.get(ref, ''),
    })
    
    return row_data

rowmakers = [build_row, build_text_row]

build_sample_tables(
    rowmakers,
    snakemake.input.sample,
    snakemake.output
)
