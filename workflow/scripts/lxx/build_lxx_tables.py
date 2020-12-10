"""
Build LXX data.
"""

import json
from pathlib import Path
from build_tables import build_sample_tables

bhsa2lxx = json.loads(Path(snakemake.input.lxx).read_text())

def build_row(node):
    """Build a row for the LXX data"""

    lxx_word = bhsa2lxx.get(str(node), {})
    row_data = {'bhsa_node': node}

    # add LXX data
    row_data.update({
        'lxx': lxx_word.get('utf8', ''),
        'lxx_tense': lxx_word.get('tense', ''),
        'lxx_voice': lxx_word.get('voice', ''),
        'lxx_mood': lxx_word.get('mood', ''),
        'lxx_person': lxx_word.get('person', ''),
        'lxx_number': lxx_word.get('number', ''),
    })
    
    # make lxx tense strings
    if lxx_word:
        row_data['lxx_tm'] = lxx_word['tense'] + ' ' + lxx_word['mood'] 
    else:
        row_data['lxx_tm'] = ''

    return row_data

rowmakers = [build_row]

build_sample_tables(
    rowmakers,
    snakemake.input.samples,
    snakemake.output
)
