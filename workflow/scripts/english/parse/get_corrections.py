"""
Retrieve necessary corrections from parsed tense data.
Apply already-corrected cases where necessary.
"""

import collections
import json
from pathlib import Path
import textwrap

verse2text = json.loads(Path(snakemake.input.verse2text).read_text())
sample = json.loads(Path(snakemake.input.sample).read_text())

trans_data = (
    ('esv', snakemake.input.esv, snakemake.params.esv_corr, snakemake.output.esv_todo),
    ('niv', snakemake.input.niv, snakemake.params.niv_corr, snakemake.output.niv_todo),
)

def export_todo(todos, data, outfile, trans):
    """Export a manual review sheet."""
    
    # cluster on verses
    verse2cases = collections.defaultdict(list)
    for case in todos:
        ref = str(tuple(data[case]['eng_ref']))
        verse2cases[ref].append(case)

    # assemble into document
    doc = f'n_cases: {len(todos)}\n\n'
    for ref, cases in verse2cases.items():
        text = '\n'.join(textwrap.wrap(verse2text[trans][ref], 80))
        ref_str = '{} {}:{}'.format(*eval(ref))
        doc += f'{ref_str}\n'
        doc += f'{text}\n'
        for case in cases:
            tense = data[case]['tense']
            span = data[case]['tense_span']
            doc += f'\t\t{case}\t{span}\t{tense}\n'
        doc += '\n'
    # export doc
    with open(outfile, 'w') as outfile:
        outfile.write(doc)

# process corrections to-do file for each translation
for trans, data_file, corr_file, todo_file in trans_data:

    trans_data = json.loads(Path(data_file).read_text())
    corr_data = json.loads(Path(corr_file).read_text())
    to_dos = []

    # go through nodes and figure out what needs to be added
    for node, data in trans_data.items():

        # only export on a sample-by-sample basis
        if int(node) not in sample:
            continue

        tense = data['tense']

        # trigger review for given example
        if '?' in tense:
            
            # trigger comparison with previous corrections
            if node in corr_data:
        
                # decide whether manual review is necessary
                orig_tense = corr_data[node]['orig_tense']

                # the underlying data has changed, delete old correction
                if orig_tense != tense:
                    del corr_data[node]
                    to_dos.append(node)
                else:
                    pass # keep correction as-is
            else:
                to_dos.append(node)

        # underlying tag has changed; delete old one
        elif node in corr_data:
            del corr_data[node]
        else:
            continue # no corrections necessary

    # export the document
    export_todo(to_dos, trans_data, todo_file, trans)

    # re-export the record
    with open(corr_file, 'w') as outfile:
        json.dump(corr_data, outfile)
