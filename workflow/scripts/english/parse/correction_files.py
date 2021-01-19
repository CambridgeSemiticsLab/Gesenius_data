"""
Open manual corrections file and apply corrections.
"""

import re
import json
import collections
import textwrap
from pathlib import Path

def export_todo(todos, word_data, versetexts, outfile):
    """Export a manual review sheet."""
   
    # cluster on verses
    verse2cases = collections.defaultdict(list)
    for case in todos:
        ref = str(tuple(word_data[case]['eng_ref']))
        verse2cases[ref].append(case)

    # assemble into document
    meta = {'n_cases': len(todos), 'complete': False}
    doc = f'{meta}\n\n'
    for ref, cases in verse2cases.items():
        text = versetexts[ref] 
        ref_str = '{} {}:{}'.format(*eval(ref))
        doc += f'{ref_str}\n'
        cases_text = ''
        span_covered = set()
        for case in cases:
            tense = word_data[case]['tense']
            span = word_data[case]['tense_span']
            gloss = word_data[case]['words']
            cases_text += f'\t\t{case}|{tense}|{span}|({gloss})\n'
            if span not in span_covered:
                text = re.sub(f'([A-Za-z ]*)(\s{span}\s)([A-Za-z ]*)([.?! ’,]*)', '\n > \g<1> [\g<2>]  \g<3>\g<4>\n', text)
                span_covered.add(span)
        lines = text.splitlines()
        lines = ['\n'.join(textwrap.wrap(line, 80)) for line in lines]
        text = '\n'.join(lines).replace('\n\n', '\n')
        doc += f'{cases_text}{text}'.replace('\n\n', '\n') + '\n\n'
    # export doc
    with open(outfile, 'w') as outfile:
        outfile.write(doc.strip())

def build_todos(sample, word_data_f, corr_file, todo_file, versetexts):
    """Gather cases to export to a review sheet."""

    word_data = json.loads(Path(word_data_f).read_text())
    corr_data = json.loads(Path(corr_file).read_text())
    to_dos = []

    # go through nodes and figure out what needs to be added
    for node, data in word_data.items():

        # only export on a sample-by-sample basis
        if int(node) not in sample:
            continue

        tense = data['tense']

        # trigger review for given example
        if ('?' in tense) and (tense not in corr_data.get(node, {})):
            to_dos.append(node)

    # export the to-do document
    export_todo(to_dos, word_data, versetexts, todo_file)

class TodoReader:
    def __init__(self, path):
        doc = Path(path).read_text()
        try:
            header, data = re.split('\n\n', doc, 1)
            self.header = eval(header)
            self.data = re.findall('\t\t(\d+)\|(.+)\|(.+)\|(.+)', data)
        except ValueError:
            self.header = eval(doc.strip())
            self.data = []

def apply_corrections(corr_file, todo_file, word_data_f, versetexts, report_f):
    """Read in a corrections record file and to-do and determine how to save changes.""" 

    # load data for processing
    word_data = json.loads(Path(word_data_f).read_text())
    corr_record = json.loads(Path(corr_file).read_text())
    corrs = TodoReader(todo_file)
    file_complete = corrs.header['complete']
    to_dos = []
    
    # iterate through all corrections and enact changes as necessary
    for node, corr_tense, span, gloss in corrs.data:

        original_tense = word_data[node]['tense']
        
        # keep original value when file complete
        if corr_tense == original_tense and file_complete:
            corr_entry = corr_record.setdefault(node, {}) 
            corr_entry[original_tense] = corr_tense.replace('?', '')

        # remap a new value
        elif corr_tense != original_tense:
            corr_entry = corr_record.setdefault(node, {}) 
            corr_entry[original_tense] = corr_tense

        # place the to-do back in the to-do set
        else:
            to_dos.append(node)

    # re-export to-do list
    export_todo(to_dos, word_data, versetexts, todo_file)

    # re-export the corrections record
    with open(corr_file, 'w') as outfile:
        json.dump(corr_record, outfile, indent=2)

    # export dummy file to give Snakemake something to chew on
    with open(report_f, 'w') as outfile:
        outfile.write('')
