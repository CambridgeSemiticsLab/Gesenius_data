"""
Retrieve necessary corrections from parsed tense data.
Apply already-corrected cases where necessary.
"""

import json
from pathlib import Path
from correction_files import build_todos

verse2text = json.loads(Path(snakemake.input.verse2text).read_text())
sample = json.loads(Path(snakemake.input.sample).read_text())

trans_data = (
    (snakemake.input.esv, snakemake.params.esv_corr, snakemake.output.esv_todo, verse2text['esv']),
    (snakemake.input.niv, snakemake.params.niv_corr, snakemake.output.niv_todo, verse2text['niv']),
)
# process corrections to-do file for each translation
for word_data_f, corr_file, todo_file, versetexts in trans_data:
    build_todos(sample, word_data_f, corr_file, todo_file, versetexts)
