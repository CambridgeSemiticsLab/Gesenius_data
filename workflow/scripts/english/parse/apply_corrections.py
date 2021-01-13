import json
from pathlib import Path
from correction_files import apply_corrections

verse2text = json.loads(Path(snakemake.input.verse2text).read_text())

trans_data = [
    (snakemake.params.esv_corr, snakemake.input.esv_todo, snakemake.output.esv_todo, snakemake.input.esv, verse2text['esv']),
    (snakemake.params.niv_corr, snakemake.input.niv_todo, snakemake.output.niv_todo, snakemake.input.niv, verse2text['niv']),
]

for corr_file, todo_infile, todo_outfile, word_data_f, versetexts in trans_data:
    apply_corrections(corr_file, todo_infile, todo_outfile, word_data_f, versetexts)
