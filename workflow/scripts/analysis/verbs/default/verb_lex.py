import sys
import pandas as pd

# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
both_df = DfLoad.eng_both()

run_analyses([
    {
        'name': 'verb_lex',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['lex', 'stem'],
        'examples': [
        ],
    },
], snakemake.output.dir)  
