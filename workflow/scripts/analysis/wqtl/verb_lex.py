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
            {
                'query': ('eng_TAMsimp == "FUT" '
                            'and lex_etcbc == "NTN["'),
            },
            {
                'query': ('eng_TAMsimp == "FUT ~ MOD shall" '
                            'and lex_etcbc == "HJH["'),
            },
            {
                'query': ('eng_TAMsimp == "FUT ~ MOD shall" '
                            'and lex_etcbc == "JD<["'),
                'spread': 10,
            },
            {
                'query': ('eng_TAMsimp == "IMPV" '
                            'and lex_etcbc == ">MR["'),
                'spread': 10,
            },
            {
                'query': ('eng_TAMsimp == "IMPV ~ MOD shall" '
                            'and lex_etcbc == "<FH["'),
                'spread': 10,
            },

        ],
    },
], snakemake.output.dir)  
