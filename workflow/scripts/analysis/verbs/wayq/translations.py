import sys
import pandas as pd

# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
agg_df = DfLoad.eng_simp_agree()
esv_df = DfLoad.esv()
niv_df = DfLoad.niv()
eng_df = DfLoad.eng_both()
disag_df_simp = DfLoad.eng_simp_disagree()
disag_df = DfLoad.eng_disagree()

run_analyses([
    {
        'name': 'eng_tenses',
        'df': eng_df,
        'index': 'eng_TAMsimp',
        'examples': [
        ],
    },
    {
        'name': 'esv_tenses',
        'df': esv_df, 
        'index': 'esv_TAM',
        'examples': [
        ]
    },
    {
        'name': 'niv_tenses',
        'df': niv_df,
        'index': 'niv_TAMsimp',
        'examples': [
            {'query': 'niv_TAMsimp == "IMPV"'}
        ]
    },
    {
        'name': 'eng_simp_disagree',
        'df': disag_df_simp,
        'index': 'eng_TAMsimp',
        'examples': [
            {
                'query': 'eng_TAMsimp == "PAST ~ PRES PART"',
                'spread': 10,
            },
            {
                'query': 'eng_TAMsimp == "PAST ~ PAST PERF"',
                'spread': 10,
            },
            {
                'query': 'eng_TAMsimp == "PAST ~ PRES PERF"',
                'spread': 10,
            },
           {
                'query': 'eng_TAMsimp == "PAST ~ PRES"',
                'spread': 10,
            },
           {
                'query': 'eng_TAMsimp == "PAST ~ TO INF"',
                'spread': 10,
            },
           {
                'query': 'eng_TAMsimp == "PRES ~ PRES PERF"',
                'spread': 10,
            },
           {
                'query': 'eng_TAMsimp == "PAST ~ PAST PROG"',
                'spread': 10,
            },
           {
                'query': 'eng_TAMsimp == "FUT ~ PRES"',
                'spread': 10,
            },
           {
                'query': 'eng_TAMsimp == "MOD could ~ PAST"',
                'spread': 10,
            },
           {
                'query': 'eng_TAMsimp == "FUT ~ PAST"',
                'spread': 10,
            },
           {
                'query': 'eng_TAMsimp.str.match(".*FUT")',
                'spread': -1,
            },
        ]
    },
], snakemake.output.dir)  


