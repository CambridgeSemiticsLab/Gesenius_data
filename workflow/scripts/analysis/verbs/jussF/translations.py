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
            {
                'query': ('esv_TAM == "PAST"')
            },

        ]
    },
    {
        'name': 'niv_tenses',
        'df': niv_df,
        'index': 'niv_TAM',
        'examples': [
            {
                'query': ('niv_TAM == "PAST"')
            },
        ]
    },
    {
        'name': 'eng_simp_disagree',
        'df': disag_df_simp,
        'index': 'eng_TAMsimp',
        'examples': [
            {
                'query': ('esv_TAM == "MOD let" or niv_TAM == "MOD let"')
            },
            {
                'query': ('esv_TAM == "MOD may" or niv_TAM == "MOD may"')
            },
            {
                'query': ('esv_TAM == "IMPV" or niv_TAM == "IMPV"')
            },
            {
                'query': ('esv_TAM == "IMPV do not" or niv_TAM == "IMPV do not"')
            },
            {
                'query': ('esv_TAM == "MOD shall" or niv_TAM == "MOD shall"')
            },
            {
                'query': ('esv_TAM == "FUT" or niv_TAM == "FUT"')
            },
            {
                'query': ('esv_TAM == "PRES" or niv_TAM == "PRES"')
            },
            {
                'query': ('esv_TAM == "TO INF" or niv_TAM == "TO INF"')
            },
            {
                'query': ('esv_TAM == "MOD would" or niv_TAM == "MOD would"')
            },
            {
                'query': ('esv_TAM == "PAST" or niv_TAM == "PAST"')
            },
            {
                'query': ('esv_TAM == "MOD must" or niv_TAM == "MOD must"')
            },
            {
                'query': ('esv_TAM == "MOD should" or niv_TAM == "MOD should"')
            },
            {
                'query': ('esv_TAM == "PRES PART" or niv_TAM == "PRES PART"')
            },
            {
                'query': ('esv_TAM == "MOD is to" or niv_TAM == "MOD is to"')
            },
            {
                'query': ('esv_TAM == "MOD can" or niv_TAM == "MOD can"')
            },


        ]
    },
], snakemake.output.dir)  
