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
        'index': 'esv_TAMsimp',
    },
    {
        'name': 'niv_tenses',
        'df': niv_df,
        'index': 'niv_TAMsimp',
    },
    {
        'name': 'trans_tam',
        'df': eng_df,
        'index': 'esv_TAMsimp',
        'columns': 'niv_TAMsimp',
        'fishers': False,
        'examples': [
        ],
    },
    {
        'name': 'disag_genre',
        'df': disag_df_simp,
        'index': 'eng_TAMsimp',
        'columns': 'genre',
    },
    {
        'name': 'disag_domain',
        'df': disag_df_simp,
        'index': 'eng_TAMsimp',
        'columns': 'domain2',
    },
    {
        'name': 'disag_gendom',
        'df': disag_df_simp,
        'index': 'eng_simp_agree',
        'columns': ['genre', 'domain2'],
    },
], snakemake.output.dir)  


