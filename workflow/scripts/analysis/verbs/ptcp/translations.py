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

def sum_top_values(df):
    """Sums the top values of dataframes."""
    top = 0.015
    pr_df = df / df.sum()
    top_pr = pr_df.loc[pr_df['sum'] >= top]
    top_ct = df.loc[top_pr.index]
    top_sum = pd.DataFrame(top_ct.sum())
    sum_pr = pd.DataFrame(top_pr.sum())
    data = {
        'top_ct': top_ct, 
        'top_pr': top_pr, 
        'top_sum': top_sum,
        'top_sum_pr': sum_pr,
    }
    return data

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
            {'query': ('esv_TAM == "MOD shall"')}
        ]
    },
    {
        'name': 'niv_tenses',
        'df': niv_df,
        'index': 'niv_TAMsimp',
    },
    {
        'name': 'both_tenses',
        'df': eng_df,
        'index': 'eng_TAMsimp',
        'special': [
            {'df': 'count' , 'do': sum_top_values}
        ],
    },
    {
        'name': 'eng_simp_agree',
        'df': agg_df,
        'index': 'eng_simp_agree',
    },
    {
        'name': 'eng_simp_disagree',
        'df': disag_df_simp,
        'index': 'eng_TAMsimp',
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
        'name': 'both_genre',
        'df': eng_df,
        'index': 'eng_TAMsimp',
        'columns': 'genre',
    },
    {
        'name': 'both_domain',
        'df': eng_df,
        'index': 'eng_TAMsimp',
        'columns': 'domain2',
    },
    {
        'name': 'both_gendom',
        'df': eng_df,
        'index': 'eng_TAMsimp',
        'columns': ['genre', 'domain2'],
    },
    {
        'name': 'disag_genre',
        'df': disag_df_simp,
        'index': 'eng_TAMsimp',
        'columns': 'genre',
    },
    {
        'name': 'disag_domain',
        'df': disag_df_simp[disag_df_simp.domain2.isin(['N', 'Q'])],
        'index': 'eng_TAMsimp',
        'columns': 'domain2',
    },
    {
        'name': 'disag_gendom',
        'df': eng_df[eng_df.domain2.isin(['N', 'Q'])],
        'index': 'eng_simp_agree',
        'columns': ['genre', 'domain2'],
    },
    {
        'name': 'inter_gendom',
        'df': disag_df_simp,
        'index': 'eng_TAMsimp',
        'columns': ['genre', 'domain2'],
        'examples': [
        ],
    },

], snakemake.output.dir)  


