"""
Analyze qatal contexts.
"""

import sys
# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analysis

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
qatal_df = DfLoad.df_safe()
eng_df = DfLoad.eng_agree()

# features needed for selections
main_genre = ['prose', 'poetry', 'prophetic']
main_dom = ['Q', 'N']

analysis_params = [
    {
        'name': 'genre',
        'df': eng_df[eng_df.genre.isin(main_genre)],
        'index': 'eng_TAM',
        'columns': 'genre',
    },
    {
        'name': 'domain',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'domain2', 
    },
    {
        'name': 'period',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'period',
    },
    {
        'name': 'period_genre',
        'df': eng_df[eng_df.genre.isin(main_genre)],
        'index': 'eng_TAM',
        'columns': ['period', 'genre'],
    },
    {
        'name': 'period_domain',
        'df': eng_df[eng_df.domain2.isin(main_dom)],
        'index': 'eng_TAM',
        'columns': ['period', 'domain2'],
    },
    {
        'name': 'genre_domain',
        'df': eng_df[(eng_df.domain2.isin(main_dom)) & (eng_df.genre.isin(main_genre))],
        'index': ['genre', 'domain2'],
        'columns': 'eng_TAM',
    },
]  

for analysis in analysis_params:
    run_analysis(analysis, snakemake.output[0])

#     [
#         'vstem',
#         [engag_df],
#         {'index': ['esv_TAM'], 'columns':['stem']}
#     ],
#     [
#         'vlex',
#         [engag_df],
#         {'index': ['esv_TAM'], 'columns':['lex', 'stem']}
#     ],
#     [
#         'vpers',
#         [engag_df],
#         {'index': ['esv_TAM'], 'columns': ['person']}
#     ],
