"""
Analyze qatal contexts.
"""

import os
import sys
import json
from pathlib import Path

# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
qatal_df = DfLoad.df_safe()
eng_df = DfLoad.eng_agree()

# features needed for selections
main_genre = ['prose', 'poetry', 'prophetic']
main_dom = ['Q', 'N']

run_analyses([

    # broader contexts
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

    # verb categories
    {
        'name': 'verb_stem',
        'df': eng_df, 
        'index': 'eng_TAM',
        'columns': 'stem',
    },
    {
        'name': 'verb_lex',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': ['lex', 'stem'],
    },
    {
        'name': 'verb_person',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'person',
    },

], snakemake.output.dir)  


