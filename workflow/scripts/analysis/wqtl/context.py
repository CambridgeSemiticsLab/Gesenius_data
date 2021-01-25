import sys

# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
both_df = DfLoad.eng_both()

run_analyses([
    {   
        'name': 'genre',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'genre',
        'examples': [
       ],
    },  
    {   
        'name': 'domain',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'domain2',
    },  
    {   
        'name': 'gendom',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['genre', 'domain2'],
        'examples': [
       ],  
    },  
    {
        'name': 'ps_gendom',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['genre', 'domain2', 'person'],
    },
    { 
        'name': 'period',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'period',
    },  
     { 
        'name': 'period_gendom',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['period', 'genre', 'domain2'],
    },  
], snakemake.output.dir)  


