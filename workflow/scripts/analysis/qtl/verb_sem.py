import sys

# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
eng_df = DfLoad.eng_agree()

run_analyses([
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
