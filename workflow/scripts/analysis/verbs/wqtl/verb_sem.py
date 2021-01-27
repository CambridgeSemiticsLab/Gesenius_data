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
        'name': 'verb_stem',
        'df': both_df, 
        'index': 'eng_TAMsimp',
        'columns': 'stem',
    },
   {
        'name': 'verb_person',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'person',
        'examples': [
            {   
                'query': ('eng_TAMsimp == "FUT" '
                            'and person == "p1" ')
            },  
            {   
                'query': ('eng_TAMsimp == "FUT ~ MOD shall" '
                            'and person == "p1" ')
            },  
       ],
    },
    {
        'name': 'is_stative',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['esv_is', 'niv_is']
    },
], snakemake.output.dir)  
