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
        'name': 'has_objc',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'has_objc',
        'examples': [
       ],
    },
    {
        'name': 'has_loca',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['has_loca'],
        'examples': [
       ],
    },
    {
        'name': 'has_time',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['has_time'],
        'examples': [
       ],
    },
], snakemake.output.dir)  

