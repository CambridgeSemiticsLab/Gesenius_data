import sys
# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
both_df = DfLoad.eng_both()


# features needed for selections
main_genre = ['prose', 'poetry', 'prophetic']
main_dom = ['Q', 'N']

run_analyses([
   {
        'name': 'args',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'cl_args',
        'examples': [
        ]
    },
], snakemake.output.dir)  

