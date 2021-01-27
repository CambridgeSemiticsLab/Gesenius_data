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
        'name': 'is_question',
        'df': both_df,
        'index': 'is_progressive',
        'columns': ['is_question'],
        'examples': [
        ],
    },
], snakemake.output.dir)  

