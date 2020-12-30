import sys
# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
eng_df = DfLoad.eng_agree()

# features needed for selections
main_genre = ['prose', 'poetry', 'prophetic']
main_dom = ['Q', 'N']

run_analyses([
   {
        'name': 'clause_type',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'clause_type',
    },
    {
        'name': 'clause_rela',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'clause_rela',
    },
    {
        'name': 'clause_rela',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'clause_rela',
    },
    {
        'name': 'rela_x',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': ['clause_rela', 'cltype_simp'],
    },
    {
        'name': 'rela_particle',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': ['clause_rela', 'prec_part'],
    },
    {
        'name': 'rela_args',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': ['clause_rela', 'cl_args'],
    },
], snakemake.output.dir)  


