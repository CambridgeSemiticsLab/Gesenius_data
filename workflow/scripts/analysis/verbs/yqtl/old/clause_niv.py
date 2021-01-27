import sys
# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
niv_df = DfLoad.niv()


# features needed for selections
main_genre = ['prose', 'poetry', 'prophetic']
main_dom = ['Q', 'N']

run_analyses([
   {
        'name': 'clause_type',
        'df': niv_df,
        'index': 'niv_TAM',
        'columns': 'clause_type',
    },
    {
        'name': 'clause_rela',
        'df': niv_df,
        'index': 'niv_TAM',
        'columns': 'clause_rela',
    },
    {
        'name': 'cltype_simp',
        'df': niv_df,
        'index': 'niv_TAM',
        'columns': 'cltype_simp',
    },
    {
        'name': 'rela_cltypesimp',
        'df': niv_df,
        'index': 'niv_TAM',
        'columns': ['clause_rela', 'cltype_simp'],
    },
   {
        'name': 'args',
        'df': niv_df,
        'index': 'niv_TAM',
        'columns': 'cl_args',
        'examples': []
    },
    {
        'name': 'rela_particle',
        'df': niv_df,
        'index': 'niv_TAM',
        'columns': ['clause_rela', 'prec_part'],
    },
    {
        'name': 'prec_part',
        'df': niv_df,
        'index': 'niv_TAM',
        'columns': 'prec_part',
        'examples': [
        ],
    },

#    {
#        'name': 'prec_part_gendom',
#        'df': niv_df,
#        'index': 'prec_part',
#        'columns': ['genre', 'domain2', 'eng_TAM'],
#    },
# 
#    {
#        'name': 'args_mother',
#        'df': niv_df,
#        'index': 'eng_TAM',
#        'columns': ['cl_args','mother_verbtype'],
#        'examples': [],
#    },
#    {
#        'name': 'has_objc',
#        'df': niv_df,
#        'index': 'eng_TAM',
#        'columns': 'has_objc',
#    },
#    {
#        'name': 'has_loca',
#        'df': niv_df,
#        'index': 'eng_TAM',
#        'columns': ['clause_rela', 'has_loca'],
#        'examples': [
#        ],
#    },
#    {
#        'name': 'has_time',
#        'df': niv_df,
#        'index': 'eng_TAM',
#        'columns': ['clause_rela', 'has_time'],
#        'examples': [
#        ],
#    },

], snakemake.output.dir)  

