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
        'name': 'clause_type',
        'df': both_df,
        'index': ['eng_TAMsimp', 'person'],
        'columns': 'clause_type',
        'examples': [
            {
                'query': ('eng_TAMsimp.isin(["FUT", "FUT ~ MOD shall"]) '
                            'and person == "p3" '
                            'and clause_type.isin(["xYqX"])')
            },
        ],
    },
    {
        'name': 'cltype_maincl',
        'df': both_df[both_df.clause_rela == 'Main'],
        'index': ['eng_TAMsimp', 'person'],
        'columns': 'clause_type',
        'examples': [
        ]
    },

   {
        'name': 'args',
        'df': both_df,
        'index': ['eng_TAMsimp', 'person'],
        'columns': 'cl_args',
        'examples': [
            {
                'query': ('eng_TAMsimp.isin(["FUT", "FUT ~ MOD shall"]) '
                            'and person == "p3" '
                            'and cl_args.isin(["_W_SV", "SV"]) '),
            },

        ]
    },
    {
        'name': 'args_maincl',
        'df': both_df[both_df.clause_rela == 'Main'],
        'index': ['eng_TAMsimp', 'person'],
        'columns': 'cl_args',
        'examples': [
        ]
    },

#
#    {
#        'name': 'main_clause_type',
#        'df': both_df[both_df.clause_rela == 'Main'],
#        'index': 'eng_TAMsimp',
#        'columns': 'clause_type',
#    },
#    {
#        'name': 'clause_rela',
#        'df': both_df,
#        'index': 'eng_TAMsimp',
#        'columns': 'clause_rela',
#    },
#    {
#        'name': 'cltype_simp',
#        'df': both_df,
#        'index': 'eng_TAMsimp',
#        'columns': 'cltype_simp',
#    },
#    {
#        'name': 'rela_cltypesimp',
#        'df': both_df,
#        'index': 'eng_TAMsimp',
#        'columns': ['clause_rela', 'cltype_simp'],
#    },
#    {
#        'name': 'prec_part',
#        'df': both_df,
#        'index': 'eng_TAMsimp',
#        'columns': 'prec_part',
#        'examples': [
#        ],
#    },

], snakemake.output.dir)  

