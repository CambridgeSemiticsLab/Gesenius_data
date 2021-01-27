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
        'index': 'eng_TAMsimp',
        'columns': 'clause_type',
    },
    {
        'name': 'main_clause_type',
        'df': both_df[both_df.clause_rela == 'Main'],
        'index': 'eng_TAMsimp',
        'columns': 'clause_type',
    },
    {
        'name': 'clause_rela',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'clause_rela',
    },
    {
        'name': 'cltype_simp',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'cltype_simp',
    },
    {
        'name': 'rela_cltypesimp',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['clause_rela', 'cltype_simp'],
    },
   {
        'name': 'args',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'cl_args',
        'examples': [
            {
                'query': ('eng_TAMsimp == "IMPV ~ MOD shall" '
                            'and cl_args.isin(["OV", "_W_OV"])'),
                'spread': 15,
            },
            {
                'query': ('eng_TAMsimp == "MOD must ~ MOD shall" '
                            'and cl_args.isin(["OV", "_W_OV"])'),
                'spread': 15,
            },
            {
                'query': ('eng_TAMsimp == "MOD can" '
                            'and cl_args == "QV" '),
                'spread': 15,
            }
        ]
    },
   {   
        'name': 'prec_part',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'prec_part',
        'examples': [
            {
                'query': ('eng_TAMsimp == "IMPV" '
                            'and prec_part == "_>L=_" ')
            },
            {
                'query': ('eng_TAMsimp == "MOD let" '
                            'and prec_part == "Ø" ')
            },
            {
                'query': ('eng_TAMsimp == "MOD may" '
                            'and prec_part == "_W_" ')
            },
        ]
    },
    {
        'name': 'mo_verbtype',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'mother_verbtype',
        'examples': [
            {
                'query': ('eng_TAMsimp == "PAST" '
                            'and mother_verbtype == "wayq"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PAST" '
                            'and mother_verbtype == "qtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
        ]
    },

], snakemake.output.dir)  

