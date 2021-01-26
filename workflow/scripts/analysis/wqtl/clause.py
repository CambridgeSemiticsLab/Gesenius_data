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
        'examples': [
            {
                'query': ('eng_TAMsimp == "FUT ~ MOD shall" '
                            'and clause_type == "WQtX"'),
            },
            {
                'query': ('eng_TAMsimp == "MOD shall" '
                            'and clause_type == "WQtX"'),
            },
        ],
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
        ]
    },
   {
        'name': 'mo_verbtype',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'mother_verbtype',
        'examples': [
            {
                'query': ('eng_TAMsimp.isin(["FUT", "FUT ~ MOD shall"]) '
                            'and mother_verbtype == "ptcp"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PRES" '
                            'and mother_verbtype == "yqtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "IMPV" '
                            'and mother_verbtype == "impv"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "MOD must ~ MOD shall" '
                            'and mother_verbtype == "Ø"'),
                'bhs_text': ['mother_intertext', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PAST" '
                            'and mother_verbtype == "wayq"'),
                'bhs_text': ['mother_intertext', 'clause_atom'],
                'spread': 35,
            },
            {
                'query': ('eng_TAMsimp == "PAST" '
                            'and mother_verbtype == "qtl"'),
                'bhs_text': ['mother_intertext', 'clause_atom'],
                'spread': 35,
            },
        ]
    },

], snakemake.output.dir)  

