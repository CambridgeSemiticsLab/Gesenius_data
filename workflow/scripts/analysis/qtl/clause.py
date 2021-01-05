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
        'name': 'cltype_simp',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'cltype_simp',
    },
    {
        'name': 'rela_cltypesimp',
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
        'name': 'args',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'cl_args',
        'examples': [
            {
                'query': 'eng_TAM == "PAST..IND" and cl_args == "_W_SV"',
                'spread': 20,
            },
            {
                'query': 'eng_TAM == "PAST..IND" and cl_args == "_W_OV"',
                'spread': 20,
            },
            {
                'query': 'eng_TAM == "PAST..IND" and cl_args == "_W_AV"',
                'spread': 20,
            },
            {
                'query': 'eng_TAM == "PAST..IND" and cl_args == "SV"',
            },
        ]
    },
    {
        'name': 'args_mother',
        'df': eng_df,
        'index': ['eng_TAM', 'cl_args'],
        'columns': 'mother_verbtype',
        'examples': [
            {
                'query': ('eng_TAM == "PAST..IND" '
                            'and cl_args.str.match("_W_[OS]V") '
                            'and mother_verbtype == "wayq" '),
                'bhs_text': ['mother_clause', 'clause_atom']
            },
            {
                'query': ('eng_TAM == "PAST..IND" '
                            'and cl_args.str.match("_W_[OS]V") '
                            'and mother_verbtype == "wayq" '
                            'and mother_verb_lex == lex_etcbc '),
                'bhs_text': ['mother_clause', 'clause_atom'],
                'spread': 10,
            },

            {
                'query': ('eng_TAM == "PAST..IND"'
                            'and cl_args.str.match("_W_[OS]V") ' 
                            'and mother_verbtype == "wayq" '
                            'and mother_verb_ps == "p3" '
                            'and mother_verb_lex == "HJH[" '),
                'bhs_text': ['mother_clause', 'clause_atom'],
                'spread': 10,
            },

        ],
    },
    {
        'name': 'has_objc',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'has_objc',
    },
    {
        'name': 'has_loca',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': ['clause_rela', 'has_loca'],
        'examples': [
            {
                'query': ('eng_TAM == "PAST..IND" ' 
                            'and has_loca == 1 '
                            'and clause_rela == "Main"')
            }
        ],
    },
    {
        'name': 'has_time',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': ['clause_rela', 'has_time'],
        'examples': [
            {
                'query': ('eng_TAM == "PAST..IND" '
                            'and has_time  == 1 '
                            'and clause_rela == "Main" ')
            }
        ],
    },

], snakemake.output.dir)  

