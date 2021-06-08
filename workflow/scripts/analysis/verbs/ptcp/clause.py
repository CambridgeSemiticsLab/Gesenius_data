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
            {
                'query': ('eng_TAMsimp == "PRES" '
                            'and cl_args == "V"')
            },
            {
                'query': ('eng_TAMsimp == "PRES" '
                            'and cl_args == "RV"')
            },
             {
                'query': ('eng_TAMsimp == "PRES" '
                            'and cl_args == "QSV"')
            },
             {
                'query': ('eng_TAMsimp == "PAST" '
                            'and cl_args == "RV"')
            },
              {
                'query': ('eng_TAMsimp == "PAST" '
                            'and cl_args == "_W_SV"')
            },
              {
                'query': ('eng_TAMsimp == "PAST" '
                            'and cl_args == "V"')
            },
            {
                'query': ('eng_TAMsimp == "PRES PART" '
                            'and cl_args == "V"')
            },
             {
                'query': ('eng_TAMsimp == "PRES PROG" '
                            'and cl_args == "RSV"')
            },
             {
                'query': ('eng_TAMsimp == "PRES PROG" '
                            'and cl_args == "ISV"')
            },
            {
                'query': ('eng_TAMsimp == "PAST PROG" '
                            'and cl_args == "_W_SV"')
            },
              {
                'query': ('eng_TAMsimp == "PAST PROG" '
                            'and cl_args == "SV"')
            },
            {
                'query': ('eng_TAMsimp == "PAST PROG" '
                            'and cl_args == "RV"')
            },
             {
                'query': ('eng_TAMsimp == "FUT" '
                            'and cl_args == "ISV"')
            },
             {
                'query': ('eng_TAMsimp == "FUT" '
                            'and cl_args == "CV"')
            },
             {
                'query': ('eng_TAMsimp == "FUT" '
                            'and cl_args == "ASV"')
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
                'query': ('eng_TAMsimp == "PRES"'
                            'and mother_verbtype == "yqtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PAST" '
                            'and mother_verbtype == "wayq"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PRES PART" '
                            'and mother_verbtype == "wayq"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PRES PROG" '
                            'and mother_verbtype == "infc"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
             {
                'query': ('eng_TAMsimp == "PRES PROG" '
                            'and mother_verbtype == "yqtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
             {
                'query': ('eng_TAMsimp == "PRES PROG" '
                            'and mother_verbtype == "wqtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
             {
                'query': ('eng_TAMsimp == "PRES PROG" '
                            'and mother_verbtype == "impv"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
             {
                'query': ('eng_TAMsimp == "PAST PROG" '
                            'and mother_verbtype == "wayq"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
 
             {
                'query': ('eng_TAMsimp == "FUT" '
                            'and mother_verbtype == "qtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
 
 
        ]
    },

], snakemake.output.dir)  

