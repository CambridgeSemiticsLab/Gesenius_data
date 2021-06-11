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
   {
        'name': 'mo_verbtype',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'mother_verbtype',
        'examples': [
            {
                'query': ('eng_TAMsimp == "PAST"'
                            'and mother_verbtype == "wayq"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PAST"'
                            'and mother_verbtype == "qtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PRES PERF"'
                            'and mother_verbtype == "qtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp.isin(["PAST PERF", "PAST ~ PAST PERF"])'
                            'and mother_verbtype == "qtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp.isin(["PAST PERF", "PAST ~ PAST PERF"])'
                            'and mother_verbtype == "wayq"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PRES"'
                            'and mother_verbtype == "yqtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PRES"'
                            'and mother_verbtype == "qtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PRES"'
                            'and mother_verbtype == "ptcp"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PRES"'
                            'and mother_verbtype == "wayq"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PAST PROG"'
                            'and mother_verbtype == "qtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PAST PROG"'
                            'and mother_verbtype == "ptcp"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PAST PROG"'
                            'and mother_verbtype == "wqtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "FUT ~ PRES"'
                            'and mother_verbtype == "qtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "FUT ~ PAST"'
                            'and mother_verbtype == "wayq"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
              {
                'query': ('eng_TAMsimp == "PAST ~ PRES PART" '
                            'and mother_verbtype == "infa"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },
            {
                'query': ('eng_TAMsimp == "PRES PART" '
                            'and mother_verbtype == "qtl"'),
                'bhs_text': ['mother_clause_atom', 'clause_atom'],
            },


        ]
    },
   {
        'name': 'verb_lex',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'lex',
        'examples': [
            {
                'query': ('eng_TAMsimp == "PAST" '
                            'and lex == "אמר"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST" '
                            'and lex == "בוא"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST" '
                            'and lex == "עשׂה"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST" '
                            'and lex == "עלה"'),
            },
            {
                'query': ('eng_TAMsimp == "PRES PERF" '
                            'and lex == "מלט"'),
            },
            {
                'query': ('eng_TAMsimp == "PRES PERF" '
                            'and lex == "מאס"'),
            },
            {
                'query': ('eng_TAMsimp == "PRES PERF" '
                            'and lex == "היה"'),
            },
            {
                'query': ('eng_TAMsimp == "PRES PERF" '
                            'and lex == "נתן"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PAST PERF" '
                            'and lex == "לקח"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PAST PERF" '
                            'and lex == "נתן"'),
            },
            {
                'query': ('eng_TAMsimp == "PRES" '
                            'and lex.isin(["ישׁב", "ירא", "ידע", "מלא"])'),
                'spread': 35,
            },
            {
                'query': ('eng_TAMsimp == "PRES" '
                            'and lex == "אמר"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PRES PART" '
                            'and lex == "חוה"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PRES PART" '
                            'and lex == "זעק"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PRES PART" '
                            'and lex == "שׁבר"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PRES PART" '
                            'and lex == "יעץ"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PAST PROG" '
                            'and lex == "הלך"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PAST PROG" '
                            'and lex == "ישׁב"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PAST PROG" '
                            'and lex == "אכל"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PAST PROG" '
                            'and lex == "היה"'),
            },

            {
                'query': ('eng_TAMsimp == "PAST ~ TO INF" '
                            'and lex == "חוה"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ TO INF" '
                            'and lex == "ישׁב"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ TO INF" '
                            'and lex == "שׁתה"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ TO INF" '
                            'and lex == "חטא"'),
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ TO INF" '
                            'and lex == "כלה"'),
            },
        ]
    },
], snakemake.output.dir)  

