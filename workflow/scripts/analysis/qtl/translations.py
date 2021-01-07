import sys
# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
eng_df = DfLoad.eng_agree()
esv_df = DfLoad.esv()
niv_df = DfLoad.niv()
both_df = DfLoad.eng_both()
disag_df = DfLoad.eng_disagree()

run_analyses([
    {
        'name': 'eng_tenses',
        'df': eng_df,
        'index': 'eng_TAM',
    },
    {
        'name': 'esv_tenses',
        'df': esv_df, 
        'index': 'esv_TAM',
    },
    {
        'name': 'niv_tenses',
        'df': niv_df,
        'index': 'niv_TAM',
    },
    {
        'name': 'eng_agree',
        'df': both_df,
        'index': 'eng_agree',
    },
    {
        'name': 'eng_disagree',
        'df': disag_df,
        'index': 'eng_TAM',
    },
    {
        'name': 'trans_tam',
        'df': both_df,
        'index': 'esv_TAM',
        'columns': 'niv_TAM',
        'fishers': False,
    },
    {
        'name': 'disag_past',
        'df': disag_df[disag_df.eng_TAM.str.match('.*PAST\.\.IND')],
        'index': 'eng_TAM',
    },
    {
        'name': 'disag_pres_perf',
        'df': disag_df[disag_df.eng_TAM.str.match('.*PRES\.PERF\.IND')],
        'index': 'eng_TAM',
        'examples': [
            {
                'query': ('eng_TAM == "FUT..IND ~ PRES.PERF.IND"'),
            },
            {
                'query': ('eng_TAM == "PRES..IND ~ PRES.PERF.IND"'),
            },
        ],
    },
    {
        'name': 'disag_domain',
        'df': disag_df[disag_df.domain2.isin(['N', 'Q'])],
        'index': 'eng_TAM',
        'columns': 'domain2',
    },

    {
        'name': 'disag_gendom',
        'df': both_df[both_df.domain2.isin(['N', 'Q'])],
        'index': 'eng_agree',
        'columns': ['genre', 'domain2'],
    },
    {
        'name': 'inter_gendom',
        'df': disag_df[disag_df.domain2.isin(['N', 'Q'])],
        'index': 'eng_TAM',
        'columns': ['genre', 'domain2'],
        'examples': [
            {
                'query': ('eng_TAM == "PAST..IND ~ PRES.PERF.IND" '
                            'and genre == "prose" '
                            'and domain2 == "Q"'),
                'spread': 10,
            },
            {
                'query': ('eng_TAM == "PAST..IND ~ PAST.PERF.IND" '
                            'and genre == "prose" '
                            'and domain2 == "N"'),

                'spread': 10,
            },
            {
                'query': ('eng_TAM == "PAST..IND ~ PRES.PERF.IND" '
                            'and genre == "prose" '
                            'and domain2 == "Q"'),

                'spread': 10,
            },
            {
                'query': ('eng_TAM == "PAST..IND ~ PRES..IND" '
                            'and genre.isin(["poetry", "prophetic"]) '
                            'and domain2 == "Q"'),

                'spread': 10,
            },
            {
                'query': ('eng_TAM == "PAST..IND ~ PRES..IND" '
                            'and genre == "instruction" '
                            'and domain2 == "Q"'),

                'spread': 2,
            },

        ],
    },

], snakemake.output.dir)  


