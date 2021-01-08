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
        'name': 'genre',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'genre',
    },
    {
        'name': 'domain',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'domain2', 
        'examples': [
            {
                'query': ('eng_TAM == "PRES..IND" '
                            'and domain2.str.match("[ND]") '
                            'and genre == "prose"'),
            }
        ],
    },
    {
        'name': 'period',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'period',
    },
    {
        'name': 'period_gendom',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': ['period', 'genre', 'domain2'],
        'examples': [
            {
                'query': ('eng_TAM == "PRES..IND" '
                            'and period == "SBH"  '
                            'and genre == "prose" '
                            'and domain2 == "Q" '),
            },
            {
                'query': ('eng_TAM == "PRES..IND" '
                            'and period == "LBH"  '
                            'and genre == "prose" '
                            'and domain2 == "Q" '),
            },
        ],
    },
    {
        'name': 'period_domain',
        'df': eng_df[eng_df.domain2.isin(main_dom)],
        'index': 'eng_TAM',
        'columns': ['period', 'domain2'],
    },
    {
        'name': 'genre_domain',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': ['genre', 'domain2'], 
        'examples': [
            {
                'query': ('eng_TAM == "PRES.PERF.IND" '
                            'and domain2 == "N" '
                            'and (~txt_type.str.match("[QD]"))')
            },
            {
                'query': ('eng_TAM == "PRES.PERF.IND" '
                            'and domain2.str.match("[ND?]") '
                            'and genre == "prose"'),
            },
            {
                'query': ('eng_TAM == "PRES.PERF.IND" '
                            'and genre == "instruction" '
                            'and domain2 == "Q" '),
                'spread': 15,
            },
            {
                'query': ('eng_TAM == "PAST.PERF.IND" '
                            'and genre == "prose" '
                            'and domain2 == "N" '),
            },
        ],
    },

], snakemake.output.dir)  


