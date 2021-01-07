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
        'df': eng_df[eng_df.genre.isin(main_genre)],
        'index': 'eng_TAM',
        'columns': 'genre',
    },
    {
        'name': 'domain',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'domain2', 
    },
    {
        'name': 'period',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'period',
    },
    {
        'name': 'period_gendom',
        'df': eng_df[eng_df.domain2.isin(main_dom)],
        'index': 'eng_TAM',
        'columns': ['period', 'genre', 'domain2'],
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
        'index': ['genre', 'domain2'],
        'columns': 'eng_TAM',
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

        ],
    },

], snakemake.output.dir)  


