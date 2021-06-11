import sys

# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
both_df = DfLoad.eng_both()

run_analyses([
    {   
        'name': 'genre',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'genre',
        'examples': [
        ],
    },  
    {   
        'name': 'domain',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'domain2',
    },  
    {   
        'name': 'gendom',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['genre', 'domain2'],
        'examples': [
            {
                'query': (
                           'eng_TAMsimp == "PAST" '
                            'and genre == "prose" '
                            'and domain2 == "N" '),
                'spread': 10,
            },
            {
                'query': (
                           'eng_TAMsimp == "PAST" '
                            'and genre == "prose" '
                            'and domain2 == "Q" '),
                'spread': 10,
            },
            {
                'query': (
                           'eng_TAMsimp == "PAST" '
                            'and genre == "prophetic" '
                            'and domain2 == "N" '),
                'spread': 10,
            },
            {
                'query': (
                           'eng_TAMsimp == "PRES PERF" '
                            'and genre == "prose" '
                            'and domain2 == "Q" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "PRES PERF" '
                            'and genre == "prophetic" '
                            'and domain2 == "Q" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "PRES PERF" '
                            'and genre == "poetry" '
                            'and domain2 == "Q" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "PRES PERF" '
                            'and genre.isin(["prose", "prophetic"]) '
                            'and domain2 == "N" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "PAST PERF"'
                            'and genre == "list" '
                            'and domain2 == "N" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "PAST PERF"'
                            'and genre == "prose" '
                            'and domain2 == "N" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "PRES"'
                            'and genre == "poetry" '
                            'and domain2 == "Q" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "PRES"'
                            'and genre == "prophetic" '
                            'and domain2 == "Q" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "PRES"'
                            'and genre == "poetry" '
                            'and domain2 == "N" '),
            },
            {
                'query': (
                           'eng_TAMsimp.isin(["PRES PART", "PAST ~ PRES PART"])'
                            'and genre == "prose" '
                            'and domain2 == "Q" '),
            },
            {
                'query': (
                           'eng_TAMsimp.isin(["PRES PART", "PAST ~ PRES PART"])'
                            'and genre == "prose" '
                            'and domain2 == "N" '),
            },
            {
                'query': (
                           'eng_TAMsimp.isin(["TO INF", "PAST ~ TO INF"])'
                            'and genre == "prose" '
                            'and domain2 == "N" '),
            },
            {
                'query': (
                           'eng_TAMsimp.isin(["TO INF", "PAST ~ TO INF"])'
                            'and genre == "prose" '
                            'and domain2 == "Q" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "PAST ~ PAST PROG"'
                            'and genre == "prose" '
                            'and domain2 == "N" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "PAST ~ PAST PROG"'
                            'and genre == "prophetic" '
                            'and domain2 == "N" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "MOD would"'
                            'and genre == "prose" '
                            'and domain2 == "N" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "MOD would"'
                            'and genre == "prophetic" '
                            'and domain2 == "Q" '),
            },
            {
                'query': (
                           'eng_TAMsimp == "MOD would"'
                            'and genre == "poetry" '
                            'and domain2 == "Q" '),
            },
        ],  
    },  
    { 
        'name': 'period',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'period',
    },  
], snakemake.output.dir)  


