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
                'query': ('genre == "poetry" '
                          'and eng_TAMsimp == "PRES"'),
                'spread': 15,
            },
            {
                'query': ('genre == "prophetic" '
                          'and eng_TAMsimp == "PRES"'),
                'spread': 15,
            },
            {
                'query': ('genre == "prose" '
                          'and domain2 == "Q" '
                          'and eng_TAMsimp == "PRES"'),
                'spread': 15,
            },
            {
                'query': ('genre == "prose" '
                          'and domain2 == "N" '
                          'and eng_TAMsimp == "PAST"'),
                'spread': 15,
            },
            {
                'query': ('genre == "prose" '
                          'and domain2 == "Q" '
                          'and eng_TAMsimp == "PAST"'),
                'spread': 15,
            },
             {
                'query': ('genre == "prophetic" '
                          'and eng_TAMsimp == "PAST"'),
                'spread': 15,
            },
             {
                'query': ('genre == "poetry" '
                          'and eng_TAMsimp == "PAST"'),
                'spread': 15,
            },
            {
                'query': ('genre == "prose" '
                          'and domain2 == "N" '
                          'and eng_TAMsimp == "PRES PART"'),
                'spread': 15,
            },
            {
                'query': ('genre == "prose" '
                          'and domain2 == "Q" '
                          'and eng_TAMsimp == "PRES PART"'),
                'spread': 15,
            },
            {
                'query': ('genre == "poetry" '
                          'and eng_TAMsimp == "PRES PART"'),
                'spread': 15,
            },
            {
                'query': ('genre == "prophetic" '
                          'and domain2 == "Q" '
                          'and eng_TAMsimp == "PRES PART"'),
                'spread': 15,
            },
            {
                'query': ('genre == "prose" '
                          'and domain2 == "Q" '
                          'and eng_TAMsimp == "PRES PROG"'),
                'spread': 15,
            },
            {
                'query': ('genre == "instruction" '
                          'and domain2 == "Q" '
                          'and eng_TAMsimp == "PRES PROG"'),
                'spread': 15,
            },
            {
                'query': ('genre == "prophetic" '
                          'and domain2 == "Q" '
                          'and eng_TAMsimp == "PRES PROG"'),
                'spread': 15,
            },
 
            {
                'query': ('genre == "prose" '
                          'and domain2 == "N" '
                          'and eng_TAMsimp == "PAST PROG"'),
                'spread': 15,
            },
            {
                'query': ('genre == "prophetic" '
                          'and eng_TAMsimp == "FUT"'),
                'spread': 15,
            },
            {
                'query': ('genre == "prose" '
                          'and domain2 == "Q" '
                          'and eng_TAMsimp == "FUT"'),
                'spread': 15,
            },
            {
                'query': ('genre == "poetry" '
                          'and domain2 == "Q" '
                          'and eng_TAMsimp == "FUT"'),
                'spread': 15,
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


