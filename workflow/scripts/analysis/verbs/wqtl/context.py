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
            {   
                'query': ('eng_TAMsimp == "PRES" '
                            'and genre == "poetry" '),
            },
            {   
                'query': ('eng_TAMsimp == "IMPV ~ MOD shall" '
                            'and genre == "instruction" '),
                'spread': 15,
            },
            {   
                'query': ('eng_TAMsimp == "MOD is to ~ MOD shall" '
                            'and genre == "instruction" '),
                'spread': 15,
            },
            {   
                'query': ('eng_TAMsimp == "MOD shall" '
                            'and genre == "instruction" '),
                'spread': 15,
            },
            {   
                'query': ('eng_TAMsimp == "MOD must ~ MOD shall" '
                            'and genre == "instruction" '),
                'spread': 15,
            },
            {   
                'query': ('eng_TAMsimp == "PAST" '
                            'and genre == "poetry" '),
            },
            {   
                'query': ('eng_TAMsimp == "PAST" '
                            'and genre == "prose" '
                            'and domain2 == "N" '),
            },
 
        ],
    },  
    {   
        'name': 'domain',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'domain2',
        'examples': [
            {   
                'query': ('eng_TAMsimp == "PRES" '
                            'and domain2 == "D"'),
            },  
        ]   
    },  
    {   
        'name': 'gendom',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['genre', 'domain2'],
        'examples': [
            {   
                'query': ('eng_TAMsimp.isin(["FUT ~ MOD shall", "FUT"]) '
                            'and genre == "prophetic" '
                            'and domain2 == "Q"'),
            },  
            {   
                'query': ('eng_TAMsimp == "IMPV"'
                            'and genre == "prose" '
                            'and domain2 == "Q"'),
            },  
            { 
                'query': ('eng_TAMsimp == "IMPV"'
                            'and genre == "prose" '
                            'and domain2 == "Q"'),
            },  

       ],  
    },  
    {
        'name': 'ps_gendom',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['genre', 'domain2', 'person'],
    },
    { 
        'name': 'period',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'period',
    },  
    { 
        'name': 'period_gendom',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['period', 'genre', 'domain2'],
    },  
#    { 
#        'name': 'period_gendom_simp',
#        'df': both_df,
#        'index': 'eng_TAMsimp2',
#        'columns': ['period', 'genre', 'domain2'],
#    },  

], snakemake.output.dir)  


