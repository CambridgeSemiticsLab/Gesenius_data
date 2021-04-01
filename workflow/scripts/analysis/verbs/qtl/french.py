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
        'name': 'french_tense',
        'df': both_df, 
        'index': 'eng_TAMsimp',
        'columns': 'french_tense',
        'examples': [
            {
                'query': ('eng_TAMsimp == "PAST" '
                          'and french_tense == "imparfait"'),
                'spread': -1,
                'extra_text': {'NBS': 'french'},
            },
            {
                'query': ('eng_TAMsimp == "PAST" '
                          'and french_tense == "passé_comp"'),
                'spread': -1,
                'extra_text': {'NBS': 'french'},
            },
            {
                'query': ('eng_TAMsimp == "PAST" '
                          'and french_tense == "passé_simp"'),
                'spread': -1,
                'extra_text': {'NBS': 'french'},
            },
            {
                'query': ('eng_TAMsimp == "PRES" '
                          'and french_tense == "passé_comp"'),
                'spread': -1,
                'extra_text': {'NBS': 'french'},
            },
        ],
    },

], snakemake.output.dir)  
