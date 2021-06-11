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
        'name': 'has_objc',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': 'has_objc',
        'examples': [
            {
                'query': ('eng_TAMsimp == "PAST ~ PRES PART" '
                            'and has_objc == 1')
            },
            {
                'query': ('eng_TAMsimp == "PAST ~ PAST PROG" '
                            'and has_objc == 1')
            },
        ],
    },
    {
        'name': 'has_loca',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['has_loca'],
        'examples': [
        ],
    },
    {
        'name': 'has_time',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['has_time'],
        'examples': [
            {
                'query': ('eng_TAMsimp == "PAST PERF" '
                            ' and has_time == 1')
            },
        ],
    },
], snakemake.output.dir)  

