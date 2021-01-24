import sys
import pandas as pd

# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
both_df = DfLoad.eng_both()

def count_bigprops(df, threshold=0.6):
    """Count proportions over a threshold in a prop table."""

    # count how many are attested
    ct_df = df > 0
    ct_df = 1 * ct_df
    ct_df = pd.DataFrame(ct_df.sum(), columns=['sum'])
    
    # count how many exceed threshold
    thr_df = df >= threshold
    thr_df = 1 * thr_df # convert to 0s and 1s
    thr_df = pd.DataFrame(thr_df.sum(), columns=['sum'])

    thr_pr = thr_df.div(ct_df)

    data = {
        'thresh_att_ct': ct_df, 
        f'thresh_ct>{threshold}': thr_df, 
        f'thresh_pr>{threshold}': thr_pr
    }
    
    return data
    

run_analyses([
    {
        'name': 'verb_lex',
        'df': both_df,
        'index': 'eng_TAMsimp',
        'columns': ['lex', 'stem'],
        'examples': [
            {
                'query': ('eng_TAMsimp == "FUT ~ MOD shall" '
                            'and lex_etcbc == "HJH[" '),
            },
            {
                'query': ('eng_TAMsimp == "FUT" '
                            'and lex_etcbc == "NTN[" '),
                'spread': 10,
            },
            {
                'query': ('eng_TAMsimp == "FUT ~ MOD shall" '
                            'and lex_etcbc == "MWT[" '),
                'spread': 10,
            },
            {
                'query': ('eng_TAMsimp == "FUT ~ MOD shall" '
                            'and lex_etcbc == "NPL[" '),
                'spread': 10,
            },
            {
                'query': ('eng_TAMsimp == "IMPV ~ MOD shall" '
                            'and lex_etcbc == "<FH[" '),
                'spread': 10,
            },
            {
                'query': ('eng_TAMsimp == "MOD may" '
                            'and lex_etcbc == ">KL[" '),
                'spread': 10,
            },
            {
                'query': ('eng_TAMsimp == "IMPV" '
                            'and lex_etcbc == "JR>[" '),
                'spread': 10,
            },
       ],
        'special': [
        ],
    },
], snakemake.output.dir)  
