import sys
import pandas as pd

# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
eng_df = DfLoad.eng_agree()

motion_verbs = ["BW>[", "HLK[", "JY>[", "JRD[", "QWM[", "NWS[", "SWR["]
punctuals = ["YWH[", "CLX[", "CB<[", "QR>[", 
            "NKH[", "JCB[", "JLD[", "KRT[", "FRP["]
statives = ["JD<[", ">HB[", "ML>[", "YRR[", "KLH[", "ZQN[", "XSH[", "RXQ[", "TMM[", "RYH["]

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
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': ['lex', 'stem'],
        'examples': [
            {
                'query': ('eng_TAM == "PAST..IND" '
                            f'and lex_etcbc.isin({motion_verbs}) '
                            'and stem == "qal"'),
            },
            {
                'query': f'eng_TAM == "PAST..IND" and lex_etcbc.isin({punctuals})',
            },
            {
                'query': 'eng_TAM == "PAST..IND" and lex_etcbc == "HJH["',
                'spread': 10,
            },
            {
                'query': 'eng_TAM == "PAST..IND" and lex_etcbc == "MLK["',
                'spread': 10,
            },
            {
                'query': 'eng_TAM == "PRES.PERF.IND" and lex_etcbc == "<FH[" and cl_args != "QV"',
                'spread': 5,
            },

            {
                'query': 'eng_TAM == "PRES.PERF.IND" and lex_etcbc == "NTN["',
                'spread': 5,
            },
            {
                'query': 'eng_TAM == "PRES.PERF.IND" and lex_etcbc == "FJM["',
                'spread': 5,
            },
            {
                'query': 'eng_TAM == "PRES.PERF.IND" and lex_etcbc == ">KL["',
                'spread': 5,
            },
            {
                'query': 'eng_TAM == "PRES.PERF.IND" and lex_etcbc == "CMR["',
                'spread': 5,
            },
            {
                'query': 'eng_TAM == "PRES.PERF.IND" and lex_etcbc == "BXR["',
                'spread': 5,
            },
            {
                'query': 'eng_TAM == "PRES.PERF.IND" and lex_etcbc == "CM<["',
                'spread': 5,
            },
            {
                'query': 'eng_TAM == "PRES.PERF.IND" and lex_etcbc == "R>H["',
                'spread': 5,
            },
            {
                'query': ('eng_TAM == "PRES..IND" '
                            f'and lex_etcbc.isin({statives}) '
                            'and stem == "qal"')
            },
            {
                'query': 'eng_TAM == "PAST.PERF.IND" and lex_etcbc == "<FH["',
            },
        ],
        'special': [
            {'df': 'prop2', 'do': count_bigprops},
        ],
    },
], snakemake.output.dir)  
