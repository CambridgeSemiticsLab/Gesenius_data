import sys

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

run_analyses([
    {
        'name': 'verb_stem',
        'df': eng_df, 
        'index': 'eng_TAM',
        'columns': 'stem',
    },
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
        ],
    },
    {
        'name': 'verb_person',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'person',
    },

], snakemake.output.dir)  
