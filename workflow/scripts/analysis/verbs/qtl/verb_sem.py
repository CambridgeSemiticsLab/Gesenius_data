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
        'examples': [
            {
                'query': ('eng_TAM == "PRES..IND" '
                            'and stem == "nif" ')
            },
            {
                'query': ('eng_TAM == "PRES..IND" '
                            'and stem == "qal" ')
            },

        ],
    },
    {
        'name': 'verb_lexst_ps',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': ['lex', 'stem', 'person'],
        'fishers': False,
    },
    {
        'name': 'verb_person',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': 'person',
        'examples': [
            {
                'query': ('eng_TAM == "PRES..IND" '
                            'and person == "p1" '),
                'spread': 10,
            },
        ],
    },
    {
        'name': 'is_stative',
        'df': eng_df,
        'index': 'eng_TAM',
        'columns': ['esv_is', 'niv_is']
    },
], snakemake.output.dir)  
