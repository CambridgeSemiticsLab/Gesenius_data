import sys
# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
eng_df = DfLoad.eng_agree()
#esv_df = DfLoad.esv()
#niv_df = DfLoad.niv()
#both_df = DfLoad.eng_both()
#disag_df = DfLoad.eng_disagree()



run_analyses([
    {
        'name': 'inchoatives',
        'df': eng_df,
        'examples': [
            {
                'query': (
                    'eng_TAM == "PAST" '
                    'and (niv.str.match(".*became") | esv.str.match(".*became"))'
                ),
                'spread': -1, # i.e. all
            } 
        ]
    },
], snakemake.output.dir)  


