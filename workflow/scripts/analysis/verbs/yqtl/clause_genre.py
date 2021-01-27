import sys
# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from load_dfs import DfLoader
from analysis import run_analyses

# load the dataframes
DfLoad = DfLoader(snakemake.input.data_dir)
both_df = DfLoad.eng_both()


# features needed for selections
main_genre = ['prose', 'poetry', 'prophetic']
main_dom = ['Q', 'N']

run_analyses([
    {
        'name': 'clause_rela',
        'df': both_df,
        'index': ['eng_TAMsimp', 'genre', 'domain2'],
        'columns': 'clause_rela',
    },
   {   
        'name': 'prec_part',
        'df': both_df,
        'index': ['eng_TAMsimp', 'genre', 'domain2'],
        'columns': 'prec_part',
        'examples': [
            {   
                'query': ('eng_TAMsimp == "PRES" '
                            'and prec_part == "_>CR_" '
                            'and genre == "prose" ')
            },  
            {   
                'query': ('eng_TAMsimp == "PRES" '
                            'and prec_part == "_KJ_" '
                            'and genre == "prose" ')
            },  
            {   
                'query': ('eng_TAMsimp == "PRES" '
                            'and prec_part.isin(["_>M_", "_W_>M_"]) '
                            'and genre == "prose" ')
            },  
        ],  
    },  

], snakemake.output.dir)  

