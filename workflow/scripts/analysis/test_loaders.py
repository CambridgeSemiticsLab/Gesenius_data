import pandas as pd
from load_dfs import DfLoader

DfLoad = DfLoader('/Users/cody/github/CambridgeSemiticsLab/Gesenius_data/results/csv/qtl')

qatal_df = DfLoad.df_safe()
eng_df = DfLoad.eng_agree()
esv_df = DfLoad.esv()
niv_df = DfLoad.niv()
disag_df = DfLoad.eng_disagree()

print(disag_df.shape)
print(disag_df.head())

#test = pd.pivot_table(
#    eng_df,
#    index='eng_TAM',
#    columns='has_objc',
#    aggfunc='size',
#    fill_value=0,
#    #dropna=False,
#)

print(test)
