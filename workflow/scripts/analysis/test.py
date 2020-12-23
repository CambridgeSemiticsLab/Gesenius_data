from load_dfs import DfLoader

DfLoad = DfLoader('/Users/cody/github/CambridgeSemiticsLab/Gesenius_data/results/csv/qtl')

qatal_df = DfLoad.df_safe()
eng_df = DfLoad.eng_agree()

print(eng_df.head())
