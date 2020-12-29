from load_dfs import DfLoader

DfLoad = DfLoader('/Users/cody/github/CambridgeSemiticsLab/Gesenius_data/results/csv/qtl')

qatal_df = DfLoad.df_safe()
eng_df = DfLoad.eng_agree()
esv_df = DfLoad.esv()
niv_df = DfLoad.niv()
