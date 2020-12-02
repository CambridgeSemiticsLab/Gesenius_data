from pathlib import Path
import pandas as pd

# path data
repo_dir = Path.home().joinpath('github/CambridgeSemiticsLab/Gesenius_data')
private_dir = repo_dir.joinpath('data/_private_')
plots_dir = repo_dir.joinpath('analysis/plots/qatal')
google_fid = private_dir.joinpath('keys/drive_folder.txt').read_text()

# dataframes specific to qatal
qatal_datapath = repo_dir.joinpath('data/_private_/verb_data/qatal_dataset.csv')
qatal_df = pd.read_csv(qatal_datapath, index_col='bhsa_node')
qatal_dfs = qatal_df[qatal_df.safe] # df data with filtered out parsing errors

# dataframes specific to all verbs
PUB_DIR = repo_dir.joinpath('data/_public_')
allverb_datapath = PUB_DIR.joinpath('verb_data/allverb_bhsa.csv')
av_df = pd.read_csv(allverb_datapath, index_col='bhsa_node')

# dataframes on verbal collocations
av_col_path = PUB_DIR.joinpath('verb_data/xverb_lexcollocations.csv')
av_col_df = pd.read_csv(av_col_path, index_col='verb_form')
