import collections
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from stats import significance as my_stats

repo_dir = Path.home().joinpath('github/CambridgeSemiticsLab/Gesenius_data')
private_dir = repo_dir.joinpath('data/_private_')
plots_dir = repo_dir.joinpath('analysis/plots/qatal')

qatal_datapath = repo_dir.joinpath('data/_private_/verb_data/qatal_dataset.csv')
qatal_df = pd.read_csv(qatal_datapath, index_col='bhsa_node')
qatal_dfs = qatal_df[qatal_df.safe] # df data with filtered out parsing errors

PUB_DIR = repo_dir.joinpath('data/_public_')
allverb_datapath = PUB_DIR.joinpath('verb_data/allverb_bhsa.csv')
av_df = pd.read_csv(allverb_datapath, index_col='bhsa_node')

av_col_path = PUB_DIR.joinpath('verb_data/xverb_lexcollocations.csv')
av_col_df = pd.read_csv(av_col_path, index_col='verb_form')

def get_props(counts_df, sum_i=1, div_i=0):
    return counts_df.div(counts_df.sum(sum_i), div_i) 

def save_fig(filename):
    filepath = plots_dir.joinpath(filename+'.svg')
    plt.savefig(filepath, format='svg', bbox_inches='tight')
    
def plot_bar_1D(data_1D, ax, title='', xlabel='', ylabel='', **plot_kwargs):
    """Make barplot from 1D data to show simple counts."""
    if not {'color', 'cmap', 'edgecolor'} & set(plot_kwargs):
        plot_kwargs.update({'color': 'lightgrey', 'edgecolor':'black'})
    data_1D.plot(kind='bar', ax=ax, **plot_kwargs)
    ax.set_title(title)
    ax.grid(axis='y')
    ax.set_axisbelow(True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

def pivot_table(data, index, columns, **pivot_kwargs):
    """Make a pivot table with standard options"""
    pivot_options = {
        'aggfunc': 'size',
        'fill_value': 0,
    }
    pivot_options.update(pivot_kwargs)
    return pd.pivot_table(
        data,
        index=index,
        columns=columns,
        **pivot_options
    )

def prop_table(count_data, sumi=1, divi=0):
    """Make a prop table from a table of counts"""
    if len(count_data.shape) == 1:
        return count_data / count_data.sum()
    else:
        return count_data.div(count_data.sum(sumi), divi)

class PivotProp:
    """Construct counts and proportion tables and store as class attribs"""
    def __init__(self, data, index, columns, sum_i=1, div_i=0, **pivot_kwargs):

        # calculuate count table
        self.ct = pivot_table(data, index, columns, **pivot_kwargs)
        # sort the ct based on biggest sums
        sort_sums = self.ct.sum().sort_values(ascending=False).index
        self.ct = self.ct[sort_sums]
        
        # calculate prop table
        self.pr = prop_table(self.ct, sumi=sum_i, divi=div_i)

def pretty_hebrew(val):
    """Render Hebrew in a dataframe."""
    return 'font-size:20px; font-family: Times New Roman; text-align: right; max-width: 500px'

def get_spread(array, n): 
    """Retrieve an even spread of indices an array/Series.
    
    https://stackoverflow.com/a/50685454/8351428
    
    Args:
        array: either numpy array or Pandas Series 
            (with multiple indices allowed)
        n: number of indices to return
    Returns:
        indexed array or series
    """
    end = len(array) - 1 
    spread = np.ceil(np.linspace(0, end, n)).astype(int)
    indices = np.unique(spread)
    try:
        return array[indices]
    except KeyError:
        return array.iloc[indices]

def show_text(df, col_default=['ref', 'sentence', 'text_full', 'lxx', 'lxx_tm', 'esv', 'esv_TAM'],
              cols=[], spread=0):
    """Display text from pandas dataframe in a readable way."""
    original_shape = df.shape
    df = df[col_default + cols]
    if spread > 0:
        spread_i = get_spread(df.index, spread)
        df = df.loc[spread_i]
    print(f'showing {df.shape[0]} of {original_shape[0]}')
    return df.style.applymap(pretty_hebrew, subset=['text_full', 'sentence'])
