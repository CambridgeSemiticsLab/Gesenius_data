from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from stats import significance as my_stats

repo_dir = Path.home().joinpath('github/CambridgeSemiticsLab/Gesenius_data')
private_dir = repo_dir.joinpath('data/_private_')
plots_dir = repo_dir.joinpath('analysis/plots/qatal')

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

def pivot_prop(data, index, columns, sum_i, div_i, **pivot_kwargs):
    """Construct pivot and prop tables at once"""
    

qatal_datapath = repo_dir.joinpath('data/_private_/verb_data/qatal_dataset.csv')
qatal_df = pd.read_csv(qatal_datapath, index_col='bhsa_node')

allverb_datapath = repo_dir.joinpath('data/_public_/verb_data/allverb_bhsa.csv')
av_df = pd.read_csv(allverb_datapath, index_col='bhsa_node')
