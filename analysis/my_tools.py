import sys
import collections
import gspread # for pushing to Google drive
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from stats import significance as my_stats

# set pandas display settings
pd.set_option('display.max_columns', 0) # show unlimited columns

repo_dir = Path.home().joinpath('github/CambridgeSemiticsLab/Gesenius_data')
private_dir = repo_dir.joinpath('data/_private_')
plots_dir = repo_dir.joinpath('analysis/plots/qatal')
google_fid = private_dir.joinpath('keys/drive_folder.txt').read_text()

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
    
def plot_bar_1D(data_1D, ax=None, title='', xlabel='', ylabel='', **plot_kwargs):
    """Make barplot from 1D data to show simple counts."""
    if not ax:
        fig, ax = plt.subplots()
    if not {'color', 'cmap', 'edgecolor'} & set(plot_kwargs):
        plot_kwargs.update({'color': 'lightgrey', 'edgecolor':'black'})
    data_1D.plot(kind='bar', ax=ax, **plot_kwargs)
    ax.set_title(title)
    ax.grid(axis='y')
    ax.set_axisbelow(True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

def pivot_table(*pivot_args, **pivot_kwargs):
    """Make a pivot table with standard options"""
    pivot_options = {
        'aggfunc': 'size',
        'fill_value': 0,
    }
    pivot_options.update(pivot_kwargs)
    return pd.pivot_table(*pivot_args, **pivot_options)

def prop_table(count_data, sumi=1, divi=0):
    """Make a prop table from a table of counts"""
    if len(count_data.shape) == 1:
        return count_data / count_data.sum()
    else:
        return count_data.div(count_data.sum(sumi), divi)

def heatmap(df, **kwargs):
    """Draw seaborne heatmap with custom settings"""
    default_kwargs = {
        'center': 1.3,
        'cmap': sns.diverging_palette(220, 10, as_cmap=True),
        'square': True,
        'linewidth': 0.5,
        'annot': True,
    }
    default_kwargs.update(kwargs)
    sns.heatmap(df, **default_kwargs)

def highlight_max(s):
    """Highlight max value in a df column."""
    is_max = s == s.max()
    return ['color: red' if v else '' for v in is_max]
    
def max_highlighter(pr_df):
    """Show proportion dataframe with highlighting."""
    return pr_df.style.apply(highlight_max, 1)

class Analyze:
    """Run standard analysis on dataset"""
    def __init__(self, *pivot_args, fishers=True, **pivot_kwargs):
        """Initialize analysis.

        Args:
            data: Pandas dataframe of observations
            index: an index on which to form the pivot table
            columns: columns for the pivot table
            fishers: whether to run Fisher's analysis
            **pivot_kwargs: kwargs for pandas pivot function
        """
        # calculuate count table
        self.ct = pivot_table(*pivot_args, **pivot_kwargs)
        # sort the ct based on biggest sums
        sort_sums = self.ct.sum().sort_values(ascending=False).index
        self.ct = self.ct.loc[:,sort_sums]
        sort_sums2 = self.ct.sum(1).sort_values(ascending=False).index
        self.ct = self.ct.loc[sort_sums2]
        
        # calculate prop table with props across rows
        self.pr = prop_table(self.ct)
        # calculate prop table with props across columns (transcribed to rows to distinguish)
        self.pr2 = prop_table(self.ct.T)
        # calculate 1-in-N odds
        # see https://math.stackexchange.com/q/1469242
        self.oneN = 1 / self.pr
        self.odds = (1 / self.pr) - 1

        # run Fisher's collocation analysis
        self.fishers = fishers
        if fishers:
            self.fish, self.fish_odds = my_stats.apply_fishers(self.ct, 0, 1)
         
    def show(self):
        """Show results of the analysis"""
        
        print('counts:')
        display(self.ct)
        print()
        
        print('proportions 1:')
        display(max_highlighter(self.pr))
        print()

        print('proportions 2:')
        display(max_highlighter(self.pr2))
        print()

        if self.fishers:
            print('Fisher\'s test with log transform:')
    
            # check for inf values for plotting, since these 
            # cannot be plotted otherwise
            if np.inf in self.fish.values or -np.inf in self.fish.values:
                display(self.fish)
                sys.stderr.write("NB: Fisher's test (+/-)np.inf replaced with 300 for plotting")
                plotfish = self.fish.replace(np.inf, 300)
                plotfish = plotfish.replace(-np.inf, -300)
                heatmap(plotfish.round())
            else:
                heatmap(self.fish.round())

class AnalysisSet:
    """Easily run a set of analyses with results stored as attribs."""

    def __init__(self, *analyze_params):
        """Initialize and run a bunch of analyses."""
        print('setting up analyses...')
        for aparam in analyze_params:
            an_name, an_args, an_kwargs = aparam
            setattr(self, an_name, Analyze(*an_args, **an_kwargs))
        print('\tdone!')
        

# avoid deprecating old variables
def PivotProp(*args, **kwargs):
    return Analyze(*args, fishers=False, **kwargs)

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

class TextShower:
    """Show stylized text examples from a Dataframe of data"""
    
    def __init__(self, default=['ref', 'sentence', 'text_full'], 
                 stylize=['sentence', 'text_full']):
        self.default = default
        self.stylize = stylize

    def show(self, df, extra=[], spread=0):
        """Display text from pandas dataframe in a readable way."""
        original_shape = df.shape
        df = df[self.default + extra]
        if spread > 0:
            spread_i = get_spread(df.index, spread)
            df = df.loc[spread_i]
        print(f'showing {df.shape[0]} of {original_shape[0]}')
        return df.style.applymap(pretty_hebrew, subset=self.stylize)

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
