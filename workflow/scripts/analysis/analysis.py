"""
This module contains code used 
for analyzing collocational
tendencies of given constructions.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from df_styles import max_highlighter
from plotting import heatmap
from stats import significance as sig
from pathlib import Path

def pivot_table(*pivot_args, **pivot_kwargs):
    """Make a pivot table with standard options"""
    pivot_options = { 
        'aggfunc': 'size',
        'fill_value': 0,
    }   
    pivot_options.update(pivot_kwargs)
    return pd.pivot_table(*pivot_args, **pivot_options)

def prop_table(count_data, sumi=1, divi=0):
    """Make a proportions table from a table of counts"""
    if len(count_data.shape) == 1:
        return count_data / count_data.sum()
    else:
        return count_data.div(count_data.sum(sumi), divi)

class Analyze:
    """Run standard analysis on verb dataset"""
    def __init__(self, data, fishers=True, **pivot_kwargs):
        """Initialize analysis.

        For the Gesenius project we examine a number of 
        contexts that co-occur with given verbs. This
        class runs the analysis on a given dataset.

        Args:
            pivot_args: arguments to feed to Pandas.pivot_table
            fishers: whether to run Fisher's analysis
            **pivot_kwargs: kwargs for pandas pivot function
        """
    
        # calculuate count table
        self.count = pivot_table(data, **pivot_kwargs)
        # sort the ct based on biggest sums
        sort_sums = self.count.sum().sort_values(ascending=False).index
        self.count = self.count.loc[:,sort_sums]
        sort_sums2 = self.count.sum(1).sort_values(ascending=False).index
        self.count = self.count.loc[sort_sums2]

        # calculate prop table with props across rows
        self.prop = prop_table(self.count)
        # calculate prop table with props across columns (transposed to rows to distinguish)
        self.prop2 = prop_table(self.count.T)
        # calculate 1-in-N odds
        # see https://math.stackexchange.com/q/1469242
        self.oneN = 1 / self.prop
        self.odds = (1 / self.prop) - 1

        # run Fisher's collocation analysis
        if fishers:
            self.fishers, self.fishers_odds = sig.apply_fishers(self.count, 0, 1)

def plot_fishers(fishers_data):
    """Plot Fisher's exact correlation scores with a heatmap."""

    # check for inf values for plotting, since these 
    # cannot be plotted otherwise
    if np.inf in fishers_data.values or -np.inf in fishers_data.values:
        msg = "NB: Fisher's test (+/-)np.inf replaced with 300 for plotting"
        plotfish = fishers_data.replace(np.inf, 300)
        plotfish = plotfish.replace(-np.inf, -300)
        heatmap(plotfish.round())
        plt.title(msg)
    else:
        heatmap(fishers_data.round())

def export_table(df, filename, styles):
    """Export a table of data to HTML."""    
    df.to_html(filename)

def run_analysis(params, outdir):
    """Execute an analysis with specified parameters"""

    do_fishers = params.get('fishers', True),

    # execute the analysis using the provided parameters
    args = {'name', 'df'}
    pivot_kwargs = {k:v for k,v in params.items() if k not in args}
    analysis = Analyze(
        params['df'], 
        **pivot_kwargs
    )

    # prepare paths for exports
    outdir = Path(outdir).joinpath(params['name'])
    if not outdir.exists():
        outdir.mkdir(parents=True)

    # export the data tables
    for name, df in analysis.__dict__.items():
        outfile = outdir.joinpath(f'{name}.html')
        export_table(df, outfile, params.get('table_styles', {}))

    # export the graphs
    if do_fishers:
        outfile = outdir.joinpath('fishers_heatmap.svg')
        plot_fishers(analysis.fishers)
        plt.savefig(outfile, bbox_inches='tight', format='svg')
