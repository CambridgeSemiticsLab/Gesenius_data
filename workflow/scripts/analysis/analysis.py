"""
This module contains code used 
for analyzing collocational
tendencies of given constructions.
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
            self.fishers, self.odds_fishers = sig.apply_fishers(self.count, 0, 1)

def get_table_headers(df):
    """Identify table headers by indices and store in a dictionary."""
    get_indices = lambda axis: [int(i) for i in np.arange(0, axis.nlevels)]
    return {
        'index_col': get_indices(df.index),
        'header': get_indices(df.columns)
    }

def run_analysis(params, outdir, round=2, table_headers={}):
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
        outfile = outdir.joinpath(f'{name}.csv')
        df = df.round(round)
        table_headers[name] = get_table_headers(df) # save table header params
        df.to_csv(str(outfile), index=True)

def run_analyses(analysis_params, out_dir):
    """Execute a set of analyses."""

    out_dir = Path(out_dir)

    # run all of the analyses in a loop
    for analysis in analysis_params:

        # track parameters for the tables
        tableheaders = {}  

        # run analysis and produce tables
        run_analysis(analysis, out_dir, table_headers=tableheaders)

        # save table parameters in the analysis directory
        headersfile = out_dir.joinpath(analysis['name']).joinpath('table_headers.json')
        with open(headersfile, 'w') as outfile:
            json.dump(tableheaders, outfile)
