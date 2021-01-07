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
from df_styles import get_spread

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
    def __init__(self, data, fishers=True, dp=True, **pivot_kwargs):
        """Initialize analysis.

        For the Gesenius project we examine a number of 
        contexts that co-occur with given verbs. This
        class runs the analysis on a given dataset.

        Args:
            pivot_args: arguments to feed to Pandas.pivot_table
            fishers: whether to run Fisher's analysis
            **pivot_kwargs: kwargs for pandas pivot function
        """
        # run multi-dimensional analyses
        if 'columns' in pivot_kwargs:
            # calculuate count table
            self.count = pivot_table(data, **pivot_kwargs)
            # sort the ct based on biggest sums
            sort_sums = self.count.sum().sort_values(ascending=False).index
            self.count = self.count.loc[:,sort_sums]

            # sums
            self.count_sum = pd.DataFrame(self.count.sum())
            self.count_sum1 = pd.DataFrame(self.count.sum(1))
            self.count_sum.columns = self.count_sum1.columns = ['sum']
            self.count_sum = self.count_sum.sort_values(by='sum', ascending=False)
            self.count_sum1 = self.count_sum1.sort_values(by='sum', ascending=False)

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
            if dp:
                self.dp = sig.apply_deltaP(self.count, 0, 1)
                self.dp2 = sig.apply_deltaP(self.count.T, 0, 1)

        # run one dimensional analyses
        else:
            self.count = pd.DataFrame(data[pivot_kwargs['index']].value_counts())
            self.count.columns = ['sum']
            self.count = self.count.sort_values(by='sum', ascending=False)
            self.prop = self.count / self.count.sum()
            self.sum = pd.DataFrame(self.count.sum())

def get_table_headers(df):
    """Identify table headers by indices and store in a dictionary."""
    get_indices = lambda axis: [int(i) for i in np.arange(0, axis.nlevels)]
    return {
        'index_col': get_indices(df.index),
        'header': get_indices(df.columns)
    }

def make_text_examples(df, ex_params):
    """Build copy and pastable text samples."""

    # execute query
    query = ex_params['query']
    df = df.query(ex_params['query'])
    n_results = df.shape[0]

    # return empty search results
    if n_results == 0:
        return [query, '0 results']

    # or process into examples:
    spread = ex_params.get('spread', 25)
    spread_i = get_spread(df.index, spread)
    df = df.loc[spread_i]

    # sort out texts
    # NB that esv and niv texts might also be similarly formatted later on
    bhs_joiner = ex_params.get('bhs_joiner', '')
    bhs_text = ex_params.get('bhs_text', ['clause_atom'])
    bhs_text = df[bhs_text].astype(str).agg(bhs_joiner.join, axis=1)
    
    exs = [query, f'{df.shape[0]} of {n_results}']

    for node in df.index:
        ref = df.loc[node]['ref_abbr']
        esv = df.loc[node]['esv'].lower()
        niv = df.loc[node]['niv'].lower()
        bhs = bhs_text[node]
            
        if niv == esv:
            ex = f'{esv} (ESV, NIV | BHS {bhs} {ref})'
        else:
            ex = f'{esv}, {niv} (ESV, NIV | BHS {bhs} {ref})'
        
        exs.append(ex)
        
    return exs

def run_analysis(params, outdir, round=2, table_headers={}):
    """Execute an analysis with specified parameters"""

    do_fishers = params.get('fishers', True),

    # execute the analysis using the provided parameters
    args = {'name', 'df', 'examples'}
    pivot_kwargs = {k:v for k,v in params.items() if k not in args}
    analysis = Analyze(
        params['df'], 
        **pivot_kwargs
    )

    # construct examples
    exs = {}
    for ex_param in params.get('examples', []):
       exs[ex_param['query']] = make_text_examples(params['df'], ex_param) 

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

    # export text examples
    ex_dir = outdir.joinpath('examples')
    for name, exs in exs.items():
        if not ex_dir.exists():
            ex_dir.mkdir()
        # text file for copy/paste
        textfile = ex_dir.joinpath(f'{name}.txt')
        textfile.write_text('\n'.join(exs))

        # TODO? df for visualizations
        
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
