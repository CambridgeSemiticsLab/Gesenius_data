"""
This module contains code used 
for analyzing collocational
tendencies of given constructions.
"""

import json
import re
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
            sort_sums2 = self.count.sum(1).sort_values(ascending=False).index
            self.count = self.count.loc[sort_sums2]

            # add counts with total column
            total = self.count.sum(1)
            total.name = 'TOTAL'
            self.count_total = pd.concat([total, self.count], 1)

            # add counts with total column
            total2 = self.count.T.sum(1)
            total2.name = 'TOTAL'
            self.count_total2 = pd.concat([total2, self.count.T], 1)

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
        elif 'index' in pivot_kwargs:
            self.count = pd.DataFrame(data[pivot_kwargs['index']].value_counts())
            self.count.columns = ['sum']
            self.count = self.count.sort_values(by='sum', ascending=False)
            self.prop = self.count / self.count.sum()
            self.sum = pd.DataFrame(self.count.sum())

        # set up an empty analysis object which can 
        # be subsequently filled with attribute assigns
        else:
            pass

def get_table_headers(df):
    """Identify table headers by indices and store in a dictionary."""
    get_indices = lambda axis: [int(i) for i in np.arange(0, axis.nlevels)]
    return {
        'index_col': get_indices(df.index),
        'header': get_indices(df.columns)
    }

def add_hit_symbols(text, search, symbol='_'):
    """Use re to find matches and insert hit symbols."""
    find = f'({search})'
    replace = f'{symbol}\g<1>{symbol}'
    return re.sub(find, replace, text)

def clean_puncts(text):
    """Fix space-separated punctuations."""
    text = re.sub(' ([,.?!;”:’])', '\g<1>', text)
    text = re.sub('([“‘]) ', '\g<1>', text)
    return text
    
def clean_roman_nums(text):
    """Replace Roman numerals in the verse ref."""
    text = re.sub('^I ', '1 ', text)
    text = re.sub('^II ', '2 ', text)
    return text

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
    if spread != -1:
        spread_i = get_spread(df.index, spread)
        df = df.loc[spread_i]

    # sort out texts
    # NB that esv and niv texts might also be similarly formatted later on
    bhs_joiner = ex_params.get('bhs_joiner', '')
    bhs_text = ex_params.get('bhs_text', ['verse'])
    bhs_text = df[bhs_text].astype(str).agg(bhs_joiner.join, axis=1)
    
    exs = [query, f'{df.shape[0]} of {n_results}\n']

    for i, node in enumerate(df.index):
        ref = df.loc[node]['ref_abbr']
        esv = df.loc[node]['esv']
        niv = df.loc[node]['niv']
        esv_verse = df.loc[node]['esv_verse']
        niv_verse = df.loc[node]['niv_verse']
        bhs = bhs_text[node].strip()
        heb = df.loc[node]['text_full']

        # add formatting to verse texts
        #ref = clean_roman_nums(ref)
        esv_verse = add_hit_symbols(esv_verse, esv)
        niv_verse = add_hit_symbols(niv_verse, niv)
        bhs = add_hit_symbols(bhs, heb)
        esv_verse = clean_puncts(esv_verse)
        niv_verse = clean_puncts(niv_verse)

        # build and add example text
        parts = [
            f'(x)\t{bhs} ({ref})',
            f'NIV\t{niv_verse}',
            f'ESV\t{esv_verse}'
        ]
        # add any extra text data 
        for name, et in ex_params.get('extra_text', {}).items():
            text = df.loc[node][et]
            parts.append(f'{name}\t{text}')

        ex = '\n'.join(parts) + '\n'
        exs.append(ex)
        
    return exs

def run_analysis(params, outdir, round=2, table_headers={}):
    """Execute an analysis with specified parameters"""

    # -- Run analysis on data tables --

    analysis_data = {}

    if 'index' in params:
        # execute the analysis using the provided parameters
        pivot_kwargs = {'index', 'columns', 'pivot_kwargs'}
        pivot_kwargs = {k:v for k,v in params.items() if k in pivot_kwargs}
        analysis = Analyze(
            params['df'], 
            **pivot_kwargs
        )

        # save analysis
        for name, df in analysis.__dict__.items():
            analysis_data[name] = df

   # run any special analyses
    for sa in params.get('special', []):

        # determine dataframe
        df = sa['df']
        if type(df) == str:
            df = analysis_data[df] 

        # run the special analysis and store as attribute
        do_funct = sa['do']
        sa_do = do_funct(df, **sa.get('kwargs', {}))
        for key, val in sa_do.items():
            analysis_data[key] = val

    # prepare paths for exports
    outdir = Path(outdir).joinpath(params['name'])
    if not outdir.exists():
        outdir.mkdir(parents=True)

    # export tables
    for name, df in analysis_data.items():
        outfile = outdir.joinpath(f'{name}.csv')
        df = df.round(round)
        table_headers[name] = get_table_headers(df) # save table header params
        df.to_csv(str(outfile), index=True)

    # -- Extract text examples --

    # construct examples
    exs = {}
    for ex_param in params.get('examples', []):
        example_df = ex_param.get('df', params['df']) 
        exs[ex_param['query']] = make_text_examples(example_df, ex_param) 

    # export text examples
    ex_dir = outdir.joinpath('examples')
    for name, exs in exs.items():
        if not ex_dir.exists():
            ex_dir.mkdir()
        # text file for copy/paste
        textfile = ex_dir.joinpath(f'{name}.txt')
        textfile.write_text('\n'.join(exs))

def run_analyses(analysis_params, out_dir):
    """Execute a set of analyses."""

    out_dir = Path(out_dir)

    # run all of the analyses in a loop
    for i, analysis in enumerate(analysis_params):

        # number the analyses
        i += 1
        name = analysis['name']
        analysis['name'] = f'{i}_{name}'

        # track parameters for the tables
        tableheaders = {}  

        # run analysis and produce tables
        run_analysis(analysis, out_dir, table_headers=tableheaders)

        # save table parameters in the analysis directory
        headersfile = out_dir.joinpath(analysis['name']).joinpath('table_headers.json')
        with open(headersfile, 'w') as outfile:
            json.dump(tableheaders, outfile)
