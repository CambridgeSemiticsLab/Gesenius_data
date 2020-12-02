"""
This module contains code used 
for analyzing collocational
tendencies of given constructions.
"""
import sys
import numpy as np
import pandas as pd
from df_styles import max_highlighter
from plotting import heatmap
from stats import significance as sig

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
    def __init__(self, *pivot_args, fishers=True, **pivot_kwargs):
        """Initialize analysis.

        For the Gesenius project we examine a number of 
        contexts that co-occur with given verbs. This
        class runs the analysis on a given dataset.

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
        # calculate prop table with props across columns (transposed to rows to distinguish)
        self.pr2 = prop_table(self.ct.T)
        # calculate 1-in-N odds
        # see https://math.stackexchange.com/q/1469242
        self.oneN = 1 / self.pr
        self.odds = (1 / self.pr) - 1

        # run Fisher's collocation analysis
        self.fishers = fishers
        if fishers:
            self.fish, self.fish_odds = sig.apply_fishers(self.ct, 0, 1)

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
