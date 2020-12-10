"""
Function for building CSV tables.
"""

import csv
import json
import collections
from pathlib import Path
from itertools import cycle

def csv_from_dict(rows, outpath):
    """Easy-write csv from list of dicts."""
    with open(outpath, 'w') as outfile:
        header = rows[0].keys()
        writer = csv.DictWriter(outfile, header)
        writer.writeheader()
        writer.writerows(rows)
   
def build_tables(table_params, samples):
    """Construct tables from a pairing of names and row makers.
    
    Args:
        table_params: dict with file paths as keys and a function
            for creating a row of data from a sample ("row maker").
        samples: objects to iterate over and calculate the rows
    Returns:
        list of dicts where each embedded dict is a row in the table
    """
    table_data = collections.defaultdict(list)
    for sample in samples:
        for table_path, row_maker in table_params.items(): 
            table_data[table_path].append(row_maker(sample)) 
    return table_data

def match_path2table(outpaths, rowmakers):
    """Match any equal/unequal number of paths with rowmaker functions.

    How to match arbitrary numbers of outpaths with finite 
    rowmaker functions?

    A simple zip won't account for when there are 
    more than 1 verbform in the sample set; we must therefore
    use itertools.cycle to ensure that the tables are repeated
    for every verbform.
    
    This works because the number of output files per verb match
    the number of tables defined (or at least they must!);
    so each verb cycles through all of the available table builds.
    """
    return dict(zip(outpaths, cycle(rowmakers)))

def build_sample_tables(rowmakers, sample_paths, outpaths):
    """Intake a path of samples and output tables for the project.
    
    Args:
        rowmakers: list of functions that return dicts of row data
            on a supplied sample
        sample_paths: paths to individual sample sets
        outdir: directory to export tables to
    Returns:
        Exports tables of calculated data.
    """
    table_params = match_path2table(outpaths, rowmakers)

    # iterate through all samples and build up the data
    for samp_set in sample_paths:
        set_file = Path(samp_set)    
        samples = sorted(json.loads(set_file.read_text()))

        # build the tables
        table_data = build_tables(table_params, samples)

        # make the exports
        for table_path, rows in table_data.items():
            csv_from_dict(rows, table_path)
