"""
Function for building CSV tables.
"""

import csv
import json
import collections
from pathlib import Path
from itertools import cycle, product

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
    outpaths = sorted(outpaths) # clusters paths with single verbs 
    path2rows = dict(zip(outpaths, cycle(rowmakers)))
    return path2rows

def build_sample_tables(rowmakers, sample_paths, outpaths):
    """Build tables given a set of rowmaker functions, inpaths, and outpaths.

    The number of unique rowmakers must equal the number of outpath
    names in order to ensure the proper alignment between the two.
    
    ! NOTE !
    The order of rowmakers MUST match the order of outpaths. Otherwise
    there will be a misalignment. It is not possible to check this 
    automatically so take care to ensure alignment!

    Args:
        rowmakers: list of functions that return dicts of row data
            on a supplied sample
        sample_paths: paths to individual json files containing samples
            in the form of a set
        outpaths: snakemake output object (iterable)
    Returns:
        Exports tables of calculated data.
    """

    # sanity check for misalignments of outpaths and rowmakers
    if len(rowmakers) != len(outpaths._names):
        raise Exception(
            "The number of rowmakers does not == number of output names! "
            "Thus alignment is impossible."
        )

    # map outpaths to their relevant rowmakers
    # NB that the order of rowmakers is crucial here!
    table_params = match_path2table(outpaths, rowmakers)

    # get a list of mappings between sample paths and outpaths;
    # solves a similar problem to mapping outpaths to rowmakers
    # namely that there is twice more outpaths than inpaths
    # and the outpaths need to be aligned in a cyclic way to
    # the various verbforms
    samp2out = list(zip(cycle(sample_paths), outpaths))

    # iterate through all samples and build up the data
    for samp_set in sample_paths:
        set_file = Path(samp_set)    
        samples = sorted(json.loads(set_file.read_text()))

        # retrieve only those outpath 2 rowmaker mappings 
        # that are relevant to this sample; this is done by 
        # crossreferencing the two outfile mappings
        samp_params = {}
        for samp_path, outpath in samp2out:
            if samp_path == samp_set:
                samp_params[outpath] = table_params[outpath]

        # build the tables
        table_data = build_tables(samp_params, samples)

        # make the exports
        for table_path, rows in table_data.items():
            csv_from_dict(rows, table_path)
