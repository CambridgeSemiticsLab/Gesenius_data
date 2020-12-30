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

def build_sample_tables(rowmakers, sample_path, outpaths):
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
        sample_path: path to a json file containing samples
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

    # align outpaths to their relevant rowmakers
    # NB that the order of rowmakers is crucial here!
    table_params = dict(zip(outpaths, rowmakers))

    # load the sample 
    sample_file = Path(sample_path)    
    samples = sorted(json.loads(sample_file.read_text()))

    # build the tables
    table_data = build_tables(table_params, samples)

    # make the exports
    for table_path, rows in table_data.items():
        csv_from_dict(rows, table_path)
