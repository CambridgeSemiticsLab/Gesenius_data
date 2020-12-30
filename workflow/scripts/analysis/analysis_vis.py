import os
import json
import collections
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# custom modules
from plotting import heatmap
from df_styles import df_highlighter

# there is a bug when using indices for rows/cols 
# when loading with pd.read_csv; see this discussion:
# https://stackoverflow.com/q/40659212
# we use warnings to silence this
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def plot_fishers(fishers_data, title=''):
    """Plot Fisher's exact correlation scores with a heatmap."""

    # try to infer the necessary figure size and set up axis
    nrows, ncols = fishers_data.shape
    figsize = (ncols, nrows)
    fig, ax = plt.subplots(figsize=figsize)

    # check for inf values for plotting, since these 
    # cannot be plotted otherwise
    if np.inf in fishers_data.values or -np.inf in fishers_data.values:
        fishers_data = fishers_data.replace(np.inf, 300)
        fishers_data = fishers_data.replace(-np.inf, -300)
        title = title + "\nNB: Fisher's test (+/-)np.inf replaced with 300 for plotting"

    # plot with specifications
    heatmap(
        fishers_data.round().astype(int), 
        ax=ax,
        #robust=True,
        fmt="d",
    )
    ax.set_title(title)
    ax.set_xlabel('')
    ax.set_ylabel('')

def table2html(df, filePath, stylesheet='', title='', highlight_rule='max'):
    """Export a table of data to HTML."""    
       
    html = """
<html>    
    <h1>{title}</h1>
    <head>
        <link rel="stylesheet" href="{stylesheet}">
    </head>
    <body>
        {data}
    </body>
</html>
    """.strip()
    table = df_highlighter(df, highlight_rule)
    table = table.set_precision(2)
    html = html.format(
        data=table.render(), 
        stylesheet=stylesheet,
        title=title,
    )

    # write to file
    filePath.write_text(html)

def visualize_analysis(analysis_dir, output_dir, table_headers, tablestyles=''):
    """Visualize the results of an analysis with HTML and SVG files

    Args:
        *TODO
    Returns:
        nothing. Outputs to a supplied output directory.
    """

    # set up output directories
    analysis_dir = Path(analysis_dir)
    outdir = Path(output_dir).joinpath(analysis_dir.name)
    if not outdir.exists():
        outdir.mkdir(parents=True)

    # output any csv files to HTML
    for csv_file in analysis_dir.glob('*.csv'):
        outfile = outdir.joinpath(f'{csv_file.stem}.html')
        headers = table_headers[csv_file.stem]
        df = pd.read_csv(csv_file, **headers)
        title = f'{analysis_dir.name}, {csv_file.stem}'
        highlight = 'fishers' if csv_file.name.startswith('fishers') else 'max'
        table2html(
            df, 
            outfile, 
            stylesheet=tablestyles, 
            title=title,
            highlight_rule=highlight
        )

        # special plots here
        if csv_file.name.startswith('fishers') and max(df.shape) < 50:
            outfile = outdir.joinpath('fishers_heatmap.svg')
            plot_fishers(df, title=analysis_dir.name)
            plt.savefig(outfile, bbox_inches='tight', format='svg')

def visualize_analyses(input_dir, output_dir, tablestyles):
    """Visualize a set of analyses."""

    # identify analyses by folder name, get their csv files
    # as a dataframe and output the data to HTML and SVG files 
    for analysis_path in Path(input_dir).glob('*'):

        # skip non-directory paths
        if not analysis_path.is_dir():
            continue

        # feed the analysis path into the visualizer which will
        # automatically identify the csv files as data to be visualized
        tablestyles = f'../../../../css/{Path(tablestyles).name}' # make tablestyles relative to the HTML file rather than snakemake
        tableheaders = json.loads(analysis_path.joinpath('table_headers.json').read_text())
        visualize_analysis(
            analysis_path, 
            output_dir,
            tableheaders,
            tablestyles=tablestyles
        )
