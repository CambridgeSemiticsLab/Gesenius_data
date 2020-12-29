import os
import json
from pathlib import Path

# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from analysis_vis import visualize_analysis

# identify analyses by folder name, get their csv files
# as a dataframe and output the data to HTML and SVG files 
for analysis_path in Path(snakemake.input.results).glob('*'):

    # skip non-directory paths
    if not analysis_path.is_dir():
        continue

    # feed the analysis path into the visualizer which will
    # automatically identify the csv files as data to be visualized
    tablestyles = f'../{Path(snakemake.input.tablestyles).name}'
    tableheaders = json.loads(analysis_path.joinpath('table_headers.json').read_text())
    visualize_analysis(
        analysis_path, 
        snakemake.output.dir, 
        tableheaders,
        tablestyles=tablestyles
    )

def build_menu(dir, html_lines, level=0):
    """Construct an HTML menu of all analyses."""
    for item in sorted(dir.glob('*')):
        indent = '&nbsp;'* 8 * level
        link = f'<a href="{item.absolute()}" target="_blank">{item.stem}</a>'
        html_lines.append(f'{indent}{link}')
        if item.is_dir():
            build_menu(item, html_lines, level=level+1)

# construct HTML menu item
out_dir = Path(snakemake.output.dir)
menu = '<html>{data}</html>'
html_lines = []
build_menu(out_dir, html_lines)
menu_html = menu.format(data='\n<br>'.join(html_lines))
out_dir.joinpath('menu.html').write_text(menu_html)

# copy the stylesheet into the analysis directory
os.system(f"cp {snakemake.input.tablestyles} {snakemake.output.dir}/.")
