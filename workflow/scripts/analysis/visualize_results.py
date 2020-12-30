from pathlib import Path

# NB snakemake runs script from /workflow directory
sys.path.append('scripts/analysis')
from analysis_vis import visualize_analyses

# run analyses on input
visualize_analyses(snakemake.input.results, snakemake.output.dir, snakemake.input.tablestyles)
