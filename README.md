# Gesenius Data

<a href="https://www.ames.cam.ac.uk"><img src="https://github.com/CambridgeSemiticsLab/nena_tf/blob/master/docs/images/CambridgeU_color.jpg" width="236" height="49"></a>

Analyzing linguistic data from the Hebrew Bible in conjunction with Gesenius. The purpose of the repository is to produce datasets used for the ongoing project to revise Gesenius, led by Prof. Geoffrey Khan.

## Directory Structure and Data

Analysis and data-production pipeline occurs in `workflow` and is output to `results`. 

This project uses [Snakemake](https://snakemake.readthedocs.io) for the pipeline. 

## License

All data which is crucial for the analysis is stored openly (MIT license) in `results/csv`. The alignment data
and the English translations are not themselves open-source and thus cannot be released. However,
we do release the data derived from those sources secondarily. For English alignments, we provide
a link between a BHSA node [from Text-Fabric](https://github.com/etcbc/bhsa/_master_/tf/c) and a 
given tense-tagging, which has been composed using a mixture of Spacy and manually input rules. Thus
we cannot provide the English string, e.g. ברא = "he created" (ESV), but we can provide the tense tag 
which is the only thing that is crucial for the analysis anyways, e.g. ברא = "simple past" (ESV).
The LXX data is under its own license, coming from the CATSS project. We only provide a small subset
here. The full dataset can be found at http://ccat.sas.upenn.edu/gopher/text/religion/biblical/parallel/.

As with any academic work, please cite when using this repository.
