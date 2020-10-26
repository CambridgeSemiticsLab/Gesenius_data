# Gesenius Data

Analyzing linguistic data from the Hebrew Bible in conjunction with Gesenius. The purpose of the repository is to produce datasets used for the ongoing project to revise Gesenius, led by Prof. Geoffrey Khan.

## Directory Structure and Data

Cambridge Semitics Lab strives to promote open science and open access wherever possible. Unfortunately, not all of the data we rely on in this repository is open. We acknowledge that this situation poses problems for reproducibility, and we commit to switch to open forms of data wherever possible, as soon as possible. 

Despite some data remaining private, we will always publish the scripts used to process that data. These scripts can be found in the `data` directory. Additionaly, all of the analysis of the data is published under `analysis`.

The .csv files used for the analyses are stored under two subdirectories:

* `data/_private_` - this directory is not pushed to Github as it contains non-open data. This includes translation alignments (NIV, ESV, KJV) from Global Bible Initiative. The LXX data stored here may be accessed separately by by running scripts found in https://github.com/codykingham/CATSS_parsers
* `data/_public_` - open data, including data based on BHSA (http://github.com/etcbc/bhsa)
