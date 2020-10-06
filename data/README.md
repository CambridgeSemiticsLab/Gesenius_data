# A note on the need for open data

NB: unfortunately most of the data we are using for this project is closed source.
This is very unfortunate. However, in the interest of openness, we are providing
the code used to process that data here in this directory. 

Furthermore, not all of the data is completely unaccessible.

The LXX data is made accessible in the CATSS_parsers repository https://github.com/codykingham/CATSS_parsers
That repo contains a tool which downloads the LXX data and saves the parallel alignments and morpholoy
as JSON files to disk. That is directly the data used in the `lxx` directory here.

The modern translation data (ESV, NIV) is unfortunately quite closed at the moment. But 
as soon as it is completed, we would like to build in an open version of an English
translation alignment produced by Unfolding Word:
https://git.door43.org/unfoldingWord/en_ult

