# Symmetric Plot for microbial differential expressed genes

A Galaxy tool to produce symmetric plots using count table, metadata table and green genes taxonomy file.

Currently produces the plots embedded in a html file for output with links to a PDF file.

Requires:

Phyloseq 1.22.3
r-getopt 1.20.0
r-doparallel 1.0.11
ghostscript 9.18


### Run symmetric_plot.R with three input files
Rscript 'symmetric_plot.r' --input.data='count.txt' --meta.data='metadata.txt' --obs.data='gg_13_5_taxonomy.txt' --taxrank='Phylum' --record='30' --norm="false" --n.column="3" --g.group="Early,Late" --outdir=outputdir --htmlfile="SymmetricPlot_out.html"


## Version history:

**XML Wrapper:**

Alpha version by Michael Thang of QFAB, Australia.

* Incorporated tests
* Requirements
* Version statement
* Citations

1.0.0: Michael Thang QFAB, Simon Gladman Melbourne Bioinformatics

* 3 input files are required 1) count.txt; 2) metadata.txt and 3) gg_13_5_taxonomy.txt
* Output symmetric plot


**R Script: symmetric_plot.r**

* Original version

1.0.0: Michael Thang QFAB

