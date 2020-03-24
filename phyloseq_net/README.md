# Phloseq Ordination Plot

A Galaxy tool to produce network plot using Phyloseq from either a BIOM1 file or 2 input tables.

Currently produces the plots embedded in a html file for output with links to a PDF file.

Requires:

Phyloseq 1.22.3
r-getopt 1.20.0
r-doparallel 1.0.11
ghostscript 9.18


### Run phyloseq_ordination_plot.R with three input files
Rscript phyloseq_net.r --infile="count.txt" --metafile="metadata.txt" --obsfile="observation" --norm="false" --xcolumn="5" --lcolumn="4" --outdir="outputdir" --htmlfile="test.html"


### Run phyloseq_nmds.R with biom file
Rscript phyloseq_net.r --biom='GP.biom' --norm="true" --xcolumn="5" --lcolumn="4" --outdir="outputdir" --htmlfile=

## Version history:

**XML Wrapper:**

Alpha version by Michael Thang of QFAB, Australia.

1.22.3: Michael Thang of QFAB, Australia

* Incorporated tests
* Requirements
* Version statement
* Citations


**R Script: phyloseq_net.r**

* Original version

0.1.0: Michael Thang QFAB

