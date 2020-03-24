# Phyloseq Richness  Plot

A Galaxy tool to produce Richness plots using Phyloseq from either a BIOM1 file or 2 input tables.

Currently produces the plots embedded in a html file for output with links to a PDF file.

Requires:

Phyloseq 1.22.3
r-getopt 1.20.0
r-doparallel 1.0.11
ghostscript 9.18

### Run phyloseq_nmds.R with biom file
Rscript 'phyloseq_richness.r' --biomfile="test.biom" --metafile="metadata.txt"--xcolumn="5" --lcolumn="4" --outdir="outputdir" --htmlfile="biom_out.html"


## Version history:

**XML Wrapper:**

Alpha version by Michael Thang of QFAB, Australia.

* Incorporated tests
* Requirements
* Version statement
* Citations


**R Script: phyloseq_richness.R**

* Original version

0.1.0: Michael Thang QFAB

