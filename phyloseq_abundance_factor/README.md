# Phyloseq Abundance Plot by factor

A Galaxy tool to produce abundance plot (by factor) using Phyloseq from either a BIOM1 file.

Currently produces the plots embedded in a html file for output with links to a PDF file.

Requires:

Phyloseq 1.22.3
r-getopt 1.20.0
r-doparallel 1.0.11
ghostscript 9.18



### Run phyloseq_abundance_factor.r with biom file
Rscript 'phyloseq_abundance_factor.r' --biomfile='test.biom' --metafile='metadata.txt' --xcolumn="5" --lcolumn="3" --factor1="4" --outdir="outputdir" --htmlfile='biom_out.html'

## Version history:

**XML Wrapper:**

Alpha version by Michael Thang of QFAB, Australia.


* Incorporated tests
* Requirements
* Version statement
* Citations


**R Script: phyloseq_abundance_factor.r**

* Original version

0.1.0: Michael Thang QFAB

