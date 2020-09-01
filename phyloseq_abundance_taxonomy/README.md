# Phyloseq Abundance Taxonomy  Plot

A Galaxy tool to produce an abundance taxonomy plot using a R package Phyloseq from either a BIOM1 file or 2 input tables.

This R script generates a plot embedded in a html file.

Requirement:

Phyloseq 1.22.3
r-getopt 1.20.0
r-doparallel 1.0.11
ghostscript 9.18


### Run phyloseq_nmds.R with biom file
Rscript 'phyloseq_abundance_taxonomy.r' --biomfile="test.biom" --metafile="metadata.txt" --xcolumn="5" --lcolumn="4" --taxonomy="Phylum" --outdir="outputdir" --htmlfile="biom_out.html"


## Version history:

**XML Wrapper:**

Alpha version by Michael Thang of QFAB, Australia.

1.22.3: Simon Gladman Melbourne Bioinformatics

* Incorporated tests
* Requirements
* Version statement
* Citations


**R Script: phyloseq_abundance_taxonomy**

* Original version

0.1.0: Michael Thang QFAB

