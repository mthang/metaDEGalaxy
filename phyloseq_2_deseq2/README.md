# Differential Expressed Analysis

A Galaxy tool to produce normalised count table and differential expressed genes table using DESeq2 from BIOM1 file and metadata file

Currently produces normalised count table and differential expressed genes table.

Requires:

DESeq2 1.14.1
Phyloseq 1.22.3
r-getopt 1.20.0

### Run phyloseq_2_deseq2.r with biom file
Rscript 'phyloseq_2_deseq2.r' --biomfile='test.biom' --metafile='metadata.txt' --factor="3" --test="Wald" --fitType="parametric" --cutoff=0.05 --result="DE_table.txt" --normalisedResult=""norm_DE_table.txt"


## Version history:

**XML Wrapper:**

Alpha version by Michael Thang of QFAB, Australia.

1.22.3.1

* Incorporated tests
* Requirements
* Version statement
* Citations


**R Script: phyloseq_2_deseq2.r**

* Original version

0.1.0: Michael Thang QFAB
