# Convert UCLUST format to tabular format 
uclust2otutable.py script is written in python and it has been modified to generate output file in tabular format. The original script can be found on the USEARCH website https://www.drive5.com/python/. The description of the script is described on this website https://www.drive5.com/python/uc2otutab_py.html.

The script takes the UC file format as an input and transform it to tabular format.

UC output file description can be found in this link https://www.drive5.com/usearch/manual/opt_uc.html.

# Field Description in UC file
* Tab-separated fields:
* 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=Label
* Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NotMatched
* For C and D types, PctId is average id with seed.
* QueryStart and SeedStart are zero-based relative to start of sequence.
* If minus strand, SeedStart is relative to reverse-complemented seed.


Requirement:

python 3.7

### Run uclust2otutable.py with uc file
python uclust2otutable.py -i uc_filename -o uc_table.txt

## Version history:

**XML Wrapper:**

Alpha version by Michael Thang of QFAB, Australia.

1.22.3.1

* Incorporated tests
* Requirements
* Version statement
* Citations


**Python script: reheader.py**

* Original version

0.1.0: Michael Thang QFAB
