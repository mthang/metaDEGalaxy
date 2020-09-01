# Rename read identifiter
reheader.py script is written in python and Biopython is used to handle fastq file format. This script is used to rename the read identifier in the fastq file by appending the sample name to the end of the read identifier. The wrapper (xml) file is include for the Galaxy use. 

The output is in fastq file format with the sample name added to the read identifier.

Requirement:

python 3.7
Biopython 1.74

### Run reheader.py with fastq file
python reheader.py -n F3D0_R1.fastq -i test-data/F3D0_R1.fastq -o test-data/test -l log -d test-data/

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
