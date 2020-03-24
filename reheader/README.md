# Rename read name

A Galaxy tool to modify read name in fastq using python and Biopython

Currently produces fastq file with new read name with samplename

Requires:

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
