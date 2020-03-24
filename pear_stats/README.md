# Generate Paired-end reads merging statistic

A Galaxy tool to produce paired-end reads merging statistic using PEAR stat log

Currently produces a paired-end read merging statistic log 

Requires:

python 3.7


### Run reheader.py with fastq file
python pear_stats.py -i "test-data/F3D0.log,test-data/F3D1.log" -s "F3D0,F3D1" -o testout

## Version history:

**XML Wrapper:**

Alpha version by Michael Thang of QFAB, Australia.


**Python script: pear_stats.py**

* Original version

0.1.0: Michael Thang QFAB
