#!/usr/bin/env python 
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from os.path import basename
import os
import re
import argparse

parser = argparse.ArgumentParser(
    description="reformat the read name (header) by appending the sample name to the read name\n" + 
				"Example:\n python reheader.py -n F3D0_R1.fastq -i test-data/F3D0_R1.fastq -o test-data/test -l mylog -d test-data/")
parser.add_argument("-v","--version",action="version",version="%(prog)s 1.0")
parser.add_argument("-n","--samplename",dest="samplename",default=False,help="input sample name")
parser.add_argument("-i","--input",dest="inputfile",default=False,help="input filename in FASTQ format")
parser.add_argument("-l","--log", dest="logfile",default=False,help="output log file")
parser.add_argument("-o","--outfile",dest="outputfile",default=False,help="output filename")
parser.add_argument("-d","--outdir",dest="outputdir",default=False,help="output directory")


if(len(sys.argv) == 1):
       parser.print_help(sys.stderr)
       sys.exit()

args = parser.parse_args()


filename = args.samplename
infile = args.inputfile
str_to_add = os.path.splitext(basename(filename))[0]
outfile = args.outputfile
outdir = args.outputdir
logfile = args.logfile


rdict = {
    '_R1': '/1',
    '_R2': '/2',
    '_1': '/1',
    '_2': '/2',
}

rdict_remove = {
    '_R1': '',
    '_R2': '',
    '_1': '',
    '_2': '',
}

def makesubs(s):
    for pattern, repl in rdict.items():
        pat1 = pattern +'_?[A-Za-z0-9]+$'
        pat2 = pattern
        combined_pat = r'|'.join((pat1, pat2))
        s = re.sub(combined_pat, repl,s)
    return s

def makesubs_remove(s):
    for pattern, repl in rdict_remove.items():
        pat1 = pattern +'_?[A-Za-z0-9]+$'
        pat2 = pattern
        combined_pat = r'|'.join((pat1, pat2))
        s = re.sub(combined_pat, repl,s)
    return s

def appendStringToSequenceHeader(inputfile,header_to_add):
    records=[]
    for seq_record in SeqIO.parse(inputfile, "fastq"):
        header =seq_record.id
        header = "{0}".format(header) + "_" +header_to_add
        record = SeqRecord(seq_record.seq,id=header,description="")
        record.letter_annotations["phred_quality"]=seq_record.letter_annotations["phred_quality"]
        records.append(record)
    return records

str_to_search = makesubs_remove(str_to_add)
str_to_add = makesubs(str_to_add)
final_records=[]
outlogfile=open(os.path.join(outdir,logfile),"w")

final_records=appendStringToSequenceHeader(infile,str_to_add)
outlogfile.write(str_to_search)
SeqIO.write(final_records, outfile , "fastq")
outlogfile.close()
