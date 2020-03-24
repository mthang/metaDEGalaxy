#!/usr/bin/env python
import sys
import argparse

parser = argparse.ArgumentParser(
    description="Parse multiple Pear statistic log to a tabular format\n" +
                "Example:\n python pear_stats.py -i \"file1.log,file2.log\" -s \"samplename1 samplename2\" -o outputfile")
parser.add_argument("-v","--version",action="version",version="%(prog)s 1.0")
parser.add_argument("-i","--input",dest="inputfilelist",default=False,help="a list of input file")
parser.add_argument("-s","--samplename", dest="samplename",default=False,help="a list of input filename")
parser.add_argument("-o","--outfile",dest="outputfile",default=False,help="Pear statistic output")


if(len(sys.argv) == 1):
       parser.print_help(sys.stderr)
       sys.exit()

args = parser.parse_args()

tags = ['Assembled reads','Discarded reads','Not assembled reads']
LINESTART=30
LINEEND  =LINESTART+2


inputfiles=args.inputfilelist.split(',')
inputfilenames=args.samplename.split(',')
outputfile=open(args.outputfile,'w')

allAssembled = 0

def processfile(instr):
	result=[]
	with open(instr,'r') as f:
		for linenum,line in enumerate(f):
			if LINESTART <= linenum <= LINEEND:
				ix = linenum-LINESTART
				if (line.startswith(tags[ix])):
					result.append(line.rstrip())
					if (ix == 0):
						token = line.strip().split('(')[1]
						token = token.replace("%)","")
						global allAssembled
						allAssembled += float(token)
					else:
						print("ARGH!:", line)
	return(result)

for element in range(0,len(inputfiles)):
    output=processfile(inputfiles[element])
    output.insert(0,inputfilenames[element])
    outputfile.write("\t".join(output))
    outputfile.write("\n")

averageAssembled = allAssembled / len(inputfiles)

averageAssembledOut=["The above assessment has been performed on 1000 randomly selected reads per sample file.\nAverage % of overlapping paired-end reads =",str(averageAssembled),"\nIf the average percentage is greater than 50%, you can consider using workflow 16S_biodiversity_for_overlap_PE.\nHowever, if the average percentage is less than 50%, use 16S_biodiversity_nonoverlap_PE."]


outputfile.write("\n\n\n")
outputfile.write("\t".join(averageAssembledOut))
outputfile.close()
