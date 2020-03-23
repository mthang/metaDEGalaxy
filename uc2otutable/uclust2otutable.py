import sys
import progress
import subprocess
import tempfile
import traceback
import argparse

parser = argparse.ArgumentParser(
    description="This script converts uclust format from vsearch to tabular format"
	)
parser.add_argument("-v","--version",action="version",version="%(prog)s 1.0")
parser.add_argument("-i","--input",dest="uclust",default=False,help="input filename in uclust format")
parser.add_argument("-o","--output",dest="otutable",default=False,help="output filename")


if(len(sys.argv) == 1):
	parser.print_help(sys.stderr)
	sys.exit()

args = parser.parse_args()

ucFileName = args.uclust
outFileName = args.otutable


# Tab-separated fields:
# 1=Type, 2=ClusterNr, 3=SeqLength or ClusterSize, 4=PctId, 5=Strand, 6=QueryStart, 7=SeedStart, 8=Alignment, 9=Label
# Record types (field 1): L=LibSeed, S=NewSeed, H=Hit, R=Reject, D=LibCluster, C=NewCluster, N=NotMatched
# For C and D types, PctId is average id with seed.
# QueryStart and SeedStart are zero-based relative to start of sequence.
# If minus strand, SeedStart is relative to reverse-complemented seed.

MaxError = -1

Type = '?'
ClusterNr = -1
Size = -1
PctId = -1.0
LocalScore = -1.0
Evalue = -1.0
Strand = '.'
QueryStart = -1
SeedStart = -1
Alignment = ""
QueryLabel = ""
TargetLabel = ""
FileName = "?"
Line = ""

TRUNC_LABELS=0

def GetSampleId(Label):
	sep=";"
	SampleID_temp = Label.split(sep,1)[0]
	SampleID = SampleID_temp.split('_',1)[-1]
	return SampleID

def OnRec():
	global OTUs, Samples, OTUTable
	if Type != 'H':
		return

	OTUId = TargetLabel
	if OTUId not in OTUIds:
		OTUIds.append(OTUId)
		OTUTable[OTUId] = {}

	SampleId = GetSampleId(QueryLabel)
	if SampleId not in SampleIds:
		SampleIds.append(SampleId)

	N = GetSizeFromLabel(QueryLabel, 1)
	try:
		OTUTable[OTUId][SampleId] += N
	except:
		OTUTable[OTUId][SampleId] = N

def Die(Msg):
	print >> sys.stderr
	print >> sys.stderr

	traceback.print_stack()
	s = ""
	for i in range(0, len(sys.argv)):
		if i > 0:
			s += " "
		s += sys.argv[i]
	print >> sys.stderr, s
	print >> sys.stderr, "**ERROR**", Msg
	print >> sys.stderr
	print >> sys.stderr
	sys.exit(1)
	print("NOTHERE!!")
	
def Warning(Msg):
	print >> sys.stderr
	print >> sys.stderr, sys.argv
	print >> sys.stderr, "**WARNING**", Msg

def isgap(c):
	return c == '-' or c == '.'

def GetSeqCount(FileName):
	Tmp = tempfile.TemporaryFile()
	try:
		TmpFile = Tmp.file
	except:
		TmpFile = Tmp
	s = subprocess.call([ "grep", "-c", "^>", FileName ], stdout=TmpFile)
	TmpFile.seek(0)
	s = TmpFile.read()
	return int(s)

def GetSeqsDict(FileName):
	return ReadSeqsFast(FileName, False)

def ReadSeqsDict(FileName, Progress = False):
	return ReadSeqsFast(FileName, Progress)

def ReadSeqsOnSeq(FileName, OnSeq, Progress = False):
	ReadSeqs3(FileName, OnSeq, Progress)

def ReadSeqsFastFile(File, Progress = False):
	Seqs = {}
	Id = ""
	N = 0
	while 1:
		if N%10000 == 0 and Progress:
			sys.stderr.write("%u seqs\r" % (N))
		Line = File.readline()
		if len(Line) == 0:
			if Progress:
				sys.stderr.write("%u seqs\n" % (N))
			return Seqs
		if len(Line) == 0:
			continue
		Line = Line.strip()
		if Line[0] == ">":
			N += 1
			Id = Line[1:]
			if TRUNC_LABELS:
				Id = Id.split()[0]
			Seqs[Id] = ""
		else:
			if Id == "":
				Die("FASTA file does not start with '>'")
			Seqs[Id] = Seqs[Id] + Line

def ReadSeqsFast(FileName, Progress = True):
	File = open(FileName)
	return ReadSeqsFastFile(File, Progress)

def ReadSeqs(FileName, toupper=False, stripgaps=False, Progress=False):
	if not toupper and not stripgaps:
		return ReadSeqsFast(FileName, False)

	Seqs = {}
	Id = ""
	File = open(FileName)
	while 1:
		Line = File.readline()
		if len(Line) == 0:
			return Seqs
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == ">":
			Id = Line[1:]
			if TRUNC_LABELS:
				Id = Id.split()[0]
			if Id in Seqs.keys():
				Die("Duplicate id '%s' in '%s'" % (Id, FileName))
			Seqs[Id] = ""
		else:
			if Id == "":
				Die("FASTA file '%s' does not start with '>'" % FileName)
			if toupper:
				Line = Line.upper()
			if stripgaps:
				Line = Line.replace("-", "")
				Line = Line.replace(".", "")
			Seqs[Id] = Seqs[Id] + Line

def ReadSeqs2(FileName, ShowProgress = True):
	Seqs = []
	Labels = []
	File = open(FileName)
	if ShowProgress:
		progress.InitFile(File, FileName)
	while 1:
		progress.File()
		Line = File.readline()
		if len(Line) == 0:
			if ShowProgress:
				print >> sys.stderr, "\n"
			return Labels, Seqs
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == ">":
			Id = Line[1:]
			if TRUNC_LABELS:
				Id = Id.split()[0]
			Labels.append(Id)
			Seqs.append("")
		else:
			i = len(Seqs)-1
			Seqs[i] = Seqs[i] + Line

def ReadSeqs3(FileName, OnSeq, ShowProgress = True):
	File = open(FileName)
	if ShowProgress:
		progress.InitFile(File, FileName)
	Label = ""
	Seq = ""
	while 1:
		Line = File.readline()
		if len(Line) == 0:
			if Seq != "":
				OnSeq(Label, Seq)
			if ShowProgress:
				print >> sys.stderr, "\n"
			return
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == ">":
			if Seq != "":
				if ShowProgress:
					progress.File()
				if TRUNC_LABELS:
					Label = Label.split()[0]
				OnSeq(Label, Seq)
			Label = Line[1:]
			Seq = ""
		else:
			Seq += Line

def WriteSeq(File, Seq, Label = ""):
	if Label != "":
		print >> File, ">" + Label
	BLOCKLENGTH = 80
	SeqLength = len(Seq)
	BlockCount = int((SeqLength + (BLOCKLENGTH-1))/BLOCKLENGTH)
	for BlockIndex in range(0, BlockCount):
		Block = Seq[BlockIndex*BLOCKLENGTH:]
		Block = Block[:BLOCKLENGTH]
		print >> File, Block

def GetSizeFromLabel(Label, Default = -1):
	Fields = Label.split(";")
	for Field in Fields:
		if Field.startswith("size="):
			return int(Field[5:])
	if Default == -1:
		Die("Missing size >" + Label)
	return Default

def StripSizeFromLabel(Label):
	Fields = Label.split(";")
	NewLabel = ""
	for Field in Fields:
		if Field.startswith("size="):
			continue
		if NewLabel != "":
			NewLabel += ";"
		NewLabel += Field
	return NewLabel

def GetQualFromLabel(Label):
	n = Label.find("qual=")
	assert n >= 0
	return Label[n+5:-1]

def StripQualFromLabel(Label):
	n = Label.find("qual=")
	assert n >= 0
	return Label[:n]

def GetField(Label, Name, Default):
	Fields = Label.split(';')
	for Field in Fields:
		if Field.startswith(Name + "="):
			n = len(Name) + 1
			return Field[n:]
	if Default == "":
		Die("Field %s= not found in >%s" % (Name, Label))
	return Default

def GetIntFieldFromLabel(Label, Name, Default):
	return int(GetField(Label, Name, Default))

def GetFieldFromLabel(Label, Name, Default):
	return GetField(Label, Name, Default)

def DeleteFieldFromLabel(Label, Name):
	NewLabel = ""
	Fields = Label.split(';')
	for Field in Fields:
		if len(Field) > 0 and not Field.startswith(Name + "="):
			NewLabel += Field + ';'
	return NewLabel

def ReplaceSize(Label, Size):
	Fields = Label.split(";")
	NewLabel = ""
	Done = False
	for Field in Fields:
		if Field.startswith("size="):
			NewLabel += "size=%u;" % Size
			Done = True
		else:
			if Field != "":
				NewLabel += Field + ";"
	if not Done:
		die.Die("size= not found in >" + Label)
	return NewLabel	

def Error(s):
	print >> sys.stderr, "*** ERROR ***", s, sys.argv
	sys.exit(1)	

def ProgressFile(File, FileSize):
#	if not sys.stderr.isatty():
#	return
	Pos = File.tell()
	Pct = (100.0*Pos)/FileSize
	Str = "%s %5.1f%%\r" % (FileName, Pct)
	sys.stderr.write(Str)

def Progress(i, N):
#	if not sys.stderr.isatty():
	return
	Pct = (100.0*i)/N
	Str = "%5.1f%%\r" % Pct
	sys.stderr.write(Str)

def PrintLine():
	print(Line)

def ParseRec(Line):
	global Type
	global ClusterNr
	global Size
	global PctId
	global Strand
	global QueryStart
	global SeedStart
	global Alignment
	global QueryLabel
	global TargetLabel
	global LocalScore
	global Evalue
	
	Fields = Line.split("\t")
	N = len(Fields)
	if N != 9 and N != 10:
		Error("Expected 9 or 10 fields in .uc record, got: " + Line)
	Type = Fields[0]
	
	try:
		ClusterNr = int(Fields[1])
	except:
		ClusterNr = -1
		
	try:	
		Size = int(Fields[2])
	except:
		Size = -1

	Fields2 = Fields[3].split('/')
	LocalScore = -1.0
	Evalue = -1.0
	if len(Fields2) == 3:
		try:
			PctId = float(Fields2[0])
			LocalScore = float(Fields2[1])
			Evalue = float(Fields2[2])
		except:
			PctId = -1.0
	else:
		try:
			PctId = float(Fields[3])
		except:
			PctId = -1.0

	Strand = Fields[4]
	
	try:
		QueryStart = int(Fields[5])
	except:
		QueryStart = -1

	try:
		SeedStart = int(Fields[6])
	except:
		SeedStart = -1

	Alignment = Fields[7]
	QueryLabel = Fields[8]
	
	if len(Fields) > 9:
		TargetLabel = Fields[9]

def GetRec(File, OnRecord):
	global Line
	while 1:
		Line = File.readline()
		if len(Line) == 0:
			return 0
		if Line[0] == '#':
			continue
		Line = Line.strip()
		if len(Line) == 0:
			return 1
		ParseRec(Line)
		Ok = OnRecord()
		if Ok != None and Ok == 0:
			return 0
		return 1

def ReadRecs(argFileName, OnRecord, ShowProgress = False):
	return ReadFile(argFileName, OnRecord, ShowProgress)

def ReadRecsOnRec(argFileName, OnRecord, ShowProgress = True):
	return ReadFile(argFileName, OnRecord, ShowProgress)

def GetRecs(argFileName, OnRecord, ShowProgress = True):
	return ReadFile(argFileName, OnRecord, ShowProgress)

def ReadFile(argFileName, OnRecord, ShowProgress = True):
	global FileName
	FileName = argFileName
	File = open(FileName)

	if ShowProgress:
		progress.InitFile(File, FileName)
	while GetRec(File, OnRecord):
		if ShowProgress:
			progress.File()
	if ShowProgress:
		progress.FileDone()

OTUIds = []
SampleIds = []
OTUTable = {}

ReadRecs(ucFileName, OnRec)

fout=open(outFileName,'w')

s = "OTUId"
for SampleId in SampleIds:
	s += "\t" + SampleId

fout.write("%s\n" % s)

for OTUId in OTUIds:
	s = OTUId
	for SampleId in SampleIds:
		try:
			n = OTUTable[OTUId][SampleId]
		except:
			n = 0
		s += "\t" + str(n)
	fout.write("%s\n" % s)

fout.close()
