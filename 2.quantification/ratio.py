import os
import sys
import pandas as pd
import math

def GeneCounts(infile2):
	in2 = open(infile2)
	headline = in2.readline()
	GeneDict = {}
	for line in in2:
		fields = line.strip().split("\t")
		gene = fields[0]
		GeneDict[gene] = [int(x) for x in fields[1:]]
	return GeneDict

def Counts(infile1,outfile,GeneDict):
	in1 = open(infile1)
	out = open(outfile,"w")
	headline = in1.readline()
	out.write(headline)
	for line in in1:
		fields = line.strip().split("\t")
		APA = fields[0]
		gene = fields[0].split("_")[0]
		counts = [float(x) for x in fields[1:]]
		total = GeneDict[gene]
		ratio = [str(counts[i]/total[i]) if total[i]>0 else str(float("nan")) for i in range(len(counts))]
		ratio.insert(0,APA)
		out.write("\t".join(ratio) + "\n")

if __name__ == "__main__":
	infile1 = sys.argv[1] # isoform counts
	infile2 = sys.argv[2] # gene counts
	outfile = sys.argv[3]
	GeneDict = GeneCounts(infile2)	
	Counts(infile1,outfile,GeneDict)

