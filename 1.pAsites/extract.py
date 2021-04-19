#!/usr/bin/python
#author: Yang Yang

import os
import sys
import optparse
from Bio import SeqIO
from Bio.Seq import Seq

def extract(inBed, out3end):
    reads = open(inBed, "r")
    ends  = open(out3end, "w")
    for read in reads:
        list = read.strip().split("\t")
        ReadName = list[3]
        ReadChr = list[0]
        ReadStart = int(list[1])
        ReadEnd = int(list[2])
        ReadStrand = list[5]
        polyA = ReadName.strip().split(",")[0]
        if polyA == "PolyA":
            if ReadStrand == "+":
                tmpStart = str(ReadEnd - 1)
                tmpEnd = str(ReadEnd)
                list[1] = tmpStart
                list[2] = tmpEnd
                ends.write("\t".join(list) + "\n")
            elif ReadStrand == "-":
                tmpStart = str(ReadStart)
                tmpEnd = str(ReadStart + 1)
                list[1] = tmpStart
                list[2] = tmpEnd
                ends.write("\t".join(list) + "\n")
            else:
                pass

def internal_polyA_priming(chrom, start, end, strand):
    STRETCH_NUMs = 6 #### consecutive A numbers
    Astretch = "A"*STRETCH_NUMs

    seq = str(genome[chrom].seq)[int(start)-20:int(end)+20]
    sequence = seq if strand == "+" else str(Seq(seq).reverse_complement())

    max_index = len(sequence) - STRETCH_NUMs +1
    windowSeqs = []
    for index in range(max_index):
        windowSeqs.append(sequence[index:(index+STRETCH_NUMs)].upper())

    stretchType = "Astretch" if Astretch in windowSeqs else "None"

    return stretchType

if __name__ == "__main__":
    # processing options
    parser = optparse.OptionParser()
    parser.add_option("-b","--input_bam",dest="inputBam",default="None",help="The bam used to extract polyA reads",type="string")
    parser.add_option("-g","--genome",dest="genome",default="None",help="The genome sequence in fasta format",type="string")
    parser.add_option("-o","--output",dest="output",default="None",help="The output bed file for polyA reads' 3ends",type="string")
    options, args = parser.parse_args()

    # Convert bam into bed using system bedtools
    os.system("bamToBed -i " + options.inputBam + " -cigar > " + options.inputBam + ".bed")
    # extract polyA reads 3'end
    extract(options.inputBam + ".bed", options.inputBam + ".3end")
    # remove internal primmed reads
    genome = SeqIO.to_dict(SeqIO.parse(open(options.genome,"r"),'fasta'))
    raw3ends = open(options.inputBam + ".3end","r")
    new3ends = open(options.output,"w")
    for line in raw3ends:
        fields = line.strip("\n").split("\t")
        chrom = fields[0]
        start = fields[1]
        end = fields[2]
        strand = fields[5]
        stretchType = internal_polyA_priming(chrom, int(start), int(end), strand)

        if stretchType == "None":
            new3ends.write(line)
    
    os.system("rm " + options.inputBam + ".bed")
    os.system("rm " + options.inputBam + ".3end")
    sys.exit()
	
