#! /usr/bin/env python

import os
import sys
import pybedtools
from Bio import SeqIO
from Bio.Seq import Seq

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

def filter(candidate_bed,output_dir):
    candidate_bed = pybedtools.BedTool(candidate_bed)

    ## intersect with protein coding genes
    coding_bed = candidate_bed.intersect("annotation/protein_coding.genes.bed", s=True, u=True)
    coding_bed.saveas("coding.bed.tmp")

    ## remove internal priming
    output = open(output_dir + "/result_coding.polyAsites.bed","w")
    for line in open("coding.bed.tmp","r"):
	fields = line.strip("\n").split("\t")
	chrom = fields[0]
	start = fields[1]
	end = fields[2]
	strand = fields[5]

	stretchType = internal_polyA_priming(chrom, int(start), int(end), strand)

	if stretchType == "None":
	    pAnames = chrom + ":" + start + "-" + end + ":" + strand
	    output.write(chrom + "\t" + start + "\t" + end + "\t" + pAnames + "\t.\t" + strand + "\n")
    output.close()
    os.system("rm coding.bed.tmp")

    ## annotate using GENCODE annotations
    result_coding_pA = pybedtools.BedTool(output_dir + "/result_coding.polyAsites.bed")	
    # exonic
    result_exonic_pA = result_coding_pA.intersect("annotation/exons.bed", s=True, u=True)
    # intronic
    result_intronic_pA = result_coding_pA.intersect("annotation/exons.bed", s=True, v=True)
    # 3UTRs
    result_3UTR_pA = result_coding_pA.intersect("annotation/3UTRs.bed", s=True, u=True)
    # store
    result_exonic_pA.saveas(output_dir + "/result_exonic.polyAsites.bed")
    result_intronic_pA.saveas(output_dir + "/result_intronic.polyAsites.bed")
    result_3UTR_pA.saveas(output_dir + "/result_3UTR.polyAsites.bed")

if __name__ == "__main__":
    # processing options
    parser = optparse.OptionParser()
    parser.add_option("-i","--input_bed",dest="inputBed",default="None",help="The bed used to filter",type="string")
    parser.add_option("-g","--genome",dest="genome",default="None",help="The genome sequence in fasta format",type="string")
    parser.add_option("-o","--output_dir",dest="outputDir",default="None",help="The directory used to store results",type="string")
    options, args = parser.parse_args()

    genome = SeqIO.to_dict(SeqIO.parse(open(options.genome,"r"),"fasta"))
    filter(options.inputBed,options.outputDir)

