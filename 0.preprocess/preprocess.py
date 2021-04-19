#!/usr/bin/python
#author: Yang Yang

# importing libraries
import subprocess
import sys
import optparse
from Levenshtein import distance
import gzip
import re
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Read1 for Reads containing UMI and cell barcode.
# Read2 for genome alignment
def attach_UMI_barcode(Read1, Read2, barcodes, mismatch_rate=1):    
    mismatch_rate = int(mismatch_rate)
    f1 = open(Read1)
    f2 = open(Read2)

    line1 = f1.readline()
    line2 = f2.readline()

    while (line1):
        line1 = f1.readline()
        target = line1[0:6]
        mismatch = [distance(target, barcodes[idx]) for idx in range(len(barcodes))] 

        if (min(mismatch) <= mismatch_rate):
            barcode = barcodes[mismatch.index(min(mismatch))]
            output_file = barcode + ".align.fastq"
            f3 = open(output_file, "a")

            UMI = line1[6:12]
            first_line = "@" + barcode + "," + UMI + "," + line2[1:]
            f3.write(first_line)

            second_line = f2.readline()
            f3.write(second_line)

            third_line = f2.readline()
            third_line = "+" + barcode + "," + UMI + "," + third_line[1:]
            f3.write(third_line)

            fourth_line = f2.readline()
            f3.write(fourth_line)

            line2 = f2.readline()

        else:
            line2 = f2.readline()
            line2 = f2.readline()
            line2 = f2.readline()
            line2 = f2.readline()

        line1 = f1.readline()
        line1 = f1.readline()
        line1 = f1.readline()

    f1.close()
    f2.close()
    f3.close()

def polyAtrimmer(barcode,numberOfAs,remainedLength):
    recordIterator = FastqGeneralIterator(open(barcode + ".align.fastq", 'r'))
    numberOfProcessedReads = 0
    numberOfCountedReads = 0
    outputFile = open(barcode + ".trimmed.fastq","w")
    
    for title, seq, qual in recordIterator:
        numberOfProcessedReads += 1
        stopLocation = re.search("A{" + str(numberOfAs) + ",}",seq)
        if stopLocation is None:
            outputFile.write("@%s\n%s\n+\n%s\n" % ("NonPolyA," + title, seq, qual))
            numberOfCountedReads += 1
        else:
            if stopLocation.start() >= remainedLength:
                outputFile.write("@%s\n%s\n+\n%s\n" % ("PolyA," + title, seq[0:stopLocation.start()], qual[0:stopLocation.start()]))
                numberOfCountedReads += 1
    outputFile.close()
    summary = open(barcode + ".summary", "w")
    print("%s\t%d\t%2.2f" % (barcode,numberOfProcessedReads,numberOfCountedReads/float(numberOfProcessedReads)*100), file=summary)

if __name__ == "__main__":
    # processing options
    parser = optparse.OptionParser()
    parser.add_option("-1","--read1_fastq",dest="read1",default="None",help="the fastq file of Read1 containing cell barcodes and UMIs",type="string")
    parser.add_option("-2","--read2_fastq",dest="read2",default="None",help="the fastq file of Reads for genome alignment",type="string")
    parser.add_option("-c","--cell_barcode",dest="barcodes",default="None",help="A list of cell barcodes",type="string")
    parser.add_option("-A","--number_of_Astretch",dest="numberOfAs",default=6,help="The number of consecutive poly(A)s",type="int")
    parser.add_option("-l","--remained_length",dest="readLength",default=20,help="The remained length of reads",type="int")
    options, args = parser.parse_args()
   
    # attach UMI and cell barcode
    with open(options.barcodes) as f:
        barcodes = f.readlines()
    barcodes = [x.strip() for x in barcodes]
    #attach_UMI_barcode(options.read1, options.read2, barcodes)
    
    # trim poly(A)
    for barcode in barcodes:
        polyAtrimmer(barcode,options.numberOfAs,options.readLength)
    
    sys.exit()
