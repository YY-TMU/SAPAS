#! /usr/bin/env Rscript

library(optparse)
library(CAGEr)

parser = OptionParser()
parser = add_option(parser, c("-i","--input"), dest="input", default="NA", help="The input file")
parser = add_option(parser, c("-o","--output"), dest="output", default="NA", help="The output file")
parser = add_option(parser, c("-g","--genome"), dest="genome", default="BSgenome.Mmusculus.UCSC.mm10", help="Specify a genome assembly that i readily available through BSgenome package")
parser = add_option(parser, c("-s","--sample"), dest="sample", default="NA", help="sample name to record")
opt = parse_args(parser)

## Creating a CAGEexp object
hs = CAGEexp(genomeName=opt$genome, inputFiles=opt$input, inputFilesType="ctss", sampleLabels= c(opt$sample))

## Reading in the data
getCTSS(hs)
CTSStagCountSE(hs)
CTSScoordinatesGR(hs)
CTSStagCountDF(hs)

## simple TPM normalization
normalizeTagCount(hs, method="simpleTpm")

## tag clustering
clusterCTSS( object = hs, 
	     threshold = 1, 
	     thresholdIsTpm = FALSE, 
	     nrPassThreshold = 1, 
	     method = "distclu", 
	     maxDist = 24, 
	     removeSingletons = FALSE)

## get tag clusters
tc = tagClusters(hs, sample=opt$sample)
tc = tc[,2:ncol(tc)]

## set TPM threshold
tc = tc[which(tc$tpm.dominant_ctss > 1),]

## write in bed format
for(i in 1:nrow(tc)){
    write(paste(tc[i,1],tc[i,7]-1,tc[i,7],".",".",tc[i,4],sep="\t"),file=opt$output,append=T)
}




