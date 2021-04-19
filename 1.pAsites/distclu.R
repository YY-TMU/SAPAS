#! /usr/bin/env Rscript

library(CAGEr)

Args = commandArgs(TRUE)
input = Args[1]
output = Args[2]

## Creating a CAGEexp object
hs = CAGEexp( genomeName = "BSgenome.Mmusculus.UCSC.mm10"
			, inputFiles = input
			, inputFilesType = "ctss"
			, sampleLabels = c("HEK293T")
)

## Reading in the data
getCTSS(hs)
CTSStagCountSE(hs)
CTSScoordinatesGR(hs)
CTSStagCountDF(hs)

## power-law distribution normalization
#normalizeTagCount(hs, method="powerLaw",fitInRange=c(5,5000),alpha=1.05,T=10^6)

## simple TPM normalization
normalizeTagCount(hs, method="simpleTpm")

## tag clustering
clusterCTSS( object = hs
           , threshold = 1
           , thresholdIsTpm = FALSE
           , nrPassThreshold = 1
           , method = "distclu"
           , maxDist = 24
           , removeSingletons = FALSE)

## get tag clusters
tc = tagClusters(hs, sample="Test")
tc = tc[,2:ncol(tc)]

## set TPM threshold
#tc = tc[which(tc$tpm > 3),]
tc = tc[which(tc$tpm.dominant_ctss > 1),]

write.table(tc,file=output,row.names=F,col.names=T,sep="\t",quote=F)

