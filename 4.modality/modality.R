Args = commandArgs(TRUE)
inputF1 = Args[1] # ratio
outputF = Args[2]

library(textmineR)

## desired distribution
distal = c(0,0,1)
proximal = c(1,0,0)
middle = c(0.1,0.8,0.1)
bimodal = c(0.45,0.1,0.45)
multimodal = c(1/3,1/3,1/3)

## binning function 
binning = function(x) {
	binsHist = hist(x,breaks=c(0,1/3,2/3,1))$counts
	bins = binsHist/sum(binsHist)
	return(bins)
}

## Jensen-Shannon Divergence
## JSD calculation for each desired distribtuion
JSDmin = function(x) {
	JSD_distal = CalcJSDivergence(x, distal)
	JSD_proximal = CalcJSDivergence(x, proximal)
	JSD_middle = CalcJSDivergence(x, middle)
	JSD_bimodal = CalcJSDivergence(x, bimodal)
	JSD_multimodal = CalcJSDivergence(x, multimodal)
	
	types = c("distal", "proximal", "middle", "bimodal", "multimodal")
	JSDs = c(JSD_distal,JSD_proximal,JSD_middle,JSD_bimodal,JSD_multimodal) 

	type = types[which(JSDs==min(JSDs))]
	res = c(JSDs,type)
}

## main
# isoform ratio
ratio = read.table(inputF1,header=T,row.names=1,sep="\t")
pAs = rownames(ratio)

for(i in pAs) {
	tmpRatio = na.omit(as.numeric(ratio[i,]))
	tmpNums = length(tmpRatio)
	
	tmpRatioBinned = binning(tmpRatio)
	tmpRes = JSDmin(tmpRatioBinned)

	tmpRes = c(i, tmpRes, tmpNums)
	write(paste(tmpRes,collapse="\t"),file=paste(outputF,".binned.txt",sep=""),append=T)
}

