Args = commandArgs(TRUE)
inputF1 = Args[1] # polyA usage
inputF2 = Args[2] # Gene Counts
inputF3 = Args[3] # k folds cross validation
outputF = Args[4]

library(caret)
set.seed(123)

HellingerDist = function(au) {
    ## funciton to calculate hellinger distance
    hellingerDist.p = function(x1, x2) {
        a = (sqrt(x1) - sqrt(x2))
        b = sqrt(sum(a*a)/2)
        return(b)
    }
    ## matrix formation
    k = dim(au)[2]
    ## calculate hellinger distance
    interdist = matrix(data=0,nrow=k,ncol=k)
    for(i in 1:k) {
        for(j in 1:i) {
            aui = au[,i]
            auj = au[,j]
            if(aui[1] != "NaN" & auj[1] != "NaN") {
                interdist[i,j] = hellingerDist.p(as.numeric(aui),as.numeric(auj))
            } else {
                interdist[i,j] = NaN
            }
        }
    }
    interdist = interdist + t(interdist)
    colnames(interdist) = rownames(interdist) = colnames(au)
    #non.na = !is.na(diag(interdist))
    #interdist = interdist[non.na, non.na]
    #return(stats::as.dist(interdist))
    return(interdist)
}

# neighbor voting method and calculate AUROC analytically
neighborVoting = function(setLabels, typeLabels, network) {

	# change network into rank based
	rank_dat = network
	rank_dat[] = rank(network, ties.method = "average", na.last = "keep")
	rank_dat[is.na(rank_dat)] = 0
	rank_dat = rank_dat/max(rank_dat)
	network = rank_dat

	n1 = dim(typeLabels)[2]
	n2 = dim(typeLabels)[1]
	s = unique(setLabels)

	test_typeLabels = matrix(typeLabels, nrow=n2, ncol=length(s)*n1)
	setCols = rep(s, each=n1)

	for(i in seq_along(s)) {
		d = which(setLabels == i)
		a = which(setCols == i)
		test_typeLabels[d,a] = 0
	}

	sum_in = (network %*% test_typeLabels)

	sum_all = matrix(apply(network, MARGIN=2, FUN=sum), nrow=dim(sum_in)[1], ncol=dim(sum_in)[2])

	predicts = sum_in/sum_all
	
	nans = which(test_typeLabels == 1, arr.ind=TRUE)
	predicts[nans] = NA

	for(i in seq_along(s)) {
		d = which(setLabels != i)
		a = which(setCols == i)
		predicts[d,a] = NA
	}

	predicts = apply(abs(predicts),MARGIN=2,FUN=rank,na.last="keep",ties.method="average")
	filter = matrix(typeLabels,nrow=n2,ncol=length(s)*n1)
	
	for(i in seq_along(s)) {
		d = which(setLabels != i)
		a = which(setCols == i)
		filter[d,a] = NA
	}

	negatives = which(filter == 0, arr.ind = TRUE)
	positives = which(filter == 1, arr.ind = TRUE)
	predicts[negatives] = 0
	
	np = colSums(filter, na.rm = TRUE) # Postives
	nn = apply(filter, MARGIN = 2, FUN = function(x) sum(x == 0, na.rm = TRUE))

	p = apply(predicts, MARGIN = 2, FUN = sum, na.rm = TRUE)

	rocNV = (p/np - (np+1)/2)/nn
	rocNV = matrix(rocNV, ncol = length(s), nrow = n1)

	colnames(rocNV) = s
	rownames(rocNV) = colnames(typeLabels)

	scores = rowMeans(rocNV, na.rm=TRUE)
	#return(as.numeric(scores[1]))
	return(as.vector(scores))
}


ratio = read.table(inputF1,header=T,row.names=1,sep="\t")
ratio = as.matrix(ratio)

pAs = rownames(ratio)
genes = c()
for(i in pAs) {
	gene = strsplit(i,"_")
	genes = c(genes,gene[[1]][1])
}
geneUniq = unique(genes)

# gene counts
geneCount = read.table(inputF2,header=T,row.names=1,sep="\t")

res = c("gene", "Number", "PVC.Num", "MNC.Num", "ISC.Num", "LPC.Num", "CHC.Num","CCKC.Num", "PVCauc", "MNCauc", "ISCauc", "LPCauc", "CHCauc", "CCKCauc", "auc")

write(paste(res,collapse="\t"),file=outputF,append=T)

# get k fold cv
k = as.numeric(inputF3)

for(gene in geneUniq) {
	idx = which(genes == gene)
	tmpCount = as.vector(as.matrix(geneCount[gene,]))
	idy = which(tmpCount>=5)
	tmpExpNum = length(idy) 
	tmpRatio = ratio[idx,idy]

	if(tmpExpNum > 2) {
		hdistMx = HellingerDist(tmpRatio)
		non.na = !is.na(diag(hdistMx))
		hdistM = hdistMx[non.na, non.na]
		hdist = hdistM	
		# convert distance to similarity
		#hsim = max(hdist) - hdist
		hsim = 1 - hdist	
	
		#hdistMxAll = hdistMxAll + hdistMx

		#cellsTsne = Rtsne(hdist, is_distance=TRUE)
		#rownames(cellsTsne$Y) = colnames(tmpRatio)
		#colnames(cellsTsne$Y) = c("tsne1", "tsne2")

		# expressed cell numbers
		Num = tmpExpNum

		anno = c()
		for(i in colnames(tmpRatio)) {
			anno = c(anno, strsplit(i,"\\.")[[1]][1])
		}

		anno[anno=="CHC1"] = "CHC"
		anno[anno=="CHC2"] = "CHC"

		## combine
		typeLabels = matrix(rep(0,Num*6),ncol=6)
		colnames(typeLabels) = c("PVC","MNC","ISC","LPC","CHC","CCKC")
		rownames(typeLabels) = anno
		typeLabels[which(anno == "PVC"),1] = 1
		typeLabels[which(anno == "MNC"),2] = 1
		typeLabels[which(anno == "ISC"),3] = 1
		typeLabels[which(anno == "LPC"),4] = 1
		typeLabels[which(anno == "CHC"),5] = 1
		typeLabels[which(anno == "CCKC"),6] = 1

		typeNums = as.vector(colSums(typeLabels))

		setLabels = rep(0,length(anno))
		kfolds = createFolds(anno,k)
		for(i in 1:length(kfolds)) {
			setLabels[kfolds[[i]]] = i
		}
		
		auc = neighborVoting(setLabels,typeLabels,hsim)
		
		res = c(gene,Num, typeNums, auc, mean(auc,na.rm=TRUE))

		write(paste(res,collapse="\t"),file=outputF,append=T)
	}
}


