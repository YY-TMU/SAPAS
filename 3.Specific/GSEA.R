Args = commandArgs(TRUE)
inputF1 = Args[1]
inputF2 = Args[2]
outputF = Args[3]

library(fgsea)
library(reactome.db)
library(data.table)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

## read ranks files
rnk.file = read.table(inputF1,header=T,row.names=1,sep="\t")
GeneIDs = rownames(rnk.file)
ranks = setNames(rnk.file[,1], GeneIDs)

## read gmt file
gmt.file = system.file("extdata", inputF2, package="fgsea")

## loading gmt file
pathways = gmtPathways(gmt.file)

## convert Entrez ID to Ensembl ID
idmap = read.table("04.AUC/mart_export.txt",header=T,sep="\t")
Bimap = as.list(org.Hs.egENSEMBL)
pathways_ens = list()

for(i in names(pathways)) {
	geneList = pathways[[i]]
	geneList_ens = as.vector(unlist(Bimap[geneList]))
	geneList_ens = geneList_ens[!is.na(geneList_ens)]
	geneList_ens = unique(as.vector(idmap[which(idmap[,1] %in% geneList_ens),2]))
	geneList_ens = intersect(geneList_ens, GeneIDs)
	pathways_ens[[i]] = geneList_ens
}
	
fgseaRes = fgsea(pathways = pathways_ens,
                 stats = ranks,
                 minSize = 15,
				 maxSize = 500,
                 nperm = 10000)

fgseaRes[, leadingEdge := lapply(leadingEdge, mapIds, x=org.Mm.eg.db, keytype="ENSEMBL", column="SYMBOL")]

res = fgseaRes[fgseaRes[,pval<0.05],]

save(ranks, pathways_ens, fgseaRes, file=paste(outputF,".fgsea.Rdata",sep=""))
fwrite(res, file=paste(outputF,".res",sep=""), sep="\t", sep2=c("", " ",""))
