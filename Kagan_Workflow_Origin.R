##Load lumi
library(lumi)
setwd("C:/Users/Courtney/Dropbox/LCL-iPSC/Origin/")

##Create object with name of data file:
data = c('YGilad-CK-Aug23-14-ProbeLevelData-NotNormalized-NoBGSubtracted-FinalReport.txt')

##Extract raw data from lumi file and preparing for removing bad probes:
data.lumi = lumiR.batch(data, lib.mapping=NULL, convertNuID=F,annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION','PROBE_ID'))

###Remove Probes
all.probes = data.lumi@featureData[[5]]
goodprobes= read.table('ht12_probes_snps_ceu_hg19_af_0.05_map_37.txt', header=T, sep='\t')
#23,890 probes
probes = goodprobes$probeID
## Convert from factor to character
probes = as.character(goodprobes$probeID)
cleanprobes = which(all.probes %in% probes)
data.lumi.clean = data.lumi[cleanprobes,]
#23,888 probes
head(data.lumi@featureData[[5]])
head(data.lumi.clean@featureData[[5]])

#Add in sample names
samplenames = read.table('covar.txt', header=T, sep ='\t')
#Re-order samplenames based on array location
samplenames = samplenames[order(samplenames$Order),]
sampleNames(data.lumi.clean) = samplenames$NewName


#Get QC Data
summary(data.lumi.clean, 'QC')

#Past output into excel and create a file called lumiQC
#Sample  mean standard.deviation detection.rate.0.01. distance.to.sample.mean
#1    YG1 7.229              1.211               0.3735                   42.24
#2    YG2 7.342              1.257               0.3543                   46.28

#Look at QC based on cell type
qcdata = read.table('lumiQC.txt', header=T, as.is=T, sep='\t')


boxplot(qcdata$mean~samplenames$Indiv, main = 'Mean Probe Intensity by Individual')
boxplot(qcdata$mean~samplenames$Batch, main = 'Mean Probe Intensity by Array')
boxplot(qcdata$mean~samplenames$Sex, main = 'Mean Probe Intensity by Sex')
boxplot(qcdata$mean~samplenames$Type, main = 'Mean Probe Intensity by Cell Type')

### NORMALIZATION: log2 stabilized and quantile normalization ###
data.norm.all <- lumiExpresso(data.lumi.clean, bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)
sampleNames(data.norm.all) = samplenames$NewName
show(data.norm.all)

##Normalize by cell type
#ipscsubset = grep ("LCL|Fib", samplenames$Name, invert=T)
#lclsubset = grep ("LCL", samplenames$Name)
#fibsubset = grep ("Fib", samplenames$Name)

#If you want to normalize by cell type - do not do for overall analysis
#data.norm.ipsc = lumiExpresso(data.lumi[, ipscsubset], bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)
#summary(data.norm.ipsc, 'QC')
#data.norm.fib = lumiExpresso(data.lumi[, fibsubset], bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)
#data.norm.lcl = lumiExpresso(data.lumi[, lclsubset], bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)
#data.norm.comb = combine(data.norm.ipsc,data.norm.fib,data.norm.lcl)
#sampleorder = c(as.character(samplenames[ipscsubset,13]),as.character(samplenames[fibsubset,13]),as.character(samplenames[lclsubset,13]))
#colnames(data.norm.comb) = sampleorder

#QC
summary(data.norm.all, 'QC')

#Summary of Samples:
#  YG1     YG2     YG3     YG4    YG5     YG6     YG7      YG8     YG9     YG10     YG11    YG12
#mean                     7.2660  7.2660  7.2660   7.266  7.266  7.2660  7.2660   7.2660  7.2660   7.2660   7.2660  7.2660
#standard deviation       1.2430  1.2430  1.2430   1.243  1.243  1.2430  1.2430   1.2430  1.2430   1.2430   1.2430  1.2430
#detection rate(0.01)     0.3735  0.3543  0.3672   0.340  0.357  0.3681  0.3657   0.3482  0.3732   0.3743   0.3535  0.3839
#distance to sample mean 43.1100 43.2500 44.4300 108.900 44.130 51.5700 43.9500 111.1000 42.8200 114.2000 108.1000 42.9800
#YG13     YG14    YG15     YG16     YG17    YG18    YG19   YG20    YG21    YG22    YG23    YG24
#mean                      7.2660   7.2660  7.2660   7.2660   7.2660  7.2660  7.2660  7.266  7.2660  7.2660  7.2660  7.2660
#standard deviation        1.2430   1.2430  1.2430   1.2430   1.2430  1.2430  1.2430  1.243  1.2430  1.2430  1.2430  1.2430
#detection rate(0.01)      0.3619   0.3296  0.3807   0.3292   0.3382  0.3848  0.3582  0.363  0.3656  0.3587  0.3666  0.3534
#distance to sample mean 112.9000 107.9000 42.4700 109.7000 111.8000 42.7200 45.2900 48.500 42.9400 42.6400 45.7400 44.4300

qcdata.n = read.table('lumiQC_norm.txt', header=T, as.is=T, sep='\t')

boxplot(qcdata.n$mean~samplenames$Indiv, main = 'Mean Probe Intensity by Individual')
boxplot(qcdata.n$mean~samplenames$Batch, main = 'Mean Probe Intensity by Array')
boxplot(qcdata.n$mean~samplenames$Sex, main = 'Mean Probe Intensity by Sex')
boxplot(qcdata.n$mean~samplenames$Type, main = 'Mean Probe Intensity by Cell Type')

#Distance to sample mean - find new qc summary
#boxplot(data.norm.all@assayData$exprs~samplenames$Indiv, main = 'Mean Probe Intensity by Individual')
#boxplot(qcdata.n$distance.to.sample.mean~samplenames$Batch, main = 'Mean Probe Intensity by Array')
#boxplot(qcdata.n$distance.to.sample.mean~samplenames$Sex, main = 'Mean Probe Intensity by Sex')
#boxplot(qcdata.n$distance.to.sample.mean~samplenames$Type, main = 'Distance to Sample Mean by Cell Type')


##Look at plots of array data (boxplot, density, etc) :
pdf(file = "QC_Normalization.pdf")
boxplot(data.lumi.clean, main= "Pre-normalization")
#plot(data.lumi.clean, what='density', main= "Density plot of intensity")
#plot(data.norm.all, what='density', main = "Density plot of intensity - Normalized ")
#plot(data.norm.comb, what='density')
boxplot(data.norm.all, main = "Post-Normalization")
dev.off()
#plot(data.norm.comb, what = 'boxplot', main = "Normalized by cell type")
plot(data.norm.all, what='density')

##Check that replicates are most related
plot(data.norm.all, what='sampleRelation')


#Subset data by cell type and look for detection in at least 4 indivduals per subset
detection = data.norm.all@assayData$detection
#colnames(detection) = sampleorder
colnames(detection) = samplenames$NewName

#iPSCdetect = detection[,type %in% selected]
#iPSCdetect = detection[,grep ("LCL|Fib", colnames(detection), invert=T)]
LiPSCdetect = detection[,grep ("L-iPSC", colnames(detection))]
Ind1detect = LiPSCdetect[,grep ("Ind1", colnames(LiPSCdetect))]
Ind2detect = LiPSCdetect[,grep ("Ind2", colnames(LiPSCdetect))]
Ind3detect = LiPSCdetect[,grep ("Ind3", colnames(LiPSCdetect))]
Ind4detect = LiPSCdetect[,grep ("Ind4", colnames(LiPSCdetect))]
FiPSCdetect = detection[,grep ("F-iPSC", colnames(detection))]
LCLdetect = detection[,grep ("LCL", colnames(detection))]
Fibdetect = detection[,grep ("Fib", colnames(detection))]

#Count the number of probes that have a detection p-value<.05
#LiPSCdetected = rowSums(LiPSCdetect<0.05)
Ind1detected = rowSums(Ind1detect <0.05)
Ind2detected = rowSums(Ind2detect <0.05)
Ind3detected = rowSums(Ind3detect <0.05)
Ind4detected = rowSums(Ind4detect <0.05)
FiPSCdetected = rowSums(FiPSCdetect<0.05)
LCLdetected = rowSums(LCLdetect<0.05)
Fibdetected = rowSums(Fibdetect<0.05)

#Here is where I select the number of indiv that need to have the probe expressed (at least in 3 people), iPSC =22,848, LCL = 15,759, Fib = 16,699
#detect.LiPSC <- which(LiPSCdetected > 2)
detect.Ind1 = which(Ind1detected > 0)
detect.Ind2 = which(Ind2detected > 0)
detect.Ind3 = which(Ind3detected > 0)
detect.Ind4 = which(Ind4detected > 0)
detect.LCL <- which(LCLdetected > 2)
detect.Fib <- which(Fibdetected > 2)
detect.FiPSC <- which(FiPSCdetected > 2)

indiv_union = c(detect.Ind1, detect.Ind2, detect.Ind3, detect.Ind4)
indiv_union_unique = unique(indiv_union)

indiv_freq = as.data.frame(table(indiv_union))
Indiv_detected = indiv_freq[which(indiv_freq[,2] >2 ),]
table(indiv_freq$Freq)

#1     2     3     4 
#1406   879   836 12773 
#write.table(detect.LiPSC, 'ProbesDetected_byIndiv_L-iPSC_test2.txt', quote=F, row.names =F, sep='\t')

detect.LiPSC = as.matrix(Indiv_detected$indiv_union)
detect.LiPSC = as.integer(detect.LiPSC)

union = c(detect.Fib, detect.LCL, detect.FiPSC, detect.LiPSC)
unionunique= unique(union)
detect.ind.all = sort.int(unionunique)
#15,293 probes detected 

norm_quant.all <- data.norm.all@assayData$exprs
###Find the column that is lumi_ID in feature data usually this column
head(data.norm.all@featureData[[5]])
##[1] "ILMN_1343291" "ILMN_1343295" "ILMN_1651199" "ILMN_1651209" "ILMN_1651210" "ILMN_1651221"
###Convert expr_quant rownames to
rownames(norm_quant.all)=data.norm.all@featureData[[5]]
#With this threshold 15,293 probes out of 47,315 are expressed
expr_quant.all <- norm_quant.all[detect.ind.all,]
#expr_iPSC_quant.all = norm_quant.all[detect.iPSC,]
#expr_quant.all = expr_iPSC_quant.all


##Load in Covariates
#Converted some factors so they are categorical 
colnames(expr_quant.all) = samplenames$NewName
#colnames(expr_quant.all) = sampleorder
#samplenames = read.table('covar_reordered.txt', header=T, sep ='\t')

batch = samplenames$Batch
mef = samplenames$MEF
type = samplenames$Type
sex =  samplenames$Sex
indiv = samplenames$Indiv
pluri =  samplenames$Pluri
novel = samplenames$Novelty
der = samplenames$Deriv

covars = list(mef, batch, type,sex,indiv,pluri,novel,der)

#Converted categorical covariates to a factor so they are levels 
batch.f = as.factor(batch)
mef.f = as.factor(mef)
type.f = as.factor(type)
sex.f = as.factor(sex)
indiv.f = as.factor(indiv)

##To look for correlations between covariates get the R^2 of the Pearson correlation
cor(novel, pluri)

#Converted numerical covariates to a numeric so they are continuous
pluri.num =as.numeric(pluri)
novel.num = as.numeric(novel)

#rm(cond, spec, tech, sex, age, pluri, novel)
#Looking at just probes expressed in iPSCs
#expr_quant.all = expr_iPSC_quant.all.clean

## Finding the unique gene names matching probes to gene names using the good probe list
gene_names=c()
for(i in 1:dim(expr_quant.all)[1]){
  gene_names=c(gene_names,as.vector(goodprobes[as.vector(goodprobes[,4])==row.names(expr_quant.all)[i],11]))
}

symbolsUniq = unique(gene_names)
length(symbolsUniq)
#[1] 11,469 - that have HGNC annotation
# 11,713 have ENSName


## This loop will give the most 3' value for multiple probes within the same gene. In the end you get a simple file with all genes that are expressed with the corresponding mean intensity expression levels across its different probes.
expr_gene = matrix(NA, ncol=24, nrow=length(unique(gene_names)))
i=0
for(gene in unique(gene_names)){
  i = i+1
  
  currRows = which(gene_names == gene)
  if(length(currRows)>1){
    if(goodprobes[currRows[1],6]=="+"){
      keepRow = currRows[which.max(goodprobes[currRows,2])]
    }
    else{
      keepRow = currRows[which.min(goodprobes[currRows,2])]
    }
  }
  else{
    keepRow=currRows[1]
  }
  expr_gene[i,] = expr_quant.all[keepRow,]
  
} 
dim(expr_gene)
rownames(expr_gene) = unique(gene_names)
colnames(expr_gene) = colnames(expr_quant.all)
#colnames(expr_gene) =samplenames$Name

write.table(expr_gene, 'OriginGeneExpression_Normalized.txt', sep='\t', row.names=T, quote=F)
#expr_gene = read.table('OriginGeneExpression_Normalized.txt', header=T, as.is=T, sep='\t', row.names=1)
#samplenames = read.table('covar.txt', header=T, sep ='\t')
#Re-order samplenames based on array location
#samplenames = samplenames[order(samplenames$Order),]
#colnames(expr_gene) = samplenames$NewName

####Data Analysis####

#To use Nick's dendogram script
#First generate the chr file
goodprobes= read.table('ht12_probes_snps_ceu_hg19_af_0.05_map_37.txt', header=T, sep='\t')
probeinfolist = cbind(goodprobes$chr, as.character(goodprobes[,4]),as.character(goodprobes[,11]))
probelist = rownames(expr_quant.all)
genelist = rownames(expr_gene)
probelist = as.matrix(probelist)
genelist = as.matrix(genelist)
colnames(probelist) = c("ILMN")
colnames(genelist) = c("GeneID")
colnames(probeinfolist) = c("Chr", "Probe", "Gene")
geneinfolist = probeinfolist[,-2]
test = unique(geneinfolist)
chrlist.gene = merge(genelist, test, by.x = "GeneID", by.y = "Gene", all.x = T, all.y = F, sort=F)
chrlist.probe = merge(probelist, probeinfolist, by.x = "ILMN", by.y = "Probe", all.x = T, all.y = F, sort=F)
chrlist.gene = merge(genelist, geneinfolist, by.x = "GeneID", by.y = "Gene", all.x = T, all.y = F, sort=F)
write.table(chrlist.gene, 'GeneList.txt', sep='\t', row.names=T, quote=F)
write.table(genelist, 'ActualGeneList.txt', sep='\t', row.names=T, quote=F)

cleangenes = which(chrlist.gene$GeneID %in% genelist)
clean.chrlist.gene = chrlist.gene[cleangenes,]
chrfinal.p = as.matrix(chrlist.probe[,2])
chrfinal.g = as.matrix(clean.chrlist.gene[,2])
xchr = chrlist.gene[which(chrlist.gene$Chr ==23 ),]
#431 X chr genes
ychr = chrlist.gene[which(chrlist.gene$Chr ==24 ),]
#7 Y chr genes
hist(as.numeric(chrfinal.g), breaks = 24, xlim = c(1,24))
sexgene = c(which(chrlist.gene$Chr ==23 ), which(chrlist.gene$Chr ==24 ))
expr_gene_nosex = expr_gene[-sexgene,]
# 11,277 genes


#If you want mean subtracted and variance divided data
#stan = apply(expr_quant.all,1, function(x) x-mean(x))
#stand = t(stan)

#Variance by probe
variance.probe.LCL = apply(expr_quant.all[,grep ("LCL", colnames(expr_quant.all))],1,var)
variance.probe.Fib = apply(expr_quant.all[,grep ("Fib", colnames(expr_quant.all))],1,var)
head(expr_quant.all[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)])
variance.probe.iPSC = apply(expr_quant.all[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)],1,var)

#variance by gene
variance.LCL = apply(expr_gene[,grep ("LCL", colnames(expr_gene))],1,var)
variance.Fib = apply(expr_gene[,grep ("Fib", colnames(expr_gene))],1,var)
head(expr_gene[,grep ("LCL|Fib", colnames(expr_gene),invert=T)])
variance.iPSC = apply(expr_gene[,grep ("LCL|Fib", colnames(expr_gene),invert=T)],1,var)
mean(variance.LCL)
# 0.03973149
mean(variance.iPSC)
# 0.02424886
mean(variance.Fib)
# 0.03112353

variance.all = cbind(variance.LCL, variance.Fib, variance.iPSC)
boxplot(variance.all)#, ylim = c(-.01,.15))
write.table(variance.all, 'Variance by cell type.txt', quote=F, sep='\t')
var.test(variance.LCL, variance.iPSC)
var.test(variance.Fib, variance.iPSC)
var.test(variance.Fib, variance.LCL)

#F test to compare two variances

#data:  variance.Fib and variance.iPSC
#F = 3.2088, num df = 11712, denom df = 11712, p-value < 2.2e-16
#alternative hypothesis: true ratio of variances is not equal to 1
#95 percent confidence interval:
#  3.094661 3.327177
#sample estimates:
#  ratio of variances 
#3.208813 

library(ggplot2)
var_all <- data.frame(var=c(variance.LCL, variance.Fib, variance.iPSC), type = rep(c("LCL","Fib", "iPSC"), times=c(length(variance.LCL))))
ggplot(var_all, aes(x=var, fill=type)) + geom_density(alpha=0.01) +xlim(-.01,.2)+xlab("Variance") + ggtitle("Gene Expression: Total variance") + theme(legend.position=c(.75,.75)) +theme(text = element_text(size=23))

var_LCL <- data.frame(var=c(variance.LCL), type = rep(c("LCL"), times=c(length(variance.LCL))))
ggplot(var_LCL, aes(x=var, fill=type)) + geom_density(alpha=0.5) +xlim(-.01,.2)+xlab("Variance") + ggtitle("Gene Expression: Total variance") + theme(legend.position=c(.75,.75)) +theme(text = element_text(size=23))
var_iPSC <- data.frame(var=c(variance.iPSC), type = rep(c("iPSC"), times=c(length(variance.iPSC))))
ggplot(var_iPSC,aes(x=var, fill=type))+ geom_density(alpha=0.5) +xlim(-.01,.2)+xlab("Variance") + ggtitle("Gene Expression: Total variance") + theme(legend.position=c(.75,.75)) +theme(text = element_text(size=23))

#Read in data
avg_beta = expr_gene
#avg_beta = expr_quant.all
#avg_beta = stand
#chr = chrfinal
chr = chrfinal.g

#Open pdf

pdf(file = "Dendrograms.pdf")

pdf(file = "Dendrograms_nosex_1k.pdf")
variance.iPSC = apply(expr_gene_nosex[,grep ("LCL|Fib", colnames(expr_gene_nosex),invert=T)],1,var)
varall = cbind(expr_gene_nosex, variance.iPSC)
varall_1k = varall[order(varall[,25], decreasing = T),]
var1k = varall_1k[1:1000,-25]


avg_beta = var1k
# Make dendogram of the data
labs = c("Ind1 F-iPSC", 
         "Ind2 L-iPSC C", "Ind3 L-iPSC B", "Ind1 LCL", "Ind4 L-iPSC C", 
         "Ind2 L-iPSC B", "Ind4 L-iPSC A", "Ind1 Fibroblast", "Ind1 L-iPSC A", 
         "Ind4 Fibroblast", "Ind2 LCL", "Ind3 F-iPSC", "Ind2 Fibroblast", 
         "Ind3 LCL", "Ind2 L-iPSC A", "Ind4 LCL", "Ind3 Fibroblast", "Ind1 L-iPSC C", 
         "Ind2 F-iPSC", "Ind3 L-iPSC A", "Ind1 L-iPSC B", "Ind4 L-iPSC B", 
         "Ind3 L-iPSC C", "Ind4 F-iPSC")
cols = c("Navy",
         "darkorange1",
         "Black",
         "ForestGreen",
         "lightcoral",
         "darkorange1",
         "lightcoral",
         "Brown",
         "Navy",
         "Brown",
         "ForestGreen",
         "Black",
         "Brown",
         "ForestGreen",
         "darkorange1",
         "ForestGreen",
         "Brown",
         "Navy",
         "darkorange1",
         "Black",
         "Navy",
         "lightcoral",
         "Black",
         "lightcoral")
#install.packages("ClassDiscovery", repos="http://R-Forge.R-project.org")
library(ClassDiscovery)
euc = hclust(dist(t(avg_beta[,1:24])))
plotColoredClusters(euc, labs, cols, cex = 1,lwd = 3, lty = 1,main = "Euclidean Distance", line = -1, xlab="", sub="", col.main = "#45ADA8")

#plot(hclust(dist(t(avg_beta[,1:24]))), xlab = "hclust/Euclidean distance", main = "Dendrogram")
#plot(hclust(dist(t(avg_beta[,1:24]), method = "manhattan")), xlab = "Manhattan", main = "Dendrogram")
#plot(hclust(dist(t(avg_beta[,1:24]), method = "canberra")), xlab = "Canberra", main = "Dendrogram")

#Make dendograms using pearson
pears = hclust(as.dist(1-cor(as.matrix(avg_beta))))
plotColoredClusters(pears, labs, cols, cex = 1,lwd = 3, lty = 1,main = "Pearson Correlation", line = -1, xlab="", sub="", col.main = "#45ADA8")
#plot(hclust(as.dist(1-cor(as.matrix(avg_beta)))),xlab = "Pearson", main = "Dendrogram using Pearson Correlation")

# Make dendogram of the data without the X chr
#plot(hclust(as.dist(1-cor(as.matrix(avg_beta[chr[,1] != "23" ,])))), xlab = "hclust/Euclidean distance", main = "Dendrogram w/o X chr")
pdf(file = "Dendrograms_nosex.pdf")
avg_beta = expr_gene_nosex
euc_nox = hclust(dist(t(avg_beta)))
plotColoredClusters(euc_nox, labs, cols, cex = 1,lwd = 3, lty = 1,main = "Euclidean Distance - No Sex Chromosomes", line = -1, xlab="", sub="", col.main = "#45ADA8")
pears_nox = hclust(as.dist(1-cor(as.matrix(avg_beta))))
plotColoredClusters(pears_nox, labs, cols, cex = 1,lwd = 3, lty = 1,main = "Peasron Correlation - No Sex Chromosomes", line = -1, xlab="", sub="", col.main = "#45ADA8")
dev.off()

# Make dendogram of iPSCs only
plot(hclust(dist(t(avg_beta[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]))), xlab = "", main = "Dendrogram with only iPSCs")

#Pearson
plot(hclust(as.dist(1-cor(avg_beta[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]))), xlab = "", main = "Dendrogram with only iPSCs")

# Make dendogram of iPSCs only w/o X chr
plot(hclust(as.dist(1-cor(as.matrix(avg_beta[chr[,1] != "23" ,grep ("LCL|Fib" , colnames(avg_beta), invert = T)][,1:16])))), xlab = "", main = "Dendrogram with only iPSCs, w/o X chr")


# Make heatmap of the data
pdf(file="Heatmaps.pdf")
library(gplots)
heatmap.2(cor(as.matrix(expr_gene[,1:24]), use = "complete"), margins=c(7,3),trace="none",main="Gene Expression Correlation", key=T, keysize=1.5,density.info="none",cexCol=0.9, labRow=NA)
heatmap.2(cor(as.matrix(expr_gene_nosex[,grep ("LCL|Fib" , colnames(expr_gene_nosex), invert = T)]), use = "complete"), margins=c(7,3),trace="none",main="Gene Expression Correlation \niPSCs Without X", key=T, keysize=1.5,density.info="none",cexCol=0.9, labRow=NA)

dev.off()

pdf(file = "Xchr_expression.pdf")

# Make histogram of all of the X chromosome 

for (i in 1:24){
  hist(avg_beta[chr == "24", i], main = colnames(avg_beta)[i])
}

dev.off()

#Test for robustness
pdf(file = "dendrogram_robustness.pdf")
plot(hclust(dist(t(avg_beta[chr != "24" ,grep ("LCL|Fib" , colnames(avg_beta), invert = T)][,1:16]))), xlab = "", main = "Dendrogram with only iPSCs, w/o X chr")
for (i in 1:20){
  tmp = sample(1:16571, 9943)
  plot(hclust(dist(t(avg_beta[chr != "24" ,grep ("LCL|Fib" , colnames(avg_beta), invert = T)][tmp,]))), xlab = "", main = paste("Robustness 60%", i))
}

for (i in 1:20){
  tmp = sample(1:16571, 4972)
  plot(hclust(dist(t(avg_beta[chr != "24" ,grep ("LCL|Fib" , colnames(avg_beta), invert = T)][tmp,1:16]))), xlab = "", main = paste("Robustness 30%", i))
}



for (i in 1:20){
  tmp = sample(1:16571, 1658)
  plot(hclust(dist(t(avg_beta[chr != "24" ,grep ("LCL|Fib" , colnames(avg_beta), invert = T)][tmp,1:16]))), xlab = "", main = paste("Robustness 10%", i))
}

dev.off()


#PCA

pdf(file = "PCA.pdf")
Leg = samplenames$Deriv
#All the data
x.pca = prcomp(na.omit(avg_beta), scale = T, center = T)
x.pca.sum = summary(x.pca)
plot(x.pca$rotation[,1], x.pca$rotation[,2], xlab = paste('PC1 (',x.pca.sum$importance[2,1], ')', sep = ''),ylab = paste('PC2 (',x.pca.sum$importance[2,2], ')', sep = ''), main = "PC1/2 all data", col = Leg, pch = 20); legend(x = "topleft", pch = 20, col = c(1:4), c("LCL origin", "Fib origin", "LCL", "Fib"))
plot(x.pca$rotation[,1], x.pca$rotation[,3], xlab = paste('PC1 (',x.pca.sum$importance[2,1], ')', sep = ''),ylab = paste('PC3 (',x.pca.sum$importance[2,3], ')', sep = ''), main = "PC1/3 all data", col = Leg, pch = 20); legend(x = "topright", pch = 20, col = c(1:4), c("LCL origin", "Fib origin", "LCL", "Fib"))
plot(x.pca$rotation[,1], x.pca$rotation[,4], xlab = paste('PC1 (',x.pca.sum$importance[2,1], ')', sep = ''),ylab = paste('PC4 (',x.pca.sum$importance[2,4], ')', sep = ''), main = "PC1/4 all data", col = Leg, pch = 20); legend(x = "topleft", pch = 20, col = c(1:4), c("LCL origin", "Fib origin", "LCL", "Fib"))


ipsc.pca = prcomp(na.omit(avg_beta[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]), scale = T, center =T)
ipsc.pca.sum = summary(ipsc.pca)

dev.off()
#Make sample covar file with only iPSCs
pdf(file = "PCA iPSC only.pdf")
rem = grep ("LCL|Fib" , samplenames$Name)
samplenames.ipsc = samplenames[-rem,]

ipsc_only_leg = c("red", "blue", "orange", "black", "blue", "black", "red", "orange", "blue", "red", "blue", "orange", "red", "black", "orange", "black")
ipsc_only_shape =  c(21, 20,20,20,20,20,20,21, 20,20,21,20,20,20,20,21)

plot(ipsc.pca$rotation[,1], ipsc.pca$rotation[,2], xlab = paste('PC1 (',ipsc.pca.sum$importance[2,1], ')', sep = ''),ylab = paste('PC2 (',ipsc.pca.sum$importance[2,2], ')', sep = ''), main = "PC1/2 iPSC only", col = ipsc_only_leg, pch = ipsc_only_shape);legend(x = "bottomleft", pch = c(20, 20, 20, 20, 20, 21), col = c("red","blue","black","orange", "black","black"), c("0961", "1194", "4280", "8126", "LCL derived", "Fib derived"))
plot(ipsc.pca$rotation[,1], ipsc.pca$rotation[,3], xlab = paste('PC1 (',ipsc.pca.sum$importance[2,1], ')', sep = ''),ylab = paste('PC3 (',ipsc.pca.sum$importance[2,3], ')', sep = ''), main = "PC1/3 iPSC only", col = ipsc_only_leg, pch = ipsc_only_shape); legend(x = "bottomleft", pch = c(20, 20, 20, 20, 20, 21), col = c("red","blue","black","orange", "black","black"), c("0961", "1194", "4280", "8126", "LCL derived", "Fib derived"))
plot(ipsc.pca$rotation[,1], ipsc.pca$rotation[,4], xlab = paste('PC1 (',ipsc.pca.sum$importance[2,1], ')', sep = ''),ylab = paste('PC4 (',ipsc.pca.sum$importance[2,4], ')', sep = ''), main = "PC1/4 iPSC only", col = ipsc_only_leg, pch = ipsc_only_shape); legend(x = "bottomleft", pch = c(20, 20, 20, 20, 20, 21), col = c("red","blue","black","orange", "black","black"), c("0961", "1194", "4280", "8126", "LCL derived", "Fib derived"))


dev.off()

covars = list(batch, type,sex,indiv,pluri,novel,der)

batch.ipsc = samplenames.ipsc$Batch
mef.ipsc = samplenames.ipsc$MEF
type.ipsc = samplenames.ipsc$Type
sex.ipsc =  samplenames.ipsc$Sex
indiv.ipsc = samplenames.ipsc$Indiv
pluri.ipsc =  samplenames.ipsc$Pluri
novel.ipsc = samplenames.ipsc$Novelty
der.ipsc = samplenames.ipsc$Deriv

covars.ipsc = list(mef.ipsc, batch.ipsc,sex.ipsc,indiv.ipsc,pluri.ipsc,novel.ipsc,der.ipsc)

#samplenames.ipsc
lmPCA = function(pca, covars, npcs)
{
  results<-c()
  for (f in covars) {
    for (i in 1:npcs)
    {
      s = summary(lm(pca$rotation[,i]~f));
      results<-c(results,pf(s$fstatistic[[1]],
                            s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
                 s$adj.r.squared)
    }
  }
  resultsM<-matrix(nrow = length(covars), ncol = 2*npcs, data =
                     results, byrow = TRUE)
  resultsM
  
  
}

pcaresults = lmPCA(x.pca,covars,4)
rownames(pcaresults) = c("batch", "type","sex","indiv","pluri","novel","der")
colnames(pcaresults) = c("PC1 pval", "PC1 adj R sqs","PC2 pval", "PC2 adj R sqs","PC3 pval", "PC3 adj R sqs","PC4 pval", "PC4 adj R sqs")

pcaresults.ipsc = lmPCA(ipsc.pca,covars.ipsc,4)
rownames(pcaresults.ipsc) = c("mef","batch","sex","indiv","pluri","novel","der")
colnames(pcaresults.ipsc) = c("PC1 pval", "PC1 adj R sqs","PC2 pval", "PC2 adj R sqs","PC3 pval", "PC3 adj R sqs","PC4 pval", "PC4 adj R sqs")

#If you want to look at PCA after subtracting the mean
avg_beta = stand
ipsc.pca = prcomp(na.omit(avg_beta[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]), scale = T, center =T)
ipsc.pca.sum = summary(ipsc.pca)

pcaresults.ipsc.stan = lmPCA(ipsc.pca,covars.ipsc,4)
rownames(pcaresults.ipsc.stan) = c("mef","batch","sex","indiv","pluri","novel","der")
colnames(pcaresults.ipsc.stan) = c("PC1 pval", "PC1 adj R sqs","P2 pval", "PC2 adj R sqs","PC3 pval", "PC3 adj R sqs","PC4 pval", "PC4 adj R sqs")

covarcor = function(covars)
{
  results<-c()
  for (f in covars) {
  
      s = cor(covar[f], covar);
      results<-c(results,pf(s$fstatistic[[1]],
                            s$fstatistic[[2]],s$fstatistic[[3]], lower.tail = FALSE),
                 s$adj.r.squared)
    }
  }
  resultsM<-matrix(nrow = length(covars), ncol = 2*npcs, data =
                     results, byrow = TRUE)
  resultsM
  
  
}

###DE Analysis Code#####
library(limma)

#meth.final = expr_quant.all
meth.final = expr_gene

labs = c("OF", "OL", "OL", "LCL", "OL", "OL", "OL", "Fib", "OL", "Fib", 
  "LCL", "OF", "Fib", "LCL", "OL", "LCL", "Fib", "OL", "OF", "OL", 
  "OL", "OL", "OL", "OF")

design <-(model.matrix(~0+labs))

colnames(design) <- c("Fib", "LCL", "OF", "OL")

fit  <- lmFit(meth.final, design)
fit <- eBayes(fit)

cm <- makeContrasts(
  OLvsOF = OL-OF,
  OLvsLCL = OL-LCL,
  OFvsFib = OF-Fib,
  LCLvsFib = LCL-Fib,
  levels=design)

fit2 <- contrasts.fit(fit, cm)
fit2 <-eBayes(fit2)

iPSC_DMR <- topTable(fit2, coef=1, adjust="BH", number=Inf, sort.by="p")
LCLs_vs_iPSC.L <- topTable(fit2, coef=2, adjust="BH", number=Inf, sort.by="p")
Fibs_vs_iPSC.F <- topTable(fit2, coef=3, adjust="BH", number=Inf, sort.by="p")
LCL_vs_Fibs <- topTable(fit2, coef=4, adjust="BH", number=Inf, sort.by="p")

overlap_iPSC_DMR_Origin_DMRs = iPSC_DMR[iPSC_DMR$adj.P.Val < 0.05 , ][rownames(iPSC_DMR[iPSC_DMR$adj.P.Val < 0.05 ,]) %in% rownames(LCL_vs_Fibs[LCL_vs_Fibs$adj.P.Val < 0.05 ,]) , ]

dist_of_Ps = iPSC_DMR[rownames(iPSC_DMR) %in% rownames(LCL_vs_Fibs[LCL_vs_Fibs$adj.P.Val < 0.05 , ]) , ]

head(iPSC_DMR[iPSC_DMR$adj.P.Val < 0.05 , ])
dim(LCL_vs_Fibs[LCL_vs_Fibs$adj.P.Val < 0.05 , ])
dim(iPSC_DMR[iPSC_DMR$adj.P.Val < 0.05 , ])
dim(LCLs_vs_iPSC.L[LCLs_vs_iPSC.L$adj.P.Val < 0.05 , ])
dim(Fibs_vs_iPSC.F[Fibs_vs_iPSC.F$adj.P.Val < 0.05 , ])

write.table(LCL_vs_Fibs,'DE_LCLvFib_FDR5.txt', quote=F, sep = '\t')
write.table(LCLs_vs_iPSC.L,'DE_LCLviPSC_FDR5.txt', quote=F, sep = '\t')
write.table(Fibs_vs_iPSC.F,'DE_FibviPSC_FDR5.txt', quote=F, sep = '\t')
write.table(iPSC_DMR, 'DE_F.iPSCvL.iPSC_FDR5.txt', quote=F, sep = "\t")

pdf(file='Histograms.pdf')
hist(iPSC_DMR$P.Value, main = "Distribution of L-iPSC vs F-iPSC DE P-values", xlab= "P-value")
hist(LCL_vs_Fibs$P.Value, main = "Distribution of LCL vs Fibroblast DE P-values", xlab= "P-value")
hist(LCLs_vs_iPSC.L$P.Value, main = "Distribution of L-iPSCs vs LCLs DE P-values", xlab= "P-value")
hist(Fibs_vs_iPSC.F$P.Value, main = "Distribution of F-iPSC vs Fibroblasts DE P-values", xlab= "P-value")
dev.off()

#Do DE with subsets of genes
LvF = LCL_vs_Fibs[LCL_vs_Fibs$adj.P.Val < 0.05 , 1:2]
LvI = LCLs_vs_iPSC.L[LCLs_vs_iPSC.L$adj.P.Val < 0.05 , 1:2]
FvI = Fibs_vs_iPSC.F[Fibs_vs_iPSC.F$adj.P.Val < 0.05 , 1:2]
LvF_genes = meth.final[ rownames(meth.final) %in% rownames(LvF),]
LvI_genes = meth.final[ rownames(meth.final) %in% rownames(LvI),]
FvI_genes = meth.final[ rownames(meth.final) %in% rownames(FvI),]

#Using only LvF genes
cm3 <- makeContrasts(
  OLvsOF = OL-OF,
  levels=design)
fit  <- lmFit(LvF_genes, design)
fit <- eBayes(fit)
fit2 <- contrasts.fit(fit, cm3)
fit2 <-eBayes(fit2)

iPSC_DMR_LvFgenes <- topTable(fit2, coef=1, adjust="BH", number=Inf, sort.by="p")

#Using only LvI genes
fit  <- lmFit(LvI_genes, design)
fit <- eBayes(fit)
fit2 <- contrasts.fit(fit, cm3)
fit2 <-eBayes(fit2)
iPSC_DMR_LvIgenes <- topTable(fit2, coef=1, adjust="BH", number=Inf, sort.by="p")

#Using only FvI genes
fit  <- lmFit(FvI_genes, design)
fit <- eBayes(fit)
fit2 <- contrasts.fit(fit, cm3)
fit2 <-eBayes(fit2)
iPSC_DMR_FvIgenes <- topTable(fit2, coef=1, adjust="BH", number=Inf, sort.by="p")

##QQ Plot
pvals = iPSC_DMR$P.Value
observed <- sort(pvals)
lobs <- -(log10(observed))

expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))

pdf("QQplot_FDRLine_FDR5DE.pdf", width=6, height=6)
plot(c(0,7), c(0,7), col="red", lwd=1, type="l", main = "QQ Plot by DE Gene Subsets",xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=19, cex=.6, col="black") 

pvalsLvF = iPSC_DMR_LvFgenes$P.Value
observedLvF <- sort(pvalsLvF)
lobsLvF <- -(log10(observedLvF))
expectedLvF <- c(1:length(observedLvF)) 
lexpLvF <- -(log10(expectedLvF / (length(expectedLvF)+1)))
points(lexpLvF, lobsLvF, pch=19, cex=.6, col="green") 

pvalsLvI = iPSC_DMR_LvIgenes$P.Value
observedLvI <- sort(pvalsLvI)
lobsLvI <- -(log10(observedLvI))
expectedLvI <- c(1:length(observedLvI)) 
lexpLvI <- -(log10(expectedLvI / (length(expectedLvI)+1)))
points(lexpLvI, lobsLvI, pch=19, cex=.6, col="cyan3") 

pvalsFvI = iPSC_DMR_FvIgenes$P.Value
observedFvI <- sort(pvalsFvI)
lobsFvI <- -(log10(observedFvI))
expectedFvI <- c(1:length(observedFvI)) 
lexpFvI <- -(log10(expectedFvI / (length(expectedFvI)+1)))
points(lexpFvI, lobsFvI, pch=19, cex=.6, col="darkorange") 

legend(x = "topleft", pch = 19, col = c("black", "green", "cyan3", "darkorange"), c("All Genes", "Genes DE LCL v Fib", "Genes DE LCL v L-iPSC", "Genes DE Fib v F-iPSC"))
lines(c(0,7), c(log10(.95/.05),(log10(.95/.05)+7)))
dev.off()

###DE Analysis Code#####
##By individual
library(limma)

meth.final = expr_gene
iPSC = meth.final[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]
iPSC_noX = iPSC[chr[,1] != "23" ,]


labs.ind = as.character(samplenames.ipsc$Indiv)

design.ind <-(model.matrix(~0+labs.ind))

colnames(design.ind) <- c("ind1", "ind2", "ind3", "ind4")

fit.ind  <- lmFit(iPSC_noX, design.ind)
fit.ind <- eBayes(fit.ind)

cm.ind <- makeContrasts(
  X1V2 = ind1-ind2,
  X1V3 = ind1-ind3,
  X1V4 = ind1-ind4,
  X2V3 = ind2-ind3,
  X2V4 = ind2-ind4,
  X3V4 = ind3-ind4,
  levels=design.ind)

fit.ind2 <- contrasts.fit(fit.ind, cm.ind)
fit.ind2 <-eBayes(fit.ind2)

X1V2 = topTable(fit.ind2, coef=1, adjust="BH", number=Inf, sort.by="p")
X1V3 = topTable(fit.ind2, coef=2, adjust="BH", number=Inf, sort.by="p")
X1V4 = topTable(fit.ind2, coef=3, adjust="BH", number=Inf, sort.by="p")
X2V3 = topTable(fit.ind2, coef=4, adjust="BH", number=Inf, sort.by="p")
X2V4 = topTable(fit.ind2, coef=5, adjust="BH", number=Inf, sort.by="p")
X3V4 = topTable(fit.ind2, coef=6, adjust="BH", number=Inf, sort.by="p")

dim(X1V2[X1V2$adj.P.Val < 0.01 , ])
#76 #72 noX
dim(X1V3[X1V3$adj.P.Val < 0.01 , ])
#87 #86 noX
dim(X1V4[X1V4$adj.P.Val < 0.01 , ])
#77 #69 noX
dim(X2V3[X2V3$adj.P.Val < 0.01 , ])
#112 #107 noX
dim(X2V4[X2V4$adj.P.Val < 0.01 , ])
#88 #83 noX
dim(X3V4[X3V4$adj.P.Val < 0.01 , ])
#119 #110 noX


##Subset out based on Sammy's variance genes##
sammy= read.table('Variance iPSCs by gene.txt', as.is=T, header=T)
filt_sammy = sammy[which(rownames(sammy) %in% rownames(expr_gene)),]
filt_sammy_allvar = filt_sammy[order(filt_sammy[,1],decreasing=TRUE),]
dim(filt_sammy)
#[1] 12549     2
#10% 
sammy10 = filt_sammy[1:1255,]
gene10 = expr_gene[which(rownames(sammy10) %in% rownames(expr_gene)),]

#Confirm the right rows were pulled
test = rbind(rownames(sammy10), rownames(gene10))
dim(unique(test))

plot(hclust(as.dist(1-cor(as.matrix(gene10[,grep ("LCL|Fib" , colnames(gene10), invert = T)])))),main = "Dendogram using only the 10% most variant genes")

#20% 
sammy20 = filt_sammy[1:2510,]
gene20 = expr_gene[which(rownames(sammy20) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene20[,grep ("LCL|Fib" , colnames(gene20), invert = T)])))), main = "Dendogram using only the 20% most variant genes")

#30% 
sammy30 = filt_sammy[1:3765,]
gene30 = expr_gene[which(rownames(sammy30) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene30[,grep ("LCL|Fib" , colnames(gene30), invert = T)])))),main = "Dendogram using only the 30% most variant genes")

#40% 
sammy40 = filt_sammy[1:5020,]
gene40 = expr_gene[which(rownames(sammy40) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene40[,grep ("LCL|Fib" , colnames(gene40), invert = T)])))),main = "Dendogram using only the 40% most variant genes")

#50% 
sammy50 = filt_sammy[1:6275,]
gene50 = expr_gene[which(rownames(sammy50) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene50[,grep ("LCL|Fib" , colnames(gene50), invert = T)])))),main = "Dendogram using only the 50% most variant genes")

#60% 
sammy60 = filt_sammy[1:7530,]
gene60 = expr_gene[which(rownames(sammy60) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene60[,grep ("LCL|Fib" , colnames(gene60), invert = T)])))),main = "Dendogram using only the 60% most variant genes")

#70% 
sammy70 = filt_sammy[1:8785,]
gene70 = expr_gene[which(rownames(sammy70) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene70[,grep ("LCL|Fib" , colnames(gene70), invert = T)])))),main = "Dendogram using only the 70% most variant genes")

#80% 
sammy80 = filt_sammy[1:10040,]
gene80 = expr_gene[which(rownames(sammy80) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene80[,grep ("LCL|Fib" , colnames(gene80), invert = T)])))),main = "Dendogram using only the 80% most variant genes")

#90% 
sammy90 = filt_sammy[1:11294,]
gene90 = expr_gene[which(rownames(sammy90) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene90[,grep ("LCL|Fib" , colnames(gene90), invert = T)])))),main = "Dendogram using only the 90% most variant genes")

#100% 
sammy100 = filt_sammy[1:12549,]
gene100 = expr_gene[which(rownames(sammy100) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene100[,grep ("LCL|Fib" , colnames(gene100), invert = T)])))),main = "Dendogram using only the 100% most variant genes")

sammy10var = filt_sammy_allvar[1:1255,]
gene10var = expr_gene[which(rownames(sammy10var) %in% rownames(expr_gene)),]

##Filter based on Sammy's variance
#Confirm the right rows were pulled
test = rbind(rownames(sammy10var), rownames(gene10var))
dim(unique(test))

plot(hclust(as.dist(1-cor(as.matrix(gene10var[,grep ("LCL|Fib" , colnames(gene10var), invert = T)])))),main = "Dendogram using only the 10var% most variant genes")

#20var% 
sammy20var = filt_sammy_allvar[1:2510,]
gene20var = expr_gene[which(rownames(sammy20var) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene20var[,grep ("LCL|Fib" , colnames(gene20var), invert = T)])))), main = "Dendogram using only the 20var% most variant genes")

#30var% 
sammy30var = filt_sammy_allvar[1:3765,]
gene30var = expr_gene[which(rownames(sammy30var) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene30var[,grep ("LCL|Fib" , colnames(gene30var), invert = T)])))),main = "Dendogram using only the 30var% most variant genes")

#40var% 
sammy40var = filt_sammy_allvar[1:5020,]
gene40var = expr_gene[which(rownames(sammy40var) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene40var[,grep ("LCL|Fib" , colnames(gene40var), invert = T)])))),main = "Dendogram using only the 40var% most variant genes")

#50var% 
sammy50var = filt_sammy_allvar[1:6275,]
gene50var = expr_gene[which(rownames(sammy50var) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene50var[,grep ("LCL|Fib" , colnames(gene50var), invert = T)])))),main = "Dendogram using only the 50var% most variant genes")

#60var% 
sammy60var = filt_sammy_allvar[1:7530,]
gene60var = expr_gene[which(rownames(sammy60var) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene60var[,grep ("LCL|Fib" , colnames(gene60var), invert = T)])))),main = "Dendogram using only the 60var% most variant genes")

#70var% 
sammy70var = filt_sammy_allvar[1:8785,]
gene70var = expr_gene[which(rownames(sammy70var) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene70var[,grep ("LCL|Fib" , colnames(gene70var), invert = T)])))),main = "Dendogram using only the 70var% most variant genes")

#80var% 
sammy80var = filt_sammy_allvar[1:10040,]
gene80var = expr_gene[which(rownames(sammy80var) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene80var[,grep ("LCL|Fib" , colnames(gene80var), invert = T)])))),main = "Dendogram using only the 80var% most variant genes")

#90var% 
sammy90var = filt_sammy_allvar[1:11294,]
gene90var = expr_gene[which(rownames(sammy90var) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene90var[,grep ("LCL|Fib" , colnames(gene90var), invert = T)])))),main = "Dendogram using only the 90var% most variant genes")

#100% 
sammy100var = filt_sammy_allvar[1:12549,]
gene100var = expr_gene[which(rownames(sammy100var) %in% rownames(expr_gene)),]

plot(hclust(as.dist(1-cor(as.matrix(gene100var[,grep ("LCL|Fib" , colnames(gene100var), invert = T)])))),main = "Dendogram using only the 100var% most variant genes")


#############

##Regress out pluri to see if this fixes dendogram##
pluri.residual.int = matrix(nrow= nrow(expr_quant.all[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]), ncol = ncol(expr_quant.all[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]))
rownames(pluri.residual.int) = rownames(expr_quant.all[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)])
colnames(pluri.residual.int) = colnames(expr_quant.all[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)])

expr_quant = expr_quant.all[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]
expr_gene.i = expr_gene[,grep ("LCL|Fib", colnames(expr_gene),invert=T)]
for (i in 1:nrow(expr_quant)) {
  model= lm(expr_quant[i,]~ samplenames.ipsc$Pluri)
  pluri.residual.int[i,] = resid(model) + model$coefficients[1]
}
plot(hclust(as.dist(1-cor(as.matrix(pluri.residual.int)))))

##Regress out array to see if this fixes dendogram##
colnames(expr_gene) = samplenames$Name
batch.residual = matrix(nrow= nrow(expr_gene), ncol = ncol(expr_gene))
rownames(batch.residual) = rownames(expr_gene)
colnames(batch.residual) = colnames(expr_gene)

for (i in 1:nrow(expr_gene)) {
  model= lm(expr_gene[i,]~ samplenames$Batch)
  batch.residual[i,] = resid(model)
  }

plot(hclust(as.dist(1-cor(as.matrix(batch.residual)))),main = "Dendrogram reg array (no intercept)", xlab="Pearson")
plot(hclust(dist(t(batch.residual))), xlab = "hclust/Euclidean distance", main = "Dendrogram reg array (no intercept)")
plot(hclust(dist(t(batch.residual), method = "manhattan")), xlab = "Manhattan", main = "Dendrogram reg array (no intercept)")

#Add intercept back
batch.residual.int = matrix(nrow= nrow(expr_gene), ncol = ncol(expr_gene))
rownames(batch.residual.int) = rownames(expr_gene)
colnames(batch.residual.int) = colnames(expr_gene)

for (i in 1:nrow(expr_gene)) {
  model= lm(expr_gene[i,]~ samplenames$Batch)
  batch.residual.int[i,] = resid(model) + model$coefficients[1]
}

plot(hclust(as.dist(1-cor(as.matrix(batch.residual.int)))),main = "Dendrogram reg array (w/intercept)", xlab="Pearson")
plot(hclust(dist(t(batch.residual.int))), xlab = "hclust/Euclidean distance", main = "Dendrogram reg array (w/intercept)")
plot(hclust(dist(t(batch.residual.int), method = "manhattan")), xlab = "Manhattan", main = "Dendrogram reg array (w/intercept)")

##Subtract mean and then regres sout batch (similar to Sammy's exp)
stan = apply(expr_gene,1, function(x) x-mean(x))
stand = t(stan)

batch.residual.int.ms = matrix(nrow= nrow(stand), ncol = ncol(stand))
rownames(batch.residual.int.ms) = rownames(stand)
colnames(batch.residual.int.ms) = colnames(stand)

for (i in 1:nrow(stand)) {
  model= lm(stand[i,]~ samplenames$Batch)
  batch.residual.int.ms[i,] = resid(model) + model$coefficients[1]
}

plot(hclust(as.dist(1-cor(as.matrix(batch.residual.int.ms)))),main = "Dendrogram reg array - ms (w/intercept)", xlab="Pearson")
plot(hclust(dist(t(batch.residual.int.ms))), xlab = "hclust/Euclidean distance", main = "Dendrogram reg array - ms (w/intercept)")
plot(hclust(dist(t(batch.residual.int.ms), method = "manhattan")), xlab = "Manhattan", main = "Dendrogram reg array - ms (w/intercept)")

#Do not add intercept back in
batch.residual.ms = matrix(nrow= nrow(stand), ncol = ncol(stand))
rownames(batch.residual.ms) = rownames(stand)
colnames(batch.residual.ms) = colnames(stand)

for (i in 1:nrow(stand)) {
  model= lm(stand[i,]~ samplenames$Batch)
  batch.residual.ms[i,] = resid(model)
  }

plot(hclust(as.dist(1-cor(as.matrix(batch.residual.ms)))),main = "Dendrogram reg array - ms (no intercept)", xlab="Pearson")
plot(hclust(dist(t(batch.residual.ms))), xlab = "hclust/Euclidean distance", main = "Dendrogram reg array - ms (no intercept)")
plot(hclust(dist(t(batch.residual.ms), method = "manhattan")), xlab = "Manhattan", main = "Dendrogram reg array - ms(no intercept)")

##Regress out PC1 to see if this fixes dendogram##
PC.residual.int = matrix(nrow= nrow(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]), ncol = ncol(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]))
rownames(PC.residual.int) = rownames(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)])
colnames(PC.residual.int) = colnames(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)])

sum.pca = prcomp(na.omit(expr_gene[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]), scale = T, center =T)
PCSf = sum.pca$rotation

for (i in 1:nrow(expr_gene)) {
  model= lm(expr_gene[i,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]~ sum.pca$rotation[,1:2])
  PC.residual.int[i,] = resid(model)
  }
plot(hclust(dist(t(PC.residual.int))), xlab = "hclust/Euclidean distance", main = "Dendrogram reg PC1 +PC2")
plot(hclust(as.dist(1-cor(as.matrix(PC.residual.int)))),main = "Dendrogram reg PC1+PC2", xlab="Pearson")
plot(hclust(dist(t(PC.residual.int), method="manhattan")), xlab = "Manhattan", main = "Dendrogram reg PC1 +PC2")

##Regress out PC1 from mean subtracted data
PC.residual.ms = matrix(nrow= nrow(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]), ncol = ncol(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]))
rownames(PC.residual.ms) = rownames(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)])
colnames(PC.residual.ms) = colnames(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)])

sum.pca = prcomp(na.omit(expr_gene[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]), scale = T, center =T)
PCSf = sum.pca$rotation

for (i in 1:nrow(expr_gene)) {
  model= lm(stand[i,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]~ sum.pca$rotation[,1])
  PC.residual.ms[i,] = resid(model)
}
plot(hclust(dist(t(PC.residual.ms))), xlab = "hclust/Euclidean distance", main = "Dendrogram reg PC1 - ms")
plot(hclust(as.dist(1-cor(as.matrix(PC.residual.ms)))),main = "Dendrogram reg PC1 - ms", xlab="Pearson")
plot(hclust(dist(t(PC.residual.ms), method="manhattan")), xlab = "Manhattan", main = "Dendrogram reg PC1 - ms")

#Regress out PC1+2 mean subtracted data
PC1_2.residual.ms = matrix(nrow= nrow(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]), ncol = ncol(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]))
rownames(PC1_2.residual.ms) = rownames(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)])
colnames(PC1_2.residual.ms) = colnames(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)])

sum.pca = prcomp(na.omit(expr_gene[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]), scale = T, center =T)
PCSf = sum.pca$rotation

for (i in 1:nrow(expr_gene)) {
  model= lm(stand[i,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]~ sum.pca$rotation[,1:2])
  PC1_2.residual.ms[i,] = resid(model)
}
plot(hclust(dist(t(PC1_2.residual.ms))), xlab = "hclust/Euclidean distance", main = "Dendrogram reg PC1 +PC2 - ms")
plot(hclust(as.dist(1-cor(as.matrix(PC1_2.residual.ms)))),main = "Dendrogram reg PC1+PC2 - ms", xlab="Pearson")
plot(hclust(dist(t(PC1_2.residual.ms), method="manhattan")), xlab = "Manhattan", main = "Dendrogram reg PC1 +PC2 - ms")


#Regress out PC1-3
PC1_3.residual = matrix(nrow= nrow(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]), ncol = ncol(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]))
rownames(PC1_3.residual) = rownames(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)])
colnames(PC1_3.residual) = colnames(expr_gene[,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)])

sum.pca = prcomp(na.omit(expr_gene[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]), scale = T, center =T)
PCSf = sum.pca$rotation

for (i in 1:nrow(expr_gene)) {
  model= lm(expr_gene[i,grep ("LCL|Fib", colnames(expr_quant.all),invert=T)]~ sum.pca$rotation[,1:2])
  PC1_3.residual[i,] = resid(model)
}
plot(hclust(dist(t(PC1_3.residual))), xlab = "hclust/Euclidean distance", main = "Dendrogram reg PC1-3")
plot(hclust(as.dist(1-cor(as.matrix(PC1_3.residual)))),main = "Dendrogram reg PC1-3", xlab="Pearson")
plot(hclust(dist(t(PC1_3.residual), method="manhattan")), xlab = "Manhattan", main = "Dendrogram reg PC1-3")


##Subset randomly down to three individuals
ipsc_expr = expr_gene[,grep ("LCL|Fib" , colnames(expr_gene), invert = T)]
S0961 = ipsc_expr[,grep ("S0961" , colnames(ipsc_expr))]
rand_S0961 = S0961[,sample(1:4, 3,replace=FALSE)]

S1194 = ipsc_expr[,grep ("S1194" , colnames(ipsc_expr))]
rand_S1194 = S1194[,sample(1:4, 3,replace=FALSE)]

S4280 = ipsc_expr[,grep ("S4280" , colnames(ipsc_expr))]
rand_4280 = S4280[,sample(1:4, 3,replace=FALSE)]

S8126 = ipsc_expr[,grep ("S8126" , colnames(ipsc_expr))]
rand_S8126 = S8126[,sample(1:4, 3,replace=FALSE)]

subset_expr_ipsc = cbind(rand_S0961, rand_S1194, rand_4280,rand_S8126)

plot(hclust(dist(t(subset_expr_ipsc))), xlab = "", main = "Dendrogram with subsetted iPSCs")
plot(hclust(dist(t(expr_gene[,grep ("LCL|Fib|.F" , colnames(expr_gene), invert = T)]))), xlab = "", main = "Dendrogram with no fibro")


#To get the correlation between indiv vs non-indiv
library(gplots)
cormatrix <- cor(expr_gene.i)
varmatrix = var(expr_gene.i)
names = c("S0961", "S1194", "S8126", "S4280", "S1194", 
          "S4280", "S0961", "S8126", "S1194", "S0961", 
          "S1194", "S8126", "S0961", "S4280", "S8126", 
          "S4280")
colnames(cormatrix) = names
rownames(cormatrix) = names

colnames(varmatrix) = names
rownames(varmatrix) = names

replicates=c()
non.replicates=c()

for(i in 1:15){
  for(j in (i+1):16){
    if(names[i]==names[j]){
      replicates=c(replicates,cormatrix[i,j])
    }
    else{
      non.replicates=c(non.replicates,cormatrix[i,j])
    }
  }
}
boxplot2(replicates,non.replicates, main = "Correlation of Samples", ylab = "Correlation", xlab = "Replicates vs Non-Replicates")

varmatrix = var(expr_gene.i)
var.m = varmatrix[,1]

boxplot2(var.m~names, main = "Variance of Individuals", ylab = "Variance", xlab = "Individuals")

varmatrix = var(pluri.residual.int)
var.m = varmatrix[,1]

col.box = c("red", "blue", "blue", 
            "blue", "blue", "blue", "blue", "red", 
            "blue", "blue", "red", "blue", "blue", 
            "blue", "blue", "red")
varorder = c("S0961", "S1194", "S8126", "S4280", 
               "S1194", "S4280", "S0961", "S8126", "S1194", "S0961", "S1194", 
               "S8126", "S0961", "S4280", "S8126", "S4280")
boxplot2(var.m~varorder, main = "Variance of individuals after regressing Pluri score", ylab = "Variance", xlab = "Individuals")
beeswarm(var.m~varorder,add=T, col=col.box, vertical=T,pch=20)

###### Irene venn diagram code
library(VennDiagram)

probes <- data.frame(rownames(iPSC_DMR))
names(probes) <- "probes"

make.venn.quad <- function(geneset1, geneset2, geneset3, geneset4, geneset1.label, geneset2.label, geneset3.label, geneset4.label, univ){
  univ$g1 <- univ$probes %in% geneset1
  univ$g2 <- univ$probes %in% geneset2
  univ$g3 <- univ$probes %in% geneset3 
  univ$g4 <- univ$probes %in% geneset4 
  #pdf(file=paste(prefix, ".pdf", sep=""), width=7, height=7)
  venn.placeholder <- draw.quad.venn(length(geneset1),length(geneset2), length(geneset3), length(geneset4), dim(univ[univ$g1 == T & univ$g2 == T , ])[1], dim(univ[univ$g1 == T & univ$g3 == T , ])[1], dim(univ[univ$g1 == T & univ$g4 == T , ])[1], dim(univ[univ$g2 == T & univ$g3 == T , ])[1], dim(univ[univ$g2 == T & univ$g4 == T , ])[1], dim(univ[univ$g3 == T & univ$g4 == T , ])[1], dim(univ[univ$g1 == T & univ$g2 == T & univ$g3 == T , ])[1], dim(univ[univ$g1 == T & univ$g2 == T & univ$g4 == T , ])[1], dim(univ[univ$g1 == T & univ$g3 == T & univ$g4 == T , ])[1], dim(univ[univ$g2 == T & univ$g3 == T & univ$g4 == T , ])[1],  dim(univ[univ$g1 == T & univ$g2 == T & univ$g3 == T & univ$g4 == T , ])[1], c(geneset1.label, geneset2.label, geneset3.label, geneset4.label), fill=c("goldenrod", "plum4", "steelblue3", "darkolivegreen3"), alpha=c(0.5, 0.5, 0.5, 0.5),col=NA, euler.d=T)
  complement.size <- dim(univ[univ$g1 == F & univ$g2 == F & univ$g3 == F & univ$g4 == F , ])[1]
  grid.text(paste(complement.size, " not DE in any", sep=""), x=0.2, y=0.08)
  #dev.off()
}

dev.off()
# Make venn of full
pdf(file = "VennDiagram_DE_FDR5.pdf")
make.venn.quad(rownames(iPSC_DMR[iPSC_DMR$adj.P.Val < 0.05 , ]), rownames(LCL_vs_Fibs[LCL_vs_Fibs$adj.P.Val < 0.05 , ]), rownames(LCLs_vs_iPSC.L[LCLs_vs_iPSC.L$adj.P.Val < 0.05 , ]), rownames(Fibs_vs_iPSC.F[Fibs_vs_iPSC.F$adj.P.Val < 0.05 , ]), paste("DEs iPSC", dim(iPSC_DMR[iPSC_DMR$adj.P.Val < 0.05 , ])[1]), paste("DEs Origins", dim(LCL_vs_Fibs[LCL_vs_Fibs$adj.P.Val < 0.05 , ])[1]), paste("DEs LCLs", dim(LCLs_vs_iPSC.L[LCLs_vs_iPSC.L$adj.P.Val < 0.05 , ])[1]), paste("DEs Fibs", dim(Fibs_vs_iPSC.F[Fibs_vs_iPSC.F$adj.P.Val < 0.05 , ])[1]), probes)
dev.off()

##Venn of probe inclusion scheme
probes <- data.frame(rownames(norm_quant.all))
names(probes) <- "probes"
make.venn.quad <- function(geneset1, geneset2, geneset3, geneset4, geneset1.label, geneset2.label, geneset3.label, geneset4.label, univ){
  univ$g1 <- univ$probes %in% geneset1
  univ$g2 <- univ$probes %in% geneset2
  univ$g3 <- univ$probes %in% geneset3 
  univ$g4 <- univ$probes %in% geneset4 
  venn.placeholder <- draw.quad.venn(cex = rep(2),cat.cex =rep(1),length(geneset1),length(geneset2), length(geneset3), length(geneset4), dim(univ[univ$g1 == T & univ$g2 == T , ])[1], dim(univ[univ$g1 == T & univ$g3 == T , ])[1], dim(univ[univ$g1 == T & univ$g4 == T , ])[1], dim(univ[univ$g2 == T & univ$g3 == T , ])[1], dim(univ[univ$g2 == T & univ$g4 == T , ])[1], dim(univ[univ$g3 == T & univ$g4 == T , ])[1], dim(univ[univ$g1 == T & univ$g2 == T & univ$g3 == T , ])[1], dim(univ[univ$g1 == T & univ$g2 == T & univ$g4 == T , ])[1], dim(univ[univ$g1 == T & univ$g3 == T & univ$g4 == T , ])[1], dim(univ[univ$g2 == T & univ$g3 == T & univ$g4 == T , ])[1],  dim(univ[univ$g1 == T & univ$g2 == T & univ$g3 == T & univ$g4 == T , ])[1], c(geneset1.label, geneset2.label, geneset3.label, geneset4.label), fill=c("goldenrod", "plum4", "steelblue3", "darkolivegreen3"), alpha=c(0.5, 0.5, 0.5, 0.5),col=NA, euler.d=T)
  complement.size <- dim(univ[univ$g1 == F & univ$g2 == F & univ$g3 == F & univ$g4 == F , ])[1]
  grid.text(paste(complement.size, " not detected in any", sep=""), x=0.2, y=0.08)
  }

dev.off()
# Make venn of full
pdf(file = "VennDiagram_ProbeDetection.pdf", width=12.5, height=12)
make.venn.quad(probes[detect.Fib,], probes[detect.LCL,], probes[detect.FiPSC,], probes[detect.LiPSC,], paste("Detected in Fibroblasts", length(detect.Fib)), paste("Detected in LCLs", length(detect.LCL)), paste("Detected in F-iPSCs", length(detect.FiPSC)), paste("Detected in L-iPSCs", length(detect.LiPSC)), probes)
dev.off()


####### Boxplot of DMPs between L-iPSCs and F-iPSCs ordered by genomic location
pdf(file = "iPSC_DE.pdf")
library(beeswarm)

iPSC_DMR_loc = expr_gene[which(rownames(expr_gene) == "TSTD1"),]
boxplot(iPSC_DMR_loc~samplenames$Type, main = "Expression of TSTD1")
beeswarm(iPSC_DMR_loc~samplenames$Type, add=T, col=2, pwcol = samplenames$Deriv, vertical=T,pch=20)
legend(x = "topright", pch = 20, col = c(1:4), c("LCL origin", "Fib origin", "LCL", "Fib"))

################## Varience explained
samplenames = read.table('covar.txt', header=T, sep ='\t')
#Re-order samplenames based on array location
samplenames = samplenames[order(samplenames$Order),]

rem = grep ("LCL|Fib" , samplenames$Name)
samplenames.ipsc = samplenames[-rem,]
origin_type= as.factor(samplenames.ipsc$Deriv)
ind = as.factor(samplenames.ipsc$Indiv)
meth.final = expr_gene_nosex

var.resid.err = matrix(ncol = 1, nrow = dim(meth.final)[1])
var.origin = matrix(ncol = 1, nrow = dim(meth.final)[1])
var.ind = matrix(ncol = 1, nrow = dim(meth.final)[1])

for (i in 1:dim(meth.final)[1]){
  #for(i in 1:100){
    tmp <- lm(unlist(meth.final[,grep ("LCL|Fib" , colnames(meth.final), invert = T)][i,]) ~ ind + origin_type)
  var.ind[i] <- anova(tmp)[1,2]/sum(anova(tmp)[,2])
  var.origin[i] <- anova(tmp)[2,2]/sum(anova(tmp)[,2])
  var.resid.err[i] <- anova(tmp)[3,2]/sum(anova(tmp)[,2])
}

p.origin = matrix(ncol = 1, nrow = dim(meth.final)[1])
p.ind = matrix(ncol = 1, nrow = dim(meth.final)[1])

for (i in 1:dim(meth.final)[1]){
  #for(i in 1:100){
  tmp <- lm(unlist(meth.final[,grep ("LCL|Fib" , colnames(meth.final), invert = T)][i,]) ~ ind + origin_type)
  p.ind[i] <- anova(tmp)[1,5]
  p.origin[i] <- anova(tmp)[2,5]
  
}
hist(p.origin)
hist(p.ind)

var.in.or = cbind(var.ind, var.origin, var.resid.err)

#Permute data - check each covariate seperatley
var.resid.err.p = matrix(ncol = 1, nrow = dim(meth.final)[1])
var.origin.p = matrix(ncol = 1, nrow = dim(meth.final)[1])
var.ind.p = matrix(ncol = 1, nrow = dim(meth.final)[1])

for (i in 1:dim(meth.final)[1]){
  #for(i in 1:100){
  perm = sample(meth.final[,grep ("LCL|Fib" , colnames(meth.final), invert = T)][i,])
  tmp <- lm(unlist(perm) ~ ind)
  var.ind.p[i] <- anova(tmp)[1,2]/sum(anova(tmp)[,2])
  var.resid.err.p[i] <- anova(tmp)[2,2]/sum(anova(tmp)[,2])
}
var.in.or.p = cbind(var.ind.p, var.resid.err.p)
mean(var.ind.p)
boxplot(var.in.or.p)

var.resid.err = matrix(ncol = 1, nrow = dim(meth.final)[1])
var.ind = matrix(ncol = 1, nrow = dim(meth.final)[1])

for (i in 1:dim(meth.final)[1]){
  #for(i in 1:100){
  tmp <- lm(unlist(meth.final[,grep ("LCL|Fib" , colnames(meth.final), invert = T)][i,]) ~ ind)
  var.ind[i] <- anova(tmp)[1,2]/sum(anova(tmp)[,2])
  var.resid.err[i] <- anova(tmp)[2,2]/sum(anova(tmp)[,2])
}

var.ind.adjr.p = matrix(ncol = 1, nrow = dim(meth.final)[1])

for (i in 1:dim(meth.final)[1]){
  #for(i in 1:100){
  perm = sample(meth.final[,grep ("LCL|Fib" , colnames(meth.final), invert = T)][i,])
  tmp <- lm(unlist(perm) ~ ind)
  a = summary(tmp)
  var.ind.adjr.p[i] <- a$adj.r.squared
}

var.ind.adjr = matrix(ncol = 1, nrow = dim(meth.final)[1])

for (i in 1:dim(meth.final)[1]){
  #for(i in 1:100){
  tmp <- lm(unlist(meth.final[,grep ("LCL|Fib" , colnames(meth.final), invert = T)][i,]) ~ ind)
  a = summary(tmp)
  var.ind.adjr[i] <- a$adj.r.squared
}

var.or.adjr = matrix(ncol = 1, nrow = dim(meth.final)[1])

for (i in 1:dim(meth.final)[1]){
  #for(i in 1:100){
  tmp <- lm(unlist(meth.final[,grep ("LCL|Fib" , colnames(meth.final), invert = T)][i,]) ~ origin_type)
  a = summary(tmp)
  var.or.adjr[i] <- a$adj.r.squared
}

var.resid.adjr = matrix(ncol = 1, nrow = dim(meth.final)[1])
for (i in 1:dim(meth.final)[1]){
   var.resid.adjr[i] <- 1 - var.or.adjr[i] - var.ind.adjr[i]
}

var.in.or.adjr = cbind(var.ind.adjr,var.or.adjr)
boxplot(var.in.or.adjr)
ind.median.adj <- c(median(var.in.or.adjr[var.in.or.adjr[,3] < .1 ,1]), median(var.in.or.adjr[var.in.or.adjr[,3] < .2 , 1]), median(var.in.or.adjr[var.in.or.adjr[,3] < .3,  1]), median(var.in.or.adjr[var.in.or.adjr[,3] < .4 , 1]), median(var.in.or.adjr[var.in.or.adjr[,3] < .5 ,1]), median(var.in.or.adjr[var.in.or.adjr[,3] < .6 , 1]), median(var.in.or.adjr[var.in.or.adjr[,3] < .7 , 1]), median(var.in.or.adjr[var.in.or.adjr[,3] < .8 , 1]), median(var.in.or.adjr[var.in.or.adjr[,3] < .9 , 1]), median(var.in.or.adjr[,1]))
ori.median.adj <- c(median(var.in.or.adjr[var.in.or.adjr[,3] < .1 ,2]), median(var.in.or.adjr[var.in.or.adjr[,3] < .2 , 2]), median(var.in.or.adjr[var.in.or.adjr[,3] < .3,  2]), median(var.in.or.adjr[var.in.or.adjr[,3] < .4 , 2]), median(var.in.or.adjr[var.in.or.adjr[,3] < .5 ,2]), median(var.in.or.adjr[var.in.or.adjr[,3] < .6 , 2]), median(var.in.or.adjr[var.in.or.adjr[,3] < .7 , 2]), median(var.in.or.adjr[var.in.or.adjr[,3] < .8 , 2]), median(var.in.or.adjr[var.in.or.adjr[,3] < .9 , 2]), median(var.in.or.adjr[,2]))
plot(ind.median.adj, pch = 16, ylim = c(0,1), xaxt= "n", main = "Proportion of Variance Explained by Individual and Tissue of Origin", ylab = "Proportion of Variance Explained", xlab = "Proportion of Variance Not Explained by Individual or Tissue of Origin")
axis(1, at=1:10, labels = c("< 10%", "< 20%","< 30%","< 40%","< 50%","< 60%","< 70%","< 80%","< 90%","< 100%"), cex.axis = .75)
points(ori.median.adj, pch = 16, col = "red")
legend("topright", c("Individual", "Tissue of Origin"), col = c("black", "red"), pch = c(16,16))

#boxplot(var.in.or, main = "Proportion of Variance Explained", xlab = "ind, origin, residual")
#most.var.in.or = var.in.or[var.in.or[,3] < .5 , ]
#boxplot(most.var.in.or, main = "Error less than 50%", xlab = "ind, origin, residual")
#most.var.in.or = var.in.or[var.in.or[,3] < .4 , ]
#boxplot(most.var.in.or, main = "Error less than 40%", xlab = "ind, origin, residual")
#most.var.in.or = var.in.or[var.in.or[,3] < .3 , ]
#boxplot(most.var.in.or, main = "Error less than 30%", xlab = "ind, origin, residual")

#var.resid.err.mvp = matrix(ncol = 1, nrow = 1000)
#var.origin.mvp = matrix(ncol = 1, nrow = 1000)
#var.ind.mvp = matrix(ncol = 1, nrow = 1000)
#for (i in 1:dim(most.variable.iPSC.probes)[1]){
#for(i in 1:1000){
#  tmp <- lm(unlist(most.variable.iPSC.probes[i,]) ~ ind + origin_type)
#  var.ind.mvp[i] <- anova(tmp)[1,2]/sum(anova(tmp)[,2])
#  var.origin.mvp[i] <- anova(tmp)[2,2]/sum(anova(tmp)[,2])
#  var.resid.err.mvp[i] <- anova(tmp)[3,2]/sum(anova(tmp)[,2])
#}

#boxplot(cbind(var.ind.mvp, var.origin.mvp, var.resid.err.mvp), main = "1,000 most variable probes", xlab = "ind, origin, residual")
#dev.off()

#Create plots to show have variance changes
ind.median <- c(median(var.in.or[var.in.or[,3] < .1 ,1]), median(var.in.or[var.in.or[,3] < .2 , 1]), median(var.in.or[var.in.or[,3] < .3,  1]), median(var.in.or[var.in.or[,3] < .4 , 1]), median(var.in.or[var.in.or[,3] < .5 ,1]), median(var.in.or[var.in.or[,3] < .6 , 1]), median(var.in.or[var.in.or[,3] < .7 , 1]), median(var.in.or[var.in.or[,3] < .8 , 1]), median(var.in.or[var.in.or[,3] < .9 , 1]), median(var.in.or[,1]))
ori.median <- c(median(var.in.or[var.in.or[,3] < .1 ,2]), median(var.in.or[var.in.or[,3] < .2 , 2]), median(var.in.or[var.in.or[,3] < .3,  2]), median(var.in.or[var.in.or[,3] < .4 , 2]), median(var.in.or[var.in.or[,3] < .5 ,2]), median(var.in.or[var.in.or[,3] < .6 , 2]), median(var.in.or[var.in.or[,3] < .7 , 2]), median(var.in.or[var.in.or[,3] < .8 , 2]), median(var.in.or[var.in.or[,3] < .9 , 2]), median(var.in.or[,2]))
ind.list <- list(var.in.or[var.in.or[,3] < .1 ,1], var.in.or[var.in.or[,3] < .2 , 1], var.in.or[var.in.or[,3] < .3,  1], var.in.or[var.in.or[,3] < .4 , 1], var.in.or[var.in.or[,3] < .5 ,1], var.in.or[var.in.or[,3] < .6 , 1], var.in.or[var.in.or[,3] < .7 , 1], var.in.or[var.in.or[,3] < .8 , 1], var.in.or[var.in.or[,3] < .9 , 1], var.in.or[,1])
ori.list <- list(var.in.or[var.in.or[,3] < .1 ,2], var.in.or[var.in.or[,3] < .2 , 2], var.in.or[var.in.or[,3] < .3,  2], var.in.or[var.in.or[,3] < .4 , 2], var.in.or[var.in.or[,3] < .5 ,2], var.in.or[var.in.or[,3] < .6 , 2], var.in.or[var.in.or[,3] < .7 , 2], var.in.or[var.in.or[,3] < .8 , 2], var.in.or[var.in.or[,3] < .9 , 2], var.in.or[,2])
#My samples contain no genes where <10% of variance is explained by individual or tissue of origin

pdf(file = "Proportion of Variance_Gene Expression.pdf")

boxplot(var.in.or, xaxt= "n", main = "Proportion of Gene Expression Variance", ylab = "Proportion of Variance Explained")
boxplot(ylim = c(0,1),ind.list, xaxt= "n",main = "Proportion of Gene Expression Variance Explained by Individual", ylab = "Proportion of Variance Explained by Individual", xlab = "Proportion of Variance Not Explained by Individual or Tissue of Origin")
axis(1, at=1:10, labels = c("<10%","< 20%","< 30%","< 40%","< 50%","< 60%","< 70%","< 80%","< 90%","< 100%"), cex.axis = .75)
#axis(side = 2, labels = c("10%", "20%","30%","40%","50%","60%","70%","80%","90%","100%"), cex.axis = .75, lwd = 0)

boxplot(ori.list, xaxt= "n", main = "Proportion of Gene Expression Variance Explained by Tissue of Origin", ylab = "Proportion of Variance Explained by Tissue of Origin", xlab = "Proportion of Variance Not Explained by Individual or Tissue of Origin")
#axis(side = 2, at=1:10, labels = c("10%", "20%","30%","40%","50%","60%","70%","80%","90%","100%"), cex.axis = .75, lwd = 0)

plot(ind.median, pch = 16, ylim = c(0,1), xaxt= "n", main = "Proportion of Variance Explained by Individual and Tissue of Origin", ylab = "Proportion of Variance Explained", xlab = "Proportion of Variance Not Explained by Individual or Tissue of Origin")
axis(1, at=1:10, labels = c("< 10%", "< 20%","< 30%","< 40%","< 50%","< 60%","< 70%","< 80%","< 90%","< 100%"), cex.axis = .75)
points(ori.median, pch = 16, col = "red")
legend("topright", c("Individual", "Tissue of Origin"), col = c("black", "red"), pch = c(16,16))
dev.off()


##Get a gene coord file
gene_coords = matrix(NA, ncol=8, nrow=length(expr_quant.all))

for(i in 1:dim(expr_quant.all)[1]){
  selection = which(goodprobes[,4]==row.names(expr_quant.all)[i])
  gene_coords[i,]= as.matrix(goodprobes[selection,])
}

gene_map = matrix(NA, ncol=8, nrow=length(unique(gene_names)))
i=0
for(gene in unique(gene_names)){
  i = i+1
  
  currRows = which(gene_names == gene)
  if(length(currRows)>1){
    if(goodprobes[currRows[1],6]=="+"){
      keepRow = currRows[which.max(goodprobes[currRows,2])]
    }
    else{
      keepRow = currRows[which.min(goodprobes[currRows,2])]
    }
  }
  else{
    keepRow=currRows[1]
  }
  gene_map[i,] = gene_coords[keepRow,1:8]
  
} 
colnames(gene_map) = colnames(goodprobes)
chr.gene = gene_map[,c(8,1)]

write.table(gene_map, 'OriginGeneCoords.txt', quote=F, sep ='\t', row.names=F)

