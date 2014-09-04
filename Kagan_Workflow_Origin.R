##Load lumi
library(lumi)
setwd("C:/Users/Courtney/Dropbox/LCL-iPSC/Origin/Array")

##Create object with name of data file:
data = c('YGilad-CK-Aug23-14-ProbeLevelData-NotNormalized-NoBGSubtracted-FinalReport.txt')

##Extract raw data from lumi file and preparing for removing bad probes:
data.lumi = lumiR.batch(data, lib.mapping=NULL, convertNuID=F,annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION','PROBE_ID'))

#Get QC Data
summary(data.lumi, 'QC')

#Past output into excel and create a file called lumiQC
#Sample  mean standard.deviation detection.rate.0.01. distance.to.sample.mean
#1    YG1 7.229              1.211               0.3735                   42.24
#2    YG2 7.342              1.257               0.3543                   46.28

#Look at QC based on cell type
qcdata = read.table('lumiQC.txt', header=T, as.is=T, sep='\t')
samplenames = read.table('covar.txt', header=T, sep ='\t')

boxplot(qcdata$mean~samplenames$Indiv, main = 'Mean Probe Intensity by Individual')
boxplot(qcdata$mean~samplenames$Batch, main = 'Mean Probe Intensity by Array')
boxplot(qcdata$mean~samplenames$MEF, main = 'Mean Probe Intensity by MEF Batch')
boxplot(qcdata$mean~samplenames$Type, main = 'Mean Probe Intensity by Cell Type')

### NORMALIZATION: log2 stabilized and quantile normalization ###
data.norm.all <- lumiExpresso(data.lumi, bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)

##Normalize by cell type
data.norm.ipsc = lumiExpresso(data.lumi[, ipscsubset], bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)
summary(data.norm.ipsc, 'QC')
data.norm.fib = lumiExpresso(data.lumi[, fibsubset], bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)
summary(data.norm.fib, 'QC')
data.norm.lcl = lumiExpresso(data.lumi[, lclsubset], bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)
summary(data.norm.lcl, 'QC')
data.norm.all = combine(data.norm.ipsc,data.norm.fib,data.norm.lcl)
sampleorder = c(as.character(samplenames[ipscsubset,3]),as.character(samplenames[fibsubset,3]),as.character(samplenames[lclsubset,3]))
colnames(data.norm.all) = sampleorder

show(data.norm.all)
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

#Subset data by cell type and look for detection in at least 4 indivduals per subset
detection = data.norm.all@assayData$detection
#colnames(detection) = sampleorder
colnames(detection) = samplenames$Name

iPSCdetect = detection[,type %in% selected]
iPSCdetect = detection[,grep ("LCL|Fib", colnames(detection), invert=T)]
LCLdetect = detection[,grep ("LCL", colnames(detection))]
Fibdetect = detection[,grep ("Fib", colnames(detection))]

#Count the number of probes that have a detection p-value<.05
iPSCdetected = rowSums(iPSCdetect<0.05)
LCLdetected = rowSums(LCLdetect<0.05)
Fibdetected = rowSums(Fibdetect<0.05)

#Here is where I select the number of indiv that need to have the probe expressed (at least in 4 people), iPSC =22,848, LCL = 15,759, Fib = 16,699
detect.iPSC <- which(iPSCdetected > 3)
detect.LCL <- which(LCLdetected > 3)
detect.Fib <- which(Fibdetected > 3)

union = c(detect.Fib, detect.LCL, detect.iPSC)
unionunique= unique(union)
detect.ind.all = sort.int(unionunique)
#24,480 probes detected 

##Look at plots of array data (boxplot, density, etc) :
plot(data.lumi, what='boxplot')
plot(data.lumi, what='density')
plot(data.norm.all, what='boxplot')
plot(data.norm.all, what='density')

##Check that replicates are most related
plot(data.norm.all, what='sampleRelation')

norm_quant.all <- data.norm.all@assayData$exprs
###Find the column that is lumi_ID in feature data usually this column
head(data.norm.all@featureData[[5]])
##[1] "ILMN_1343291" "ILMN_1343295" "ILMN_1651199" "ILMN_1651209" "ILMN_1651210" "ILMN_1651221"
###Convert expr_quant rownames to
rownames(norm_quant.all)=data.norm.all@featureData[[5]]
#With this threshold 24,480 probes out of 47,315 are expressed
expr_quant.all <- norm_quant.all[detect.ind.all,]
#expr_iPSC_quant.all = norm_quant.all[detect.iPSC,]
expr_quant.all = expr_iPSC_quant.all

###Subset expression by Darren's good probes
goodprobes= read.table('HT-12v4_Probes_inhg19EnsemblGenes_NoCEUHapMapSNPs_Stranded.txt', header=T)
probes = goodprobes$probeID
## Convert from factor to character
probes = as.character(goodprobes$probeID)
expr_quant.all.clean = expr_quant.all[rownames(expr_quant.all) %in% probes, ]
expr_iPSC_quant.all.clean = expr_iPSC_quant.all[rownames(expr_iPSC_quant.all) %in% probes, ]
dim(expr_quant.all.clean)  
# 16,593 probes of the original 26,953 "good" probes are in this data set
expr_quant.all= expr_quant.all.clean
remove(expr_quant.all.clean)
probelist = rownames(expr_quant.all)
probelist_iPSC = rownames(expr_iPSC_quant.all.clean)

##Load in Covariates
#Converted some factors so they are categorical 
samplenames = read.table('covar.txt', header=T, sep ='\t')
colnames(expr_quant.all) = samplenames$Name
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

## Finding the unique gene names matching probes to gene names using Darren's good probe list
gene_names=c()
for(i in 1:dim(expr_quant.all)[1]){
  gene_names=c(gene_names,as.vector(goodprobes[as.vector(goodprobes[,4])==row.names(expr_quant.all)[i],8]))
}

symbolsUniq = unique(gene_names)
length(symbolsUniq)
#[1] 12872


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

#write.table(gene_map, 'OriginGeneCoords.txt', quote=F, sep ='\t', row.names=F)


####Data Analysis####

#To use Nick's dendogram script
#First generate the chr file
probeinfolist = cbind(goodprobes$chr, as.character(goodprobes[,4]))
probelist = rownames(expr_quant.all)
probelist = as.matrix(probelist)
colnames(probelist) = c("ILMN")
colnames(probeinfolist) = c("Chr", "Probe")
chrlist = merge(probelist, probeinfolist, by.x = "ILMN", by.y = "Probe", all.x = T, all.y = F, sort=F)
chrfinal = as.matrix(chrlist[,2])
hist(as.numeric(chrfinal))

#If you want mean subtracted and variance divided data
stan = apply(expr_quant.all,1, function(x) x-mean(x))
stand = t(stan)

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
# 0.03724039
mean(variance.iPSC)
# 0.02375778
mean(variance.Fib)
# 0.03143111

variance.all = cbind(variance.LCL, variance.Fib, variance.iPSC)
boxplot(variance.all, ylim = c(-.01,.05))
var.test(variance.LCL, variance.iPSC)

#F test to compare two variances

#data:  variance.LCL and variance.iPSC
#F = 3.581, num df = 12871, denom df = 12871, p-value < 2.2e-16
#alternative hypothesis: true ratio of variances is not equal to 1
#95 percent confidence interval:
#  3.459422 3.706946
#sample estimates:
#  ratio of variances 
#3.581046 

library(ggplot2)
var_all <- data.frame(var=c(variance.LCL, variance.Fib, variance.iPSC), type = rep(c("LCL","Fib", "iPSC"), times=c(length(variance.LCL))))
ggplot(var_all, aes(x=var, fill=type)) + geom_density(alpha=0.1) +xlim(-.01,.2)+xlab("Variance") + ggtitle("Gene Expression: Total variance") + theme(legend.position=c(.75,.75)) +theme(text = element_text(size=23))

var_LCL <- data.frame(var=c(variance.LCL), type = rep(c("LCL"), times=c(length(variance.LCL))))
ggplot(var_LCL, aes(x=var, fill=type)) + geom_density(alpha=0.5) +xlim(-.01,.2)+xlab("Variance") + ggtitle("Gene Expression: Total variance") + theme(legend.position=c(.75,.75)) +theme(text = element_text(size=23))
var_iPSC <- data.frame(var=c(variance.iPSC), type = rep(c("iPSC"), times=c(length(variance.iPSC))))
ggplot(var_iPSC,aes(x=var, fill=type))+ geom_density(alpha=0.5) +xlim(-.01,.2)+xlab("Variance") + ggtitle("Gene Expression: Total variance") + theme(legend.position=c(.75,.75)) +theme(text = element_text(size=23))

#Read in data
avg_beta = expr_gene
#avg_beta = expr_quant.all
#avg_beta = stand
#chr = chrfinal
chr = chr.gene

#Open pdf

pdf(file = "heatmap_dendrogram.pdf")

# Make dendogram of the data
plot(hclust(dist(t(avg_beta[,1:24]))), xlab = "hclust/Euclidean distance", main = "Dendrogram")

plot(hclust(dist(t(avg_beta[,1:24]), method = "manhattan")), xlab = "Manhattan", main = "Dendrogram")
plot(hclust(dist(t(avg_beta[,1:24]), method = "canberra")), xlab = "Canberra", main = "Dendrogram")

#Make dendograms using pearson
plot(hclust(as.dist(1-cor(as.matrix(avg_beta)))),xlab = "Pearson", main = "Dendrogram using Pearson Correlation")

# Make dendogram of the data without the X chr
plot(hclust(dist(t(avg_beta[chr[,2] != "chrX" ,]))), xlab = "hclust/Euclidean distance", main = "Dendrogram w/o X chr")

# Make dendogram of iPSCs only
plot(hclust(dist(t(avg_beta[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]))), xlab = "", main = "Dendrogram with only iPSCs")

#Pearson
plot(hclust(as.dist(1-cor(avg_beta[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]))), xlab = "", main = "Dendrogram with only iPSCs")

# Make dendogram of iPSCs only w/o X chr
plot(hclust(dist(t(avg_beta[chr[,2] != "chrX" ,grep ("LCL|Fib" , colnames(avg_beta), invert = T)][,1:16]))), xlab = "", main = "Dendrogram with only iPSCs, w/o X chr")


# Make heatmap of the data
heatmap(cor(as.matrix(avg_beta[,1:24]), use = "complete"),main="Gene expression correlation",)

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
plot(x.pca$rotation[,2], x.pca$rotation[,3], xlab = paste('PC2 (',x.pca.sum$importance[2,2], ')', sep = ''),ylab = paste('PC3 (',x.pca.sum$importance[2,3], ')', sep = ''), main = "PC2/3 all data", col = Leg, pch = 20); legend(x = "topleft", pch = 20, col = c(1:4), c("LCL origin", "Fib origin", "LCL", "Fib"))

ipsc.pca = prcomp(na.omit(avg_beta[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]), scale = T, center =T)
ipsc.pca.sum = summary(ipsc.pca)

dev.off()
#Make sample covar file with only iPSCs
pdf(file = "PCA iPSC only.pdf")
rem = grep ("LCL|Fib" , samplenames$Name)
samplenames.ipsc = samplenames[-rem,]

ipsc_only_leg = c("red", "blue", "orange", "black", "blue", "black", "red", "orange", "blue", "red", "blue", "orange", "red", "black", "orange", "black")

ipsc_only_shape =  c(18, 20,20,20,20,20,20,18, 20,20,18,20,20,20,20,18)

plot(ipsc.pca$rotation[,1], ipsc.pca$rotation[,2], xlab = paste('PC1 (',ipsc.pca.sum$importance[2,1], ')', sep = ''),ylab = paste('PC2 (',ipsc.pca.sum$importance[2,2], ')', sep = ''), main = "PC1/2 iPSC only", col = ipsc_only_leg, pch = ipsc_only_shape); legend(x = "topright", pch = c(20, 20, 20, 20, 20, 18), col = c("red","blue","black","orange", "black","black"), c("0961", "1194", "4280", "8126", "LCL derived", "Fib derived"))


plot(ipsc.pca$rotation[,2], ipsc.pca$rotation[,3], xlab = paste('PC2 (',ipsc.pca.sum$importance[2,2], ')', sep = ''),ylab = paste('PC3 (',ipsc.pca.sum$importance[2,3], ')', sep = ''), main = "PC2/3 iPSC only", col = ipsc_only_leg, pch = ipsc_only_shape); legend(x = "topright", pch = c(20, 20, 20, 20, 20, 18), col = c("red","blue","black","orange", "black","black"), c("0961", "1194", "4280", "8126", "LCL derived", "Fib derived"))


plot(ipsc.pca$rotation[,3], ipsc.pca$rotation[,4], xlab = paste('PC3 (',ipsc.pca.sum$importance[2,3], ')', sep = ''),ylab = paste('PC4 (',ipsc.pca.sum$importance[2,4], ')', sep = ''), main = "PC3/4 iPSC only", col = ipsc_only_leg, pch = ipsc_only_shape); legend(x = "topright", pch = c(20, 20, 20, 20, 20, 18), col = c("red","blue","black","orange", "black","black"), c("0961", "1194", "4280", "8126", "LCL derived", "Fib derived"))


plot(ipsc.pca$rotation[,4], ipsc.pca$rotation[,5], xlab = paste('PC4 (',ipsc.pca.sum$importance[2,4], ')', sep = ''),ylab = paste('PC5 (',ipsc.pca.sum$importance[2,5], ')', sep = ''), main = "PC4/5 iPSC only", col = ipsc_only_leg, pch = ipsc_only_shape); legend(x = "topright", pch = c(20, 20, 20, 20, 20, 18), col = c("red","blue","black","orange", "black","black"), c("0961", "1194", "4280", "8126", "LCL derived", "Fib derived"))

plot(ipsc.pca$rotation[,5], ipsc.pca$rotation[,6], xlab = paste('PC5 (',ipsc.pca.sum$importance[2,5], ')', sep = ''),ylab = paste('PC6 (',ipsc.pca.sum$importance[2,6], ')', sep = ''), main = "PC5/6 iPSC only", col = ipsc_only_leg, pch = ipsc_only_shape); legend(x = "topright", pch = c(20, 20, 20, 20, 20, 18), col = c("red","blue","black","orange", "black","black"), c("0961", "1194", "4280", "8126", "LCL derived", "Fib derived"))


dev.off()

covars = list(mef, batch, type,sex,indiv,pluri,novel,der)

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
rownames(pcaresults) = c("mef","batch", "type","sex","indiv","pluri","novel","der")
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
# Make heatmap of the data
heatmap(cor(as.matrix(meth.final), use = "complete"))

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

overlap_iPSC_DMR_Origin_DMRs = iPSC_DMR[iPSC_DMR$adj.P.Val < 0.01 , ][rownames(iPSC_DMR[iPSC_DMR$adj.P.Val < 0.01 ,]) %in% rownames(LCL_vs_Fibs[LCL_vs_Fibs$adj.P.Val < 0.01 ,]) , ]

dist_of_Ps = iPSC_DMR[rownames(iPSC_DMR) %in% rownames(LCL_vs_Fibs[LCL_vs_Fibs$adj.P.Val < 0.01 , ]) , ]

head(iPSC_DMR[iPSC_DMR$adj.P.Val < 0.01 , ])
dim(LCL_vs_Fibs[LCL_vs_Fibs$adj.P.Val < 0.01 , ])
dim(iPSC_DMR[iPSC_DMR$adj.P.Val < 0.01 , ])
dim(LCLs_vs_iPSC.L[LCLs_vs_iPSC.L$adj.P.Val < 0.01 , ])
dim(Fibs_vs_iPSC.F[Fibs_vs_iPSC.F$adj.P.Val < 0.01 , ])

write.table(LCL_vs_Fibs,'DE_LCLvFib.txt', quote=F, sep = '\t')
write.table(LCLs_vs_iPSC.L,'DE_LCLviPSC.txt', quote=F, sep = '\t')
write.table(Fibs_vs_iPSC.F,'DE_FibviPSC.txt', quote=F, sep = '\t')

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
library(VennDiagram, lib="~/R_libs")

make.venn.pair <- function(geneset1, geneset2, prefix, geneset1.label, geneset2.label,universe){
  pdf(file=paste(prefix, ".pdf", sep=""), width=7, height=7)
  venn.placeholder <- draw.pairwise.venn(length(geneset1),length(geneset2), length(which(geneset1 %in% geneset2) == TRUE), c(geneset1.label, geneset2.label), fill=c("goldenrod", "plum4"), alpha=c(0.5, 0.5),col=NA, euler.d=T)
  complement.size <- dim(universe)[1] - length(geneset1) - length(geneset2) + length(which(geneset1 %in% geneset2) == TRUE)
  grid.text(paste(complement.size, " not DE\nin either", sep=""), x=0.1, y=0.1)
  dev.off()
}

make.venn.pair(ran.array.ipsc[ran.array.ipsc$adj.P.Val < 0.01,]$genes, ran.array.ipsc[ran.array.ipsc$DE.liver > 10,]$genes, "de_ipsc_vs_liver_2008", "DE iPSC\n FDR 1%", "DE liver,\nBlekhman 2008\n FDR 0.6%", ran.array.ipsc)
make.venn.pair(ran.array.ipsc[ran.array.ipsc$adj.P.Val < 0.01,]$genes, ran.array.ipsc[ran.array.ipsc$DE.heart > 10,]$genes, "de_ipsc_vs_heart_2008", "DE iPSC\n FDR 1%", "DE heart,\nBlekhman 2008\n FDR 0.6%", ran.array.ipsc)
make.venn.pair(ran.array.ipsc[ran.array.ipsc$adj.P.Val < 0.01,]$genes, ran.array.ipsc[ran.array.ipsc$DE.kidney > 10,]$genes, "de_ipsc_vs_kidney_2008", "DE iPSC\n FDR 1%", "DE kidney,\nBlekhman 2008\n FDR 0.6%", ran.array.ipsc)

make.venn.pair(ran.exprs.ipsc[ran.exprs.ipsc$adj.P.Val < 0.01,]$genes, ran.exprs.ipsc[ran.exprs.ipsc$deHC == T,]$genes, "de_ipsc_vs_liver_2010", "DE iPSC\n FDR 1%", "DE liver,\nBlekhman 2010\n FDR 5%", ran.exprs.ipsc)



make.venn.triple <- function(geneset1, geneset2, geneset3, prefix, geneset1.label, geneset2.label, geneset3.label, universe){
  universe$g1 <- universe$genes %in% geneset1
  universe$g2 <- universe$genes %in% geneset2
  universe$g3 <- universe$genes %in% geneset3
  pdf(file=paste(prefix, ".pdf", sep=""), width=7, height=7)
  venn.placeholder <- draw.triple.venn(length(geneset1), length(geneset2), length(geneset3), dim(universe[universe$g1 == T & universe$g2 == T,])[1], dim(universe[universe$g2 == T & universe$g3 == T,])[1], dim(universe[universe$g1 == T & universe$g3 == T,])[1], dim(universe[universe$g1 == T & universe$g2 == T & universe$g3 == T,])[1], c(geneset1.label, geneset2.label, geneset3.label), fill=c("goldenrod", "plum4", "steelblue3"), alpha=c(0.5, 0.5, 0.5),col=NA, euler.d=T)
  complement.size <- dim(universe[universe$g1 == F & universe$g2 == F & universe$g3 == F,][1])
  grid.text(paste(complement.size, " not DE\nin any", sep=""), x=0.1, y=0.1)
  dev.off()
}


make.venn.triple(ran.array.ipsc[ran.array.ipsc$DE.kidney > 10,]$genes, ran.array.ipsc[ran.array.ipsc$DE.heart > 10,]$genes, ran.array.ipsc[ran.array.ipsc$DE.liver > 10,]$genes, "de_triple_blekhman_2008", "DE Kidney", "DE Heart", "DE Liver", ran.array.ipsc)