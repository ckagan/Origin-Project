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
colnames(detection) = samplenames$Name
type = samplenames$Type
selected = c('1')
iPSCdetect = detection[,type %in% selected]
selected = c('2')
LCLdetect = detection[,type %in% selected]
selected = c('3')
Fibdetect = detection[,type %in% selected]

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

###Subset expression by Darren good probes
goodprobes= read.table('HT-12v4_Probes_inhg19EnsemblGenes_NoCEUHapMapSNPs_Stranded.txt', header=T)
probes = goodprobes$probeID
## Convert from factor to character
probes = as.character(goodprobes$probeID)
expr_quant.all.clean = expr_quant.all[rownames(expr_quant.all) %in% probes, ]
dim(expr_quant.all.clean)  
# 16,593 probes of the original 26,953 "good" probes are in this data set
expr_quant.all= expr_quant.all.clean
remove(expr_quant.all.clean)
probelist = rownames(expr_quant.all)

##Load in Covariates
#Converted some factors so they are categorical 
samplenames = read.table('covar.txt', header=T, sep ='\t')
colnames(expr_quant.all) = samplenames$Name

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

#To use Nick's dendogram script
#First generate the chr file
probeinfolist = cbind(goodprobes$chr, as.character(goodprobes[,4]))
probelist = rownames(expr_quant.all)
probelist = as.matrix(probelist)
colnames(probelist) = c("ILMN")
colnames(probeinfolist) = c("Chr", "Probe")
chrlist = merge(probelist, probeinfolist, by.x = "ILMN", by.y = "Probe", all.x = T, all.y = F, sort=F)
chrfinal = as.matrix(chrlist[,2])

#If you want mean subtracted and variance divided data
#stan = apply(expr_quant.all,1, function(x) (x-mean(x))/sd(x))
stan = apply(expr_quant.all,1, function(x) x-mean(x))
stand = t(stan)
#Read in data
avg_beta = expr_quant.all
avg_beta = stand
chr = chrfinal

hist(as.numeric(chrfinal))
#Open pdf

pdf(file = "heatmap_dendrogram.pdf")

# Make dendogram of the data
plot(hclust(dist(t(avg_beta[,1:24]))), xlab = "", main = "Dendrogram")

#Make dendograms using pearson
plot(hclust(as.dist(1-cor(as.matrix(avg_beta)))))

# Make dendogram of the data without the X chr
plot(hclust(dist(t(avg_beta[chr != "24" ,1:24]))), xlab = "", main = "Dendrogram w/o X chr")

# Make dendogram of iPSCs only
plot(hclust(dist(t(avg_beta[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]))), xlab = "", main = "Dendrogram with only iPSCs")

#Pearson
plot(hclust(as.dist(1-cor(avg_beta[,grep ("LCL|Fib" , colnames(avg_beta), invert = T)]))), xlab = "", main = "Dendrogram with only iPSCs")

# Make dendogram of iPSCs only w/o X chr
plot(hclust(dist(t(avg_beta[chr != "24" ,grep ("LCL|Fib" , colnames(avg_beta), invert = T)][,1:16]))), xlab = "", main = "Dendrogram with only iPSCs, w/o X chr")


# Make heatmap of the data
heatmap(cor(as.matrix(avg_beta[,1:24]), use = "complete"))

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

#Make sample covar file with only iPSCs
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

samplenames.ipsc
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
colnames(pcaresults) = c("PC1 pval", "PC1 adj R sqs","P2 pval", "PC2 adj R sqs","PC3 pval", "PC3 adj R sqs","PC4 pval", "PC4 adj R sqs")

pcaresults.ipsc = lmPCA(ipsc.pca,covars.ipsc,4)
rownames(pcaresults.ipsc) = c("mef","batch","sex","indiv","pluri","novel","der")
colnames(pcaresults.ipsc) = c("PC1 pval", "PC1 adj R sqs","P2 pval", "PC2 adj R sqs","PC3 pval", "PC3 adj R sqs","PC4 pval", "PC4 adj R sqs")


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



#Make a heatmap with raw data
library("gplots")
cor <- cor(expr_quant.all,method="pearson", use="complete.obs")
heatmap.2(cor, main="Probe expression correlation", key=T, revC=T, density.info="none", trace="none")

#Make PC plots with raw data
library(TeachingDemos)
op <- par(mfrow = c(3,3), ## split region
          oma = c(5,0,4,0) + 0.1, ## create outer margin
          mar = c(5,4,2,2) + 0.1) ## shrink some margins
tmp1 <- cnvrt.coords( 0.5, 0, input='plt' )$tdev
cor <- cor(expr_quant.all, method="pearson", use="complete.obs")
title.PC = "PCA of Gene Exp"
sum.PC <- prcomp(na.omit(cor))
sumsum <- summary(sum.PC)
#prints out plots in c(#rows, #columns)
color = samplenames[,5]
par(mfrow = c(2,2),oma=c(0,0,2,0)) 
plot(c(1:12),sum.PC$rotation[,1],col=color, xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:12),sum.PC$rotation[,1], samplenames[,1], cex = .5, pos=3)   
for(i in 2:4) {
  plot(sum.PC$rotation[,1], sum.PC$rotation[,i], col=color,pch=20,main=title.PC, xlab=paste("PC 1 -", (sumsum$importance[2,1]*100),"% of variance", sep=" "), ylab=paste("PC",i,"-",(sumsum$importance[2,i]*100),"% of variance", sep=" "))
  }  
par(op)
#All plotting is clipped to the device region
par(xpd=NA)

##To see if covariates are correlated with a PC (looking at PC1-7)
sum.PC <- prcomp(na.omit(cor))
summary(lm(sum.PC$rotation[, 1:7] ~ age.num))
summary(lm(sum.PC$rotation[, 1:7] ~ cond.f))
summary(lm(sum.PC$rotation[, 1:7] ~ novel.num))
summary(lm(sum.PC$rotation[, 1:7] ~ sex.f))
summary(lm(sum.PC$rotation[, 1:7] ~ spec.f))
summary(lm(sum.PC$rotation[, 1:7] ~ tech.f))
summary(lm(sum.PC$rotation[, 1:7] ~ pluri.num))

#To get the correlation between replicates vs non-replicates
replicates=c()
non.replicates=c()
names = as.character(samplenames[1:11,12])
for(i in 1:10){
  for(j in (i+1):11){
    if(names[i]==names[j]){
      replicates=c(replicates,cor[i,j])
    }
    else{
      non.replicates=c(non.replicates,cor[i,j])
    }
  }
}
boxplot.n(replicates,non.replicates, main = "Correlation of iPSC Samples", ylab = "Correlation", xlab = "Replicates                                              Non-Replicates")


## Finding the unique gene names matching probes to gene names using Darren's good probe list
gene_names=c()
for(i in 1:dim(expr_quant.all)[1]){
  gene_names=c(gene_names,as.vector(goodprobes[as.vector(goodprobes[,4])==row.names(expr_quant.all)[i],8]))
}

symbolsUniq = unique(gene_names)
length(symbolsUniq)
#[1] 13698


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

expr_gene_ord = expr_gene[order(colnames(expr_gene))]


## To make heatmap (with key and title)
library(gplots)
cor.q <- cor(expr_gene,method="pearson", use="complete.obs")
heatmap.2(cor.q, main="Gene Expression", key=T, revC=T, density.info="none", trace="none")

## To make another version of heatmap with color codes 
library(gplots)
covarmi= as.character(samplenames[-12,7])
heatmap.2(cor.q, main= "Correlation", key=T, revC=T,ColSideColors=covarmi, density.info="none", trace="none")

##Make a heatmap using only iPSC
expr_gene_iPSC = expr_gene[,-12]
cor.i <- cor(expr_gene_iPSC,method="pearson", use="complete.obs")
heatmap.2(cor.i, main= "Correlation of iPSCs", key=T, revC=T,ColSideColors=covarmi, density.info="none", trace="none")

##PCA all probes
plot_colors<-c("black","purple")
library(TeachingDemos)
op <- par(mfrow = c(3,3), ## split region
          oma = c(5,0,4,0) + 0.1, ## create outer margin
          mar = c(5,4,2,2) + 0.1) ## shrink some margins
tmp1 <- cnvrt.coords( 0.5, 0, input='plt' )$tdev
title.PC = "PCA of Gene Exp"
sum.PC.i <- prcomp(na.omit(cor.i))
sumsum.i <- summary(sum.PC.i)
color = covarmi
#prints out plots in c(#rows, #columns)
par(mfrow = c(2,2),oma=c(0,0,2,0)) 
plot(c(1:11),sum.PC.i$rotation[,1],xlab="Index of Samples",pch = 20, ylab=paste("PC 1 -",(sumsum.i$importance[2,1]*100),"% of variance",sep=" "),main=title.PC)
text(c(1:11),sum.PC.i$rotation[,1], samplenames[-12,1], cex = .5, pos=3)   
for(i in 2:4) {
  plot(sum.PC.i$rotation[,1], sum.PC.i$rotation[,i], pch=20,main=title.PC, xlab=paste("PC 1 -", (sumsum.i$importance[2,1]*100),"% of variance", sep=" "), ylab=paste("PC",i,"-",(sumsum.i$importance[2,i]*100),"% of variance", sep=" "))
  text(sum.PC.i$rotation[,1], sum.PC.i$rotation[,i], samplenames[-12,1], cex = .5, pos=3)}  
par(op)


##Extraction PC data from PCA
PCS1 = sum.PC$rotation
View(PCS1)
PCS2 =sum.PC2$rotation
View(PCS2)
PCS3 =sum.PC3$rotation
View(PCS3)
#Check this
#PCs are columns, inds are rows, like this:
#               PC1          PC2          PC3           PC4          PC5
#87_MH  -0.08176976 -0.132745281 -0.069088998  0.1525655116 -0.106039703
#125_MH -0.18216685 -0.083451452  0.079887299 -0.1936390876 -0.204922215
#168_ML -0.01938816 -0.205549632 -0.027594568 -0.0324146775 -0.223454860
#209_FL -0.12189288 -0.129940372 -0.052735639  0.0120791405 -0.178501057


#Looking for corellation of PC's to covariates
#Look at p-value
summary(lm(sum.PC.final$rotation[1:9, 1:7] ~ prot.f))
summary(lm(sum.PC.final$rotation[-2, 1:7] ~ bmi))

#To pull out expression information
#For the expression of CDKL3 gene
expr_gene["MESP1",]
finalstatus = samplenames[,13]
#To make a boxplot
boxplot(expr_gene["MYL3",]~finalstatus, main = "Expression of MYL3", ylab = "expression")