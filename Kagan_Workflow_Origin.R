##Load lumi
library(lumi)

##Create object with name of data file:
data = c('YGilad-CK-Aug23-14-ProbeLevelData-NotNormalized-NoBGSubtracted-FinalReport.txt')

##Extract raw data from lumi file and preparing for removing bad probes:
data.lumi = lumiR.batch(data, lib.mapping=NULL, convertNuID=F,annotationColumn=c('ACCESSION', 'SYMBOL', 'PROBE_SEQUENCE', 'PROBE_START', 'CHROMOSOME', 'PROBE_CHR_ORIENTATION', 'PROBE_COORDINATES', 'DEFINITION','PROBE_ID'))


### NORMALIZATION: log2 stabilized and quantile normalization ###
data.norm.all <- lumiExpresso(data.lumi, bg.correct=TRUE, bgcorrect.param=list(method='forcePositive'), variance.stabilize=TRUE, varianceStabilize.param = list(method="log2"), normalize=TRUE, normalize.param=list(method="quantile"), QC.evaluation=TRUE, QC.param=list(), verbose=TRUE)
#Take only the probes that have a detection p-value<.05 in at least one individual
detect_quant.all= rowSums(data.norm.all@assayData$detection<0.05)
norm_quant.all <- data.norm.all@assayData$exprs
detect.ind.all <- which(detect_quant.all > 2)



##Look at plots of array data (boxplot, density, etc) :
plot(data.lumi, what='boxplot')
plot(data.lumi, what='density')
plot(data.norm.all, what='boxplot')
plot(data.norm.all, what='density')

##Check that replicates are most related
plot(data.norm.all, what='sampleRelation')

###Find the column that is lumi_ID in feature data usually this column
head(data.norm.all@featureData[[5]])
##[1] "ILMN_1343291" "ILMN_1343295" "ILMN_1651199" "ILMN_1651209" "ILMN_1651210" "ILMN_1651221"
###Convert expr_quant rownames to
rownames(norm_quant.all)=data.norm.all@featureData[[5]]
#With this threshold 26,953 probes out of 47,315 are expressed
expr_quant.all <- norm_quant.all[detect.ind.all,]


###Subset expression by Darren good probes
goodprobes= read.table('HT-12v4_Probes_inhg19EnsemblGenes_NoCEUHapMapSNPs_Stranded.txt', header=T)
probes = goodprobes$probeID
## Convert from factor to character
probes = as.character(goodprobes$probeID)
expr_quant.all.clean = expr_quant.all[rownames(expr_quant.all) %in% probes, ]
# dim(expr_quant.all.clean)  
# 17,867 probes of the original 26,953 "good" probes are in this data set
expr_quant.all= expr_quant.all.clean
remove(expr_quant.all.clean)

##Load in Covariates
#Converted chip and batch to a factor so they are categorical 
samplenames = read.table('sample_names_corrected.txt', header=F, sep ="")
colnames(expr_quant.all) = samplenames[,1]

cond = samplenames[,2]
spec = samplenames[,4]
tech = samplenames[,6]
sex =  samplenames[,8]
age = samplenames[,9]
pluri =  samplenames[,10]
novel = samplenames[,11]
#Converted categorical covariates to a factor so they are levels 
cond.f = as.factor(cond)
spec.f = as.factor(spec)
tech.f = as.factor(tech)
sex.f = as.factor(sex)
##To look for correlations between covariates get the R^2 of the Pearson correlation
cor(age, pluri)^2
cor(age, novel)^2
cor(sex, pluri)^2
cor(sex, novel)^2
cor(spec, pluri)^2
cor(spec, novel)^2
cor(tech, pluri)^2
cor(tech, novel)^2
#Converted numerical covariates to a numeric so they are continuous
age.num = as.numeric(age)
pluri.num =as.numeric(pluri)
novel.num = as.numeric(novel)
rm(cond, spec, tech, sex, age, pluri, novel)

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