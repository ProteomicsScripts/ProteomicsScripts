library(mice)
library(GAMBoost)
library(DESeq)
library(openxlsx)
library(limma)
library(edgeR)
library(sva)
cols = c('steelblue1','red')
source("~/Dropbox/scripts/heatmap3.R")
source("~/Dropbox/scripts/biplot_custom.R")
source("~/Dropbox/scripts/var_components.R")
# comment
matcolor =  colorRampPalette(c('white','black'))(n = 20)
matcolorposneg =  colorRampPalette(c('blue','white','red'))(n = 50)
todaysdate = gsub("-","",Sys.Date())

## This is an R script for the conversion of mzTab to a better readable tsv format.

# clear entire workspace
rm(list = ls())
setwd("/home/hess/Studien/022_Schilling_Proteomics/Harvoni/WhitelistLFQ_Braunschweig")
input.file <- "E5023_HCV_HRW_H_1_0.mzTab"
whitelist.input.file <- "E5023_HCV_HRW_H_1_0__WHITELISTED.mzTab"
output.file <- "E5023_HCV_HRW_H_1_0.tsv"

# find start of the section
startSection <- function(file, section.identifier) {
  data <- file(file, "r")
  row = 0
  while (TRUE) {
    row = row + 1
    line = readLines(data, n=1)
    if (substr(line, 1, 3)==section.identifier) {
      break
    }
  }
  close(data)
  return (row)
}

# count the occurences of character c in string s
countOccurrences <- function(char,s) {
	s2 <- gsub(char,"",s)
	return (nchar(s) - nchar(s2))
}

# check that the protein accession is of the format *|*|*
# Note that NA returns TRUE.
checkAccessionFormat <- function(accession) {
  if (is.na(accession)) {
    return (TRUE)
  }
  n <- length(accession)
  count <- countOccurrences("[|]",accession)
  m <- length(which(count==2))
  return (n==m)
}

# Extracts the second entry from a string of the form *|*|*.
getAccession <- function(string) {
  if (is.na(string)) {
    return (NA)
  }
  return (unlist(strsplit(string, "[|]"))[2])
}

# Extracts the third entry from a string of the form *|*|*.
getGene <- function(string) {
  if (is.na(string)) {
    return (NA)
  }
  return (unlist(strsplit(string, "[|]"))[3])
}

# read the PEP section of an mzTab file
readMzTabPEP <- function(file) {
  
  # find start of the PEP section
  first.row <- startSection(file, "PEH")
  
  # read entire mzTab
  data <- read.table(file, sep="\t", skip=first.row-1, fill=TRUE, header=TRUE, quote="", na.strings=c("null","NA"), stringsAsFactors=FALSE, check.names=FALSE)
  
  # extract PEP data
  peptide.data <- data[which(data[,1]=="PEP"),]
  peptide.data$PEH <- NULL
  
  # In case the accession column is of the format *|*|*, we split this column into an accession and a gene column.
  if (all(sapply(peptide.data$accession, checkAccessionFormat))) {
    peptide.data$gene <- sapply(peptide.data$accession, getGene)
    peptide.data$accession <- sapply(peptide.data$accession, getAccession)
  }
  
  return (peptide.data)
}

peptide.data <- readMzTabPEP(input.file)
## write.table(peptide.data, output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

whitelist.peptide.data <- readMzTabPEP(whitelist.input.file)
wlpep = grepl("peptide_abundance_study_variable",colnames(whitelist.peptide.data))
wlpep = whitelist.peptide.data[,wlpep]



peptabund = grepl("peptide_abundance_study_variable",colnames(peptide.data))
peptdata = peptide.data[,peptabund]

nolabel = is.na(peptide.data[,1])

nonan <- rowSums(is.na(peptdata))
nonan <- nonan ==0
x = peptdata[nonan,]
setwd("/home/hess/Studien/022_Schilling_Proteomics/Harvoni/")

samples  <- readWorkbook("20180418_sample_list_lars.xlsx",sheet=1,startRow = 1,rowNames=TRUE)

assignments <- read.csv("20180418_labels.txt",stringsAsFactors=FALSE,sep="\t",header=FALSE)
assignments[,3] <- gsub("file:///","",assignments[,3],fixed=TRUE)
colnames(x) <- assignments[,3]

samples <- samples[colnames(x),]

group = samples$Group
batch = samples$Batch
batch2 = gsub(".*HCV_(.*?)_.*","\\1",rownames(samples)
xtabs(~batch2+batch)
xtabs(~group+batch2)
# b2 = gsub(".*HCV_(.*?)_.*","\\1",rownames(samples))
setwd("/home/hess/Studien/022_Schilling_Proteomics/Harvoni/out/")

###################################################################
########## Begin with Analysis ####################################
###################################################################
## image(abs(cor(peptdata[nonan,])))
## heatmap(cor(peptdata[nonan,]))


fname = sprintf("%s_distribution_per_sample_all_peptides.png",todaysdate)
png(fname,width=1500,height=2000,pointsize=24)
par(mfrow=c(3,1),mar=c(3,3,3,3))
boxplot(log10(peptdata+1),col=c('red','steelblue1','olivedrab3')[as.factor(group)])
boxplot(log10(peptdata+1),col=c('red','steelblue1','olivedrab3')[as.factor(batch)])
dev.off()


fname = sprintf("%s_distribution_per_sample_na_removed.png",todaysdate)
png(fname,width=1500,height=2000,pointsize=24)
par(mfrow=c(3,1),mar=c(3,3,3,3))
boxplot(log10(x+1),col=c('red','steelblue1','olivedrab3')[as.factor(group)])
boxplot(log10(x+1),col=c('red','steelblue1','olivedrab3')[as.factor(batch)])
dev.off()

fname = sprintf("%s_distribution_per_sample_whitelist.png",todaysdate)
png(fname,width=1500,height=2000,pointsize=24)
par(mfrow=c(3,1),mar=c(3,3,3,3))
boxplot(log10(wlpep+1),col=c('red','steelblue1','olivedrab3')[as.factor(group)])
boxplot(log10(wlpep+1),col=c('red','steelblue1','olivedrab3')[as.factor(batch)])
dev.off()

nf = calcNormFactors(x)
model <- model.matrix(~group)
fname = sprintf("%s_voom_transformation.png",todaysdate)
png(fname,width=1500,height=1500,pointsize=24)
datavoom = voom(x,model,lib.size=colSums(x)*nf,plot=TRUE)$E
dev.off()

fname = sprintf("%s_distribution_per_sample_all_peptides_voom_whitelist_correction_factors.png",todaysdate)
png(fname,width=1500,height=2000,pointsize=24)
par(mfrow=c(3,1),mar=c(3,3,3,3))
boxplot(datavoom,col=c('red','steelblue1','olivedrab3')[as.factor(group)])
boxplot(datavoom,col=c('red','steelblue1','olivedrab3')[as.factor(batch)])
dev.off()


mod = model.matrix(~group)
# boxplot(,col=c('red','steelblue1','olivedrab3')[as.factor(group)])
cb = ComBat(datavoom,batch=batch,mod=mod)

fname = sprintf("%s_distribution_per_sample_all_peptides_voom_combat_na_removed.png",todaysdate)
png(fname,width=1500,height=2000,pointsize=24)
par(mfrow=c(3,1),mar=c(3,3,3,3))
boxplot(cb,col=c('red','steelblue1','olivedrab3')[as.factor(group)])
boxplot(cb,col=c('red','steelblue1','olivedrab3')[as.factor(batch)])
dev.off()






pca = prcomp(scale(t(datavoom)),scale.=FALSE,center=FALSE)
fname = sprintf("%s_pca_voom_whitelist_correction_factors.pdf",todaysdate)
pdf(fname,width=12,height=12,pointsize=24)
par(mfrow=c(1,1),mar=c(3,3,3,3))
f.upper <- function(x,y)
{
    points(x,y,col=c('red','steelblue1','olivedrab3')[as.factor(group)],pch=19)
    text(x,y,b2,cex=.3)
}
f.lower <- function(x,y)
{
    text(x,y,b2,cex=.3)
}
pairs(pca$x[,1:5],upper.panel=f.upper,lower.panel=f.lower)
pairs(pca$x[,1:5],col=c('red','steelblue1','olivedrab3')[as.factor(batch)],pch=19)
dev.off()

pca = prcomp(scale(t(cb)),scale.=FALSE,center=FALSE)
fname = sprintf("%s_pca_voom_combat_whitelist_correction_factors.pdf",todaysdate)
pdf(fname,width=12,height=12,pointsize=24)
par(mfrow=c(1,1),mar=c(3,3,3,3))
f.upper <- function(x,y)
{
    points(x,y,col=c('red','steelblue1','olivedrab3')[as.factor(group)],pch=19)
    text(x,y,b2,cex=.3)
}
f.lower <- function(x,y)
{
    text(x,y,b2,cex=.3)
}
pairs(pca$x[,1:5],upper.panel=f.upper,lower.panel=f.lower)
pairs(pca$x[,1:5],col=c('red','steelblue1','olivedrab3')[as.factor(batch)],pch=19)
dev.off()



###################################################################
########## Variable Selection #####################################
###################################################################


## BOosting
## BL vs. fu48

set.seed(1234)
ybinbool = group!="H"
ybin = ifelse(group[ybinbool]=="BL",1,0)
xbin = scale(t(datavoom[,ybinbool]))
## bootstrap
k = 100
bootstrap <- lapply(1:k, function(arg) sample(1:nrow(xbin),replace=TRUE))
bsres <- matrix(0,ncol(xbin),k)
iter = 0
for (bo in bootstrap)
{
    iter = iter + 1
    xsub = xbin[bo,]
    ysub = ybin[bo]
    boostres <- cv.GLMBoost(xsub,ysub,penalty=100,maxstepno=100,trace=TRUE,family=binomial(),K=10,type="error")
    boostres <- GLMBoost(xsub,ysub,penalty=100,stepno=boostres$selected,family=binomial())
    
    selgenes <- unique(boostres$selected)
    bsres[selgenes,iter] = 1

}
inclfreq = rowSums(bsres)/k
selgenes <- (1:length(inclfreq))[inclfreq>=0.01]
selrows = rownames(x)[selgenes]

ycol = matrix(cols[ybin+1],length(ybin),1)
fname = sprintf("%s_heatmap_selected_peptides.pdf",todaysdate)
pdf(fname,width=15,height=15,pointsize=24)
heatmap.3(t(xbin[,selgenes]),col=matcolorposneg,ColSideColors=ycol,margins=c(10,10))
dev.off()

pca = prcomp(xbin[,selgenes],scale.=FALSE,center=FALSE)
fname = sprintf("%s_pca_selected_features.pdf",todaysdate)
pdf(fname,width=12,height=12,pointsize=24)
par(mfrow=c(1,1),mar=c(3,3,3,3))
pairs(pca$x[,1:5],col=c('red','steelblue1','olivedrab3')[as.factor(ybin)],pch=19)
pairs(pca$x[,1:5],col=c('red','steelblue1','olivedrab3')[as.factor(batch)][ybinbool],pch=19)
dev.off()




write.csv(file="selected_peptides_without_label_complete.csv",peptide.data[as.numeric(selrows),])

