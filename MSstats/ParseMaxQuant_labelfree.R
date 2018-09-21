# install BioConductor
# source("http://bioconductor.org/biocLite.R")

# install MSstats
# biocLite("MSstats")

library(MSstats)

# Try to follow MaxQtoMSstatsFormat example in http://127.0.0.1:13881/library/MSstats/doc/MSstats.html

# Read in MaxQuant files
proteinGroups <- read.table("txt_labelfree/proteinGroups.txt", sep="\t", header=TRUE)
infile <- read.table("txt_labelfree/evidence.txt", sep="\t", header=TRUE)

# Read in annotation including condition and biological replicates per run.
# Users should make this annotation file. It is not the output from MaxQuant.
annot <- read.csv("txt_labelfree/annotation.csv", header=TRUE)

raw <- MaxQtoMSstatsFormat(evidence=infile, annotation=annot, proteinGroups=proteinGroups)
