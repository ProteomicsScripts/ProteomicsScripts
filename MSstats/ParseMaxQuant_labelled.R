# install BioConductor
# source("http://bioconductor.org/biocLite.R")

# install MSstats
# biocLite("MSstats")

library(MSstats)

# Try to follow MaxQtoMSstatsFormat example in http://127.0.0.1:13881/library/MSstats/doc/MSstats.html

# Read in MaxQuant files
proteinGroups <- read.table("txt_TMT/proteinGroups.txt", sep="\t", header=TRUE)
infile <- read.table("txt_TMT/evidence.txt", sep="\t", header=TRUE)

# Read in annotation including condition and biological replicates per run.
# Users should make this annotation file. It is not the output from MaxQuant.
annot <- read.csv("txt_TMT/annotation.csv", header=TRUE)

raw <- MaxQtoMSstatsFormat(evidence=infile, annotation=annot, proteinGroups=proteinGroups)

# resulting error:
# Error in MaxQtoMSstatsFormat(evidence = infile, annotation = annot, proteinGroups = proteinGroups) : 
# ** Run is not provided in Annotation. Please check the annotation file.

# The example in 4.3.1 here http://msstats.org/wp-content/uploads/2017/01/MSstats_v3.7.3_manual.pdf
# does not list a run column.
