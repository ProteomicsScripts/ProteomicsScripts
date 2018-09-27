# install BioConductor
# source("http://bioconductor.org/biocLite.R")
# install.packages("tidyr")

# install MSstats
# biocLite("MSstats")

library(MSstats)
library(tidyr)

# clear entire workspace
rm(list = ls())

# Try to follow MaxQtoMSstatsFormat example in https://bioconductor.org/packages/release/bioc/vignettes/MSstats/inst/doc/MSstats.html

# Read in MaxQuant files
proteinGroups <- read.table("txt/proteinGroups.txt", sep="\t", header=TRUE)
infile <- read.table("txt/evidence.txt", sep="\t", header=TRUE)

# Read in annotation including condition and biological replicates per run.
# Users should make this annotation file. It is not the output from MaxQuant.
annot <- read.csv("annotation.csv", header=TRUE)

# Read MaxQuant results
raw <- MaxQtoMSstatsFormat(evidence=infile, annotation=annot, proteinGroups=proteinGroups)

# Keep only relevant columns
raw <- raw[,c("PeptideSequence", "PrecursorCharge", "Run", "Intensity")]

# long (tidy) to wide format using spread() from tidyr package
raw.wide <- spread(data = raw, key = Run, value = Intensity)

# Write imported data
write.table(raw.wide, "output.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
