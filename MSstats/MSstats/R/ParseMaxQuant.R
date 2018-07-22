library("reshape2")
source("MaxQtoMSstatsFormat.R")

# Read in MaxQuant files
proteinGroups <- read.table("./txt_TMT/proteinGroups.txt", sep="\t", header=TRUE)

infile <- read.table("./txt_TMT/evidence.txt", sep="\t", header=TRUE)

# Read in annotation including condition and biological replicates per run.
# Users should make this annotation file. It is not the output from MaxQuant.
annot <- read.csv("./txt_TMT/dummy_annotation.csv", header=TRUE)

print(MaxQtoMSstatsFormat(evidence=infile, 
                          annotation=annot, 
                          proteinGroups=proteinGroups,
                          useUniquePeptide=FALSE, 
                          summaryforMultipleRows=max, 
                          fewMeasurements="keep", 
                          removeMpeptides=FALSE,
                          removeOxidationMpeptides=FALSE,
                          removeProtein_with1Peptide=FALSE))
