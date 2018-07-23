library("reshape2")  # this is a necessary dependency to execute MaxQtoMSstatsFormat.R
source("MSstats/R/MaxQtoMSstatsFormat.R")


# this closely follows the example in MSstats docu
printMaxQtoMSstats <- function(foldername) {
    # Read in MaxQuant files
    proteinGroups <- read.table(sprintf("%s/proteinGroups.txt", foldername), sep="\t", header=TRUE)

    infile <- read.table(sprintf("%s/evidence.txt", foldername), sep="\t", header=TRUE)

    # Read in annotation including condition and biological replicates per run.
    # Users should make this annotation file. It is not the output from MaxQuant.
    annot <- read.table(sprintf("%s/annotation.csv", foldername), sep=",",  header=TRUE)
    print(annot)

    print(MaxQtoMSstatsFormat(evidence=infile, 
                              annotation=annot, 
                              proteinGroups=proteinGroups,
                              useUniquePeptide=FALSE, 
                              summaryforMultipleRows=max, 
                              fewMeasurements="keep", 
                              removeMpeptides=FALSE,
                              removeOxidationMpeptides=FALSE,
                              removeProtein_with1Peptide=FALSE))

}

printMaxQtoMSstats("./txt_TMT")
