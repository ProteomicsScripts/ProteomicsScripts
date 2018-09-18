## This is an R script for the conversion of mzTab to a better readable tsv format

# clear entire workspace
rm(list = ls())

# input data as mzTab
input.file <- 'misc/example.mzTab'

# output data as tsv
output.file <- 'misc/example.tsv'

# count the occurences of character c in string s
countOccurrences <- function(c,s) {
  s2 <- gsub(c,"",s)
  return (nchar(s) - nchar(s2))
}

# check that all protein accessions are of the format *|*|*
checkAccessionFormat <- function(accessions) {
  n <- length(accessions)
  count <- countOccurrences("[|]",accessions)
  m <- length(which(count==2))
  return (n==m)
}

# collapse rows
# (In mzTab, PSMs with multiple protein accessions are reported in multiple rows. This function collapses them to a single row.)
collapseRows <- function(psm.data) {
  
  # generate index vector idx
  tmp.psm.id <- 0
  idx <- c()
  accessions.tmp <- c()
  accessions.strings <- c()
  for (i in 1:length(psm.data$PSM_ID)) {
    
    if (psm.data$PSM_ID[i] == tmp.psm.id) {
      if (length(accessions.tmp) > 0) {
        accessions.strings <- c(accessions.strings, toString(accessions.tmp, sep=','))
        accessions.tmp <- c()
      }
      
      idx <- c(idx,i)
      tmp.psm.id <- tmp.psm.id + 1
    }
    
    accessions.tmp <- c(accessions.tmp, psm.data$accession[i])
  }
  accessions.strings <- c(accessions.strings, toString(accessions.tmp, sep=','))
  
  psm.data <- psm.data[idx,]
  psm.data$accession <- accessions.strings
  
  return (psm.data)
}

# read the PSM section of an mzTab file
readMzTabPSM <- function(file) {
  
  # read entire mzTab
  no.col <- max(count.fields(file, sep = "\t", quote=""))
  data <- read.table(file, sep="\t", fill=TRUE, quote="", col.names=1:no.col)
  
  # extract PSM data
  psm.data <- data[which(data[,1]=="PSM"),]
  colnames(psm.data) <- unlist(data[which(data[,1]=="PSH")[1],])
  psm.data$PSH <- NULL
  
  # simplify accession (in case it is of the format *|*|* )
  psm.data$accession <- as.character(psm.data$accession)
  if (checkAccessionFormat(psm.data$accession)) {
    list <- strsplit(psm.data$accession, "[|]")
    psm.data$accession <- unlist(lapply(list, '[[', 2))
    psm.data$gene <- unlist(lapply(list, '[[', 3))
  }
  
  psm.data <- collapseRows(psm.data)
  
  return (psm.data)
}

# check whether the column search_engine_score[n] exists
checkSearchEngineScoreExists <- function(table, n) {
  column <- paste("search_engine_score[", as.character(n), "]", sep="")
  return (column %in% colnames(table))
}

# check whether the column opt_global_SpecEValue_score exists
checkEValueExists <- function(table) {
  column <- "opt_global_SpecEValue_score"
  return (column %in% colnames(table))
}

# returns index of the best quantification with this sequence
# requires a 'opt_global_modified_sequence' and 'search_engine_score[1]' column in data frame 'data'
indexBest <- function(sequence, data) {
  idx <- which(data$opt_global_modified_sequence==sequence)
  min <- min(as.numeric(as.character(data$`search_engine_score[1]`[idx])))
  idx.m <- which(data[idx,]$`search_engine_score[1]`==min)
  return(idx[idx.m])
}

# makes the sequences unique by picking the quants with maximum intensity
makeSequencesUnique <- function(peptide.data) {
  unique.sequences <- unique(peptide.data$"opt_global_modified_sequence")
  idx <- unlist(lapply(unique.sequences, FUN=indexBest, data = peptide.data))
  return(peptide.data[idx,])
}





























# read mzTab data
peptide.data <- readMzTabPSM(input.file)

# make sequence unique
peptide.data <- makeSequencesUnique(peptide.data)

# write unqiue data as tsv
write.table(peptide.data, output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
