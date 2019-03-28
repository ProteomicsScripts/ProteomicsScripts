## This R script generates a set of basic plots and tables which summarise the PSM section of mzTab files.
## It is divided into three parts:
## (1) definition of global options and parameters
## (2) definition of functions
## (3) main part i.e. generation of plots and tables
##
## To install dependencies, please run in R:
## install.packages("xtable")       # for mod summary table

library(xtable)

# clear entire workspace
rm(list = ls())

# options
options(digits=10)

#input.file <- 'analysis.mzTab'
input.file <- 'example.mzTab'

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

# function describing how to collapse rows
# In the case of string columns (e.g. accessions), the row entries are first made unique and then written to a comma-separated string.
# In all other cases, the entry of the first row is returned.
collapseRows <- function(x) {
  if (is.character(x)) {
    x <- paste(unique(x[!is.na(x)]), collapse=",")
    if (x=="") {
      return(NA)
    }
    else {
      return (x)
    }
  }
  else {
    return (x[1])
  }
}

# count the occurences of character c in string s
countOccurrences <- function(c,s) {
  s2 <- gsub(c,"",s)
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

# read the PSM section of an mzTab file
readMzTabPSM <- function(file) {
  
  # find start of the PSM section
  first.row <- startSection(file, "PSH")
  
  # read entire mzTab
  data <- read.table(file, sep="\t", skip=first.row-1, fill=TRUE, header=TRUE, quote="", na.strings=c("null","NA"), stringsAsFactors=FALSE, check.names=FALSE)
  
  # extract PSM data
  psm.data <- data[which(data[,1]=="PSM"),]
  psm.data$PSH <- NULL
  
  # In case the accession column is of the format *|*|*, we split this column into an accession and a gene column.
  if (all(sapply(psm.data$accession, checkAccessionFormat))) {
    psm.data$gene <- sapply(psm.data$accession, getGene)
    psm.data$accession <- sapply(psm.data$accession, getAccession)
  }
  
  # In the mzTab format, PSMs with multiple protein accessions are written to multiple rows.
  # Here we collapse these rows and separate multiple accessions/genes by comma.
  psm.data <- aggregate(psm.data, by=list(temp = psm.data$PSM_ID), FUN=collapseRows)
  psm.data$temp <- NULL
  
  return (psm.data)
}

# read the MTD section of an mzTab file
readMzTabMTD <- function(file) {
  
  # find start of the MTD section
  first.row <- startSection(file, "MTD")
  
  # read entire mzTab
  data <- read.table(file, sep="\t", skip=first.row-1, fill=TRUE, header=TRUE, quote="", na.strings=c("null","NA"), stringsAsFactors=FALSE, check.names=FALSE)
  
  # extract MTD data
  mtd.data <- data[which(data[,1]=="MTD"),]
  mtd.data$MTD <- NULL
  
  return (mtd.data)
}

# returns the search engine of the search engine score n
# and NA if it does not exists
getSearchEngine <- function(meta.data, n) {
  
  row <- paste("psm_search_engine_score[", as.character(n), "]", sep="")
  
  # extract the third entry
  x <- meta.data[which(meta.data[,1]==row)[1],2]
  x <- gsub("\\[", "", x)
  x <- gsub("\\]", "", x)
  x <- unlist(strsplit(x, ","))[3]
  x <- gsub(" ", "", x)
  
  return(x)
}

# create a summary table of all modifications and their specificities
# required input is a dataframe with a "sequence" and "modifications" column in mzTab standard
createModsSummary <- function(data)
{
  # extract relevant data
  data <- data[,c("sequence","modifications")]
  data <- data[!is.na(data$modifications),]
  data <- data[data$modifications != c(""),]
  
  # check if any mods are reported
  if (dim(data)[1] == 0)
  {
    stats <- t(data.frame(c("no mods reported","","")))
    colnames(stats) <- c("modification","specificity","number")
    rownames(stats) <- c()
    return(stats)
  }
  
  # split comma-separted mods into multiple columns
  all.mods <- strsplit(data$modifications, split=",")
  l <- sapply(all.mods, length)
  idx <- rep(1:length(l), l)
  data <- data[idx,]
  data$modifications <- unlist(all.mods)
  
  # extract specificity
  getSiteIndex <- function(mod)
  {
    return(unlist(strsplit(mod, split="-"))[1])
  }
  data$idx <- sapply(data$modifications, getSiteIndex)
  getSiteAA <- function(idx, sequence)
  {
    if (idx == 0)
    {
      return("N-term")
    }
    else
    {
      return(substr(sequence, idx, idx))
    }
  }
  data$specificity <- mapply(getSiteAA, idx = data$idx, sequence = data$sequence)
  
  # extract mod accession
  getModAccession <- function(mod)
  {
    return(unlist(strsplit(mod, split=":"))[2])
  }
  data$accession <- sapply(data$modifications, getModAccession)
  
  # create summary statistics
  stats <- aggregate(data$accession, by=list(data$accession, data$specificity), FUN=length)
  colnames(stats) = c("mod","specificity","number")
  stats <- stats[order(stats$number, decreasing = TRUE),]
  
  # replace mod accession by mod name
  Accession2Mod <- rep("",3000)
  Accession2Mod[1] <- "Acetyl"
  Accession2Mod[4] <- "Carbamidomethyl"
  Accession2Mod[34] <- "Methyl"
  Accession2Mod[35] <- "Oxidation"
  Accession2Mod[36] <- "Dimethyl"
  Accession2Mod[39] <- "Methylthio"
  Accession2Mod[188] <- "Label:13C(6)"
  Accession2Mod[199] <- "Dimethyl:2H(4)"
  Accession2Mod[259] <- "Label:13C(6)15N(2)"
  Accession2Mod[267] <- "Label:13C(6)15N(4)"
  Accession2Mod[284] <- "Methyl:2H(2)"
  Accession2Mod[329] <- "Methyl:2H(3)13C(1)"
  Accession2Mod[330] <- "Dimethyl:2H(6)13C(2)"
  Accession2Mod[425] <- "Dioxidation"
  Accession2Mod[510] <- "Dimethyl:2H(4)13C(2)"
  Accession2Mod[737] <- "TMT6plex"
  stats$mod <- Accession2Mod[as.numeric(stats$mod)]
  
  return(stats)
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

# plot score distribution
plotScoreDistribution <- function(scores, pdf.file, breaks.histogram=80, score.name=NULL)
{
  if (is.factor(scores))
  {
    scores <- as.numeric(as.character(scores))
  }
  
  # construct x-axis label
  if (is.null(score.name))
  {
    x.label <- expression('log'[10]*' score')
    scores <- log10(scores)
  }
  else if ((score.name=="OMSSA") || (score.name=="MS-GF+"))
  {
    x.label <- bquote('-log'[10]*' score   ' ~ '(' * .(score.name) * ')')
    scores <- -log10(scores)
  }
  else
  {
    x.label <- bquote('log'[10]*' score   ' ~ '(' * .(score.name) * ')')
    scores <- log10(scores)
  }
  
  pdf(file=pdf.file, height=4)
  hist(scores, xlab=x.label, ylab="frequency", freq=TRUE, main="", col="grey", breaks=breaks.histogram)
  dev.off()
}





psm.data <- readMzTabPSM(input.file)

mtd.data <- readMzTabMTD(input.file)

# create mod summary statistics
stats <- createModsSummary(psm.data)

n.total <- dim(psm.data)[1]
n.unique <- length(which(psm.data$unique==1))
n.nonredundant <- length(unique(psm.data$sequence))
n.unique.nonredundant <- length(unique(psm.data[which(psm.data$unique==1),]$sequence))

n.target <- length(which(psm.data$opt_global_target_decoy=="target"))
n.decoy <- length(which(psm.data$opt_global_target_decoy=="decoy"))
n.target.decoy <- length(which(psm.data$opt_global_target_decoy=="target+decoy"))
n.neither <- length(which(psm.data$opt_global_target_decoy=="null"))
n.target.nonredundant <- length(unique(psm.data[which(psm.data$opt_global_target_decoy=="target"),]$sequence))
n.decoy.nonredundant <- length(unique(psm.data[which(psm.data$opt_global_target_decoy=="decoy"),]$sequence))
n.target.decoy.nonredundant <- length(unique(psm.data[which(psm.data$opt_global_target_decoy=="target+decoy"),]$sequence))
n.neither.nonredundant <- length(unique(psm.data[which(psm.data$opt_global_target_decoy=="null"),]$sequence))

if (checkSearchEngineScoreExists(psm.data, 1))
{
  scores <- psm.data$`search_engine_score[1]`
  score.name <- getSearchEngine(mtd.data, 1)

  if (length(scores) > 0)
  {
    plotScoreDistribution(scores, "plot__search_engine_score.pdf", score.name=score.name)
  }

  scores <- psm.data[which(psm.data$opt_global_target_decoy=="target"),]$`search_engine_score[1]`

  if (length(scores) > 0)
  {
    plotScoreDistribution(scores, "plot__search_engine_score__target.pdf", score.name=score.name)
  }

  scores <- psm.data[which(psm.data$opt_global_target_decoy=="decoy"),]$`search_engine_score[1]`

  if (length(scores) > 0)
  {
    plotScoreDistribution(scores, "plot__search_engine_score__decoy.pdf", score.name=score.name)
  }
}

if (checkEValueExists(psm.data))
{
  scores <- psm.data$`opt_global_SpecEValue_score`

  if (length(scores) > 0)
  {
    plotScoreDistribution(scores, "plot__e_value_score.pdf")
  }

  scores <- psm.data[which(psm.data$opt_global_target_decoy=="target"),]$`opt_global_SpecEValue_score`

  if (length(scores) > 0)
  {
    plotScoreDistribution(scores, "plot__e_value_score__target.pdf")
  }

  scores <- psm.data[which(psm.data$opt_global_target_decoy=="decoy"),]$`opt_global_SpecEValue_score`

  if (length(scores) > 0)
  {
    plotScoreDistribution(scores, "plot__e_value_score__decoy.pdf")
  }
}
