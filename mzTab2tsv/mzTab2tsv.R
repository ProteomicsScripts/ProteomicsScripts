## This is an R script for the conversion of mzTab to a better readable tsv format
## To install dependencies, run in R:
## install.packages("corrplot")

library("corrplot")

# clear entire workspace
rm(list = ls())

# input data as mzTab
input.file <- 'misc/example.mzTab'
#input.file <- '/Users/lars/Dropbox/UrinaryBiomarkers/UK181-209_ALL.mzTab'

# output data as tsv
output.file <- 'misc/example.tsv'
#output.file <- '/Users/lars/Dropbox/UrinaryBiomarkers/UK181-209_ALL.tsv'

# options and parameters
options(digits=10)
FcCutoff <- 8    # fold change cutoff, i.e. infinite fc values are mapped to +/-FcCutoff

# use the stripped or modified sequence in makeSequencesUnique()
stripped = FALSE                        

# Biognosys iRT spike-in peptides
peptides.of.interest <- c("LGGNEQVTR", "GAGSSEPVTGLDAK", "VEATFGVDESNAK", "YILAGVENSK", "TPVISGGPYEYR", "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR", "LFLQFGAQGSPFLK")

# proteins of interest
proteins.of.interest <- c("P0A853", "P37647")

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

# check if the vector/column is empty
# i.e. all entries are NA or "" etc. or the the vector is of length 0
isEmpty <- function(column)
{
  column <- column[!is.na(column)]
  column <- column[column != c("")]
  return(length(column) == 0)
}

# count the occurences of character c in string s
countOccurrences <- function(char,s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2))
}

# check that the protein accession is of the format *|*|*
# Note that NA returns TRUE.
checkAccessionFormat <- function(accession) {
  if (isEmpty(accession)) {
    return(FALSE)
  }
  
  if (all(is.na(accession))) {
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

# splits fasta protein accession into UniProt accession and gene name
splitAccession <- function(peptide.data) {
  # simplify accession (in case it is of the format *|*|* )
  peptide.data$accession <- as.character(peptide.data$accession)
  if (checkAccessionFormat(peptide.data$accession)) {
    list <- strsplit(peptide.data$accession,"[|]")
    peptide.data$accession <- unlist(lapply(list, '[[', 2))
    peptide.data$gene <- unlist(lapply(list, '[[', 3))
  }
  
  return (peptide.data)
}

# check if a specific "peptide_abundance_study_variable[n]" column exists
abundanceExists <- function(data, n)
{
  column <- paste("peptide_abundance_study_variable[",as.character(n),"]",sep="")
  return (column %in% colnames(data))
}

# returns the number of quantification channels i.e. the number of "peptide_abundance_study_variable[*]" columns
numberOfAbundances <- function(data)
{
  columns <- colnames(data)
  return(length(which(grepl("peptide_abundance_study_variable", columns))))
}

# determine fold changes and map to finite numbers
# fc = log2(abundances1/abundances2)
calculateFoldChange <- function(abundances1, abundances2) {
  offset <- 1e-10       # avoids devisions by zero
  max.fc <- FcCutoff    # map knock-out fold changes to finite values 
  
  # calculate fold changes
  abundances1 <- abundances1 + offset
  abundances2 <- abundances2 + offset
  fc <- log2(abundances1/abundances2)
  
  # map fc to finite values
  fc[which(fc < -FcCutoff)] <- -FcCutoff
  fc[which(fc > +FcCutoff)] <- +FcCutoff
  
  return (fc)
}

# plot fold change vs log intensity
plotFcLogIntensity <- function(fc.vector, intensity.vector, fc.label, pdf.file) {
  pdf(file=pdf.file)
  x <- fc.vector
  y <- log10(intensity.vector)
  df <- data.frame(x,y)
  x <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  #cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  cols <-  colorRampPalette(c("#2166AC", "#3F8EC0", "#80B9D8", "#BCDAEA", "#E6EFF3", "#F9EAE1", "#FAC8AF", "#ED9576", "#D25749", "#B2182B"))(256)
  df$col <- cols[df$dens]
  plot(y~x, data=df[order(df$dens),], pch=20, col=col, xlab=fc.label, ylab=expression('log'[10]*' intensity'))
  abline(v=0, col = "gray", lty=1)
  abline(v=median(fc.vector, na.rm=TRUE), col = "gray", lty=2)
  dev.off()
}

# plot distribution
plotDistribution <- function(vector, label, pdf.file) {
  pdf(file=pdf.file, height=4)
  density <- density(vector, na.rm=TRUE, bw="nrd0")
  plot(density$x, density$y, xlab=label, type="n", ylab="density", main="", yaxt='n')
  lines(density$x, density$y, col="gray", lwd=2)
  abline(v=0, col = "gray", lty=1)
  abline(v=median(vector, na.rm=TRUE), col = "gray", lty=2)
  dev.off()
}

# Kendrick nominal fractional mass plot
plotKendrick <- function(mass, pdf.file) {
  nominal = floor(mass)
  fractional = mass - nominal
  pdf(file=pdf.file)
  x <- nominal
  y <- fractional
  df <- data.frame(x,y)
  x <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  #cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  cols <-  colorRampPalette(c("#2166AC", "#3F8EC0", "#80B9D8", "#BCDAEA", "#E6EFF3", "#F9EAE1", "#FAC8AF", "#ED9576", "#D25749", "#B2182B"))(256)
  df$col <- cols[df$dens]
  plot(y~x, data=df[order(df$dens),], pch=20, col=col, xlab='nominal mass [Da]', ylab='fractional mass [Da]')
  dev.off()
}

# returns a dataframe containing only the peptide quantification columns
# required input is a dataframe with columns "peptide_abundance_study_variable[*]"
getPeptideQuants <- function(data)
{
  idx <- grepl("peptide_abundance_study_variable", colnames(data))
  return(data[,idx])
}

plotCorrelations <- function(data, pdf.file) {
    # extract study variables
    study_variables.n <- numberOfAbundances(data)
    study_variables.data = getPeptideQuants(data)
    
    # (optional) z-score normalisation
    #study_variables.data <- scale(study_variables.data, center = TRUE, scale = TRUE)

    corr = cor(study_variables.data[complete.cases(study_variables.data),], method="pearson")    # possible methods: "pearson", "spearman", "kendall"
    
    # rename columns and rows
    colnames(corr) <- 1: study_variables.n
    rownames(corr) <- 1: study_variables.n
    cols <- colorRampPalette(c("#2166AC", "#3F8EC0", "#80B9D8", "#BCDAEA", "#E6EFF3", "#F9EAE1", "#FAC8AF", "#ED9576", "#D25749", "#B2182B"))(256)
    
    pdf(file=pdf.file)
    if (study_variables.n < 12)
    {
      # use combined "number/circle" plotting method for SILAC, TMT and small LFQ analyses
      corrplot(corr, cl.lim=c(min(corr),max(corr)), col = cols, is.corr=FALSE, method = "number")
      #corrplot.mixed(corr, cl.lim=c(min(corr), max(corr)), cols=cols, is.corr=FALSE, lower="number", upper="circle")
    }
    else
    {
      # use "circle" plotting method for LFQ analyses
      corrplot(corr, cl.lim=c(min(corr),max(corr)), col = cols, is.corr=FALSE, method = "circle")
    }
    dev.off()
    
    return(corr)
}

plotBoxplot <- function(data, pdf.file) {
  # extract study variables
  study_variables.data = getPeptideQuants(data)
  colnames(study_variables.data) <- as.character(1:(dim(study_variables.data)[2]))
  
  # (optional) z-score normalisation
  #study_variables.data <- scale(study_variables.data, center = TRUE, scale = TRUE)

  pdf(file=pdf.file, height = 6, width = 10)
  boxplot(study_variables.data, log="y", ylab="expression", xlab="samples", las=2)
  dev.off()
}

findPeptidesOfInterest <- function(data)
{
  retain.columns=c("sequence", "accession", "charge", "retention_time", "mass_to_charge")
  new.column.names=c("Sequence", "Accession", "Charge", "Retention Time", "m/z" )
  
  # check if sequence column is non-empty
  if (isEmpty(data$sequence))
  {
    df <- t(data.frame(c("no sequences reported", rep("", 4))))
    colnames(df) <- c("modified sequence", "accession", "charge", "retention time", "m/z" )
    rownames(df) <- c()
    return(df)
  }
  
  pattern = paste(peptides.of.interest, collapse="|")
  df <- as.data.frame(data[grepl(pattern, data$sequence),])
  
  # check if results are empty
  if (dim(df)[1] == 0)
  {
    df <- t(data.frame(c("no matching sequences found", rep("", 4))))
    colnames(df) <- c("modified sequence", "accession", "charge", "retention time", "m/z" )
    rownames(df) <- c()
    return(df)
  }
  
  # sort in the same order as peptides.of.interest vector
  df <- df[order(match(df$sequence, peptides.of.interest)),]
  
  df <- df[,c("opt_global_modified_sequence", "accession", "charge", "retention_time", "mass_to_charge")]
  colnames(df) <- c("modified sequence", "accession", "charge", "retention time", "m/z" )
  
  return(df)
}

findProteinsOfInterest <- function(data) {
    # check if protein accession column is non-empty
  if (isEmpty(data$accession))
  {
    df <- t(data.frame(c("", "no accessions reported", rep("", 3))))
    colnames(df) <- c("modified sequence", "accession", "charge", "retention time", "m/z" )
    rownames(df) <- c()
    return(df)
  }

  pattern = paste(proteins.of.interest, collapse="|")
  df <- as.data.frame(data[grepl(pattern, data$accession),])
  
  # check if results are empty
  if (dim(df)[1] == 0)
  {
    df <- t(data.frame(c("", "no matching accessions found", rep("", 3))))
    colnames(df) <- c("modified sequence", "accession", "charge", "retention time", "m/z" )
    rownames(df) <- c()
    return(df)
  }
  
  # sort sequences in alphabetic order
  df <- df[order(df$sequence),]
  
  # sort in the same order as proteins.of.interest vector
  df <- df[order(match(df$accession, proteins.of.interest)),]
  df <- df[,c("opt_global_modified_sequence", "accession", "charge", "retention_time", "mass_to_charge")]
  colnames(df) <- c("modified sequence", "accession", "charge", "retention time", "m/z" )
  
  return(df)
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

# returns index of the best quantification with this sequence
# requires a 'sequence' and 'intensity' column in data frame 'data'
indexMaxIntensity <- function(sequence, data) {
  if (stripped)
  {
    idx <- which(data$sequence==sequence)
  }
  else
  {
    idx <- which(data$opt_global_modified_sequence==sequence)
  }
  max <- max(data$intensity[idx])
  idx.m <- which(data[idx,]$intensity==max)
  return(idx[idx.m])
}

# makes the sequences unique by picking the quants with maximum intensity
makeSequencesUnique <- function(peptide.data) {
  if (stripped)
  {
    unique.sequences <- unique(peptide.data$sequence)
  }
  else
  {
    unique.sequences <- unique(peptide.data$opt_global_modified_sequence)
  }
  idx <- unlist(lapply(unique.sequences, FUN=indexMaxIntensity, data = peptide.data))
  return(peptide.data[idx,])
}





























# read mzTab data
peptide.data <- readMzTabPEP(input.file)

# remove decoy and contaminant hits and split accession
if (!isEmpty(peptide.data$accession))
{
  peptide.data <- peptide.data[which(substr(peptide.data$accession,1,4)!="dec_"),]
  peptide.data <- peptide.data[which(substr(peptide.data$accession,1,4)!="CON_"),]

  # Note that decoys and contaminants might not be of the form *|*|* and accessions might not have been split in readMzTabPEP().
  # Hence we split the accessions here again after removing decoys and accessions.
  peptide.data <- splitAccession(peptide.data)
}

# calculate peptide intensity as mean of all study variables
peptide.data$intensity <- (peptide.data$"peptide_abundance_study_variable[1]" + peptide.data$"peptide_abundance_study_variable[2]" + peptide.data$"peptide_abundance_study_variable[3]" + peptide.data$"peptide_abundance_study_variable[4]")/4

# get identified peptides
# make them unique i.e. if multiple quants for the same sequence exist, keep only the one with the highest intensity
peptide.data.identified <- peptide.data[which(!is.na(peptide.data$sequence)),]
peptide.data.unique <- makeSequencesUnique(peptide.data.identified)

# write unqiue data as tsv
write.table(peptide.data.unique, output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# get peptides of interest
peptides.of.interest <- findPeptidesOfInterest(peptide.data)

# total number of quantified and identified peptides
n.peptides <- dim(peptide.data)[1]
n.peptides.identified <- dim(peptide.data.identified)[1]
n.peptides.identified.modified.unique <- length(unique(peptide.data.identified$opt_global_modified_sequence))
n.peptides.identified.stripped.unique <- length(unique(peptide.data.identified$sequence))

# create mod summary statistics
mod.stats <- createModsSummary(peptide.data)

# plot correlation matrix of peptide abundances
corr <- plotCorrelations(peptide.data, "misc/plot_Correlations.pdf")
corr.min <- min(corr)
corr.median <- median(corr)
corr.max <- max(corr)

# plot boxplot of peptide abundances
plotBoxplot(peptide.data, "misc/plot_Boxplot.pdf")
