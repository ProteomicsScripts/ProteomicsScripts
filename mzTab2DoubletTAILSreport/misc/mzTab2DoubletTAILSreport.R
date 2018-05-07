## R script for triple TAILS experiments, input mzTab

# install packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("UniProt.ws")
#install.packages('rasterVis')

rm(list = ls())

library(UniProt.ws)
library(rasterVis)

input.file <- "data.mzTab"
output.file <- "data.tsv"

# options
options(digits=10)

# fc cutoff, i.e. infinite fc values are mapped to +/-FcCutoff
FcCutoff <- 8

# amino acid vicinity for the frequency plots
aa.vicinity <- 6

plot.HL <- "FcLogIntensity_HL.pdf"

frequency.heatmap.HL.mInf <- "frequency_heatmap_HL_mInf.pdf"
frequency.heatmap.HL.zero <- "frequency_heatmap_HL_zero.pdf"
frequency.heatmap.HL.pInf <- "frequency_heatmap_HL_pInf.pdf"

all.amino.acids <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")

species <- 9606    # homo sapiens
columns <- c("SEQUENCE","GO", "SUBCELLULAR-LOCATIONS", "PROTEIN-NAMES", "GENES", "KEGG")

# count the occurences of character c in string s
countOccurrences <- function(char,s) {
	s2 <- gsub(char,"",s)
	return (nchar(s) - nchar(s2))
}

# check that all protein accessions are of the format *|*|*
checkAccessionFormat <- function(accessions) {
	n <- length(accessions)
	count <- countOccurrences("[|]",accessions)
	m <- length(which(count==2))
	return (n==m)
}

# read the PEP section of an mzTab file
readMzTabPEP <- function(file) {
  
  # read entire mzTab
  no.col <- max(count.fields(file, sep = "\t", quote=""))
  data <- read.table(file,sep="\t",fill=TRUE, quote="", col.names=1:no.col, stringsAsFactors=FALSE)
  
  # extract PEP data
  peptide.data <- data[which(data[,1]=="PEP"),]
  colnames(peptide.data) <- unlist(data[which(data[,1]=="PEH")[1],])
  peptide.data$PEH <- NULL
  
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

# make the input unique and combine them in a space separated string
UniqueSpaceSeparated <- function(x) {
  u <- unique(x)
  r <- paste(u, collapse=" ")
  return(r)
}

# annotate protein accessions with columns from UniProt
annotateAccession <- function(data) {
  keytype <- "UNIPROTKB"
  keys <- data$accession
  a <- UniProt.ws::select(up, keys, columns, keytype)
  b <- aggregate(a, by=list(a[,keytype]), FUN=UniqueSpaceSeparated)
  b[,1] <- NULL
  idx <- match(keys,b[,keytype])
  annotations <- b[idx,]
  annotations$UNIPROTKB <- NULL
  data <- cbind(data, annotations)
  return (data)
}

# determine intensities and fold changes and map to finite numbers
# H=0 L=finite              =>  mapped to fc=-FcCutoff
# H=finite L=0              =>  mapped to fc=+FcCutoff
# H=0 L=0 (i.e. M=finite)   =>  mapped to fc=0
calculateIntensityFC <- function(peptide.data) {
  offset <- 1    # avoids devisions by zero
  max.fc <- FcCutoff    # map knock-out fold changes to finite values 
  
  # convert to numeric
  peptide.data$"peptide_abundance_study_variable[1]" <- as.numeric(peptide.data$"peptide_abundance_study_variable[1]")
  peptide.data$"peptide_abundance_study_variable[2]" <- as.numeric(peptide.data$"peptide_abundance_study_variable[2]")

  # calculate intensity of doublets (simple mean of the two abundances)
  peptide.data$intensity = (peptide.data$"peptide_abundance_study_variable[1]" + peptide.data$"peptide_abundance_study_variable[2]")/2
  
  # calculate fold changes
  peptide.data$"peptide_abundance_study_variable[1]" <- peptide.data$"peptide_abundance_study_variable[1]" + offset
  peptide.data$"peptide_abundance_study_variable[2]" <- peptide.data$"peptide_abundance_study_variable[2]" + offset
  peptide.data$fc.H.L <- log2(peptide.data$"peptide_abundance_study_variable[2]"/peptide.data$"peptide_abundance_study_variable[1]")

  # map fc to finite values
  peptide.data[which(peptide.data$fc.H.L < -FcCutoff),]$fc.H.L <- -FcCutoff
  peptide.data[which(peptide.data$fc.H.L > FcCutoff),]$fc.H.L <- FcCutoff

  return (peptide.data)
}

# returns index of the best quantification with this sequence
indexMaxIntensity <- function(sequence) {
  idx <- which(peptide.data$sequence==sequence)
  max <- max(peptide.data$intensity[idx])
  idx.m <- which(peptide.data[idx,]$intensity==max)
  return(idx[idx.m])
}

# makes the sequences unique by picking the quants with maximum intensity
makeSequencesUnique <- function(peptide.data) {
  unique.sequences <- unique(peptide.data$sequence)
  idx <- unlist(lapply(unique.sequences, FUN=indexMaxIntensity))
  return(peptide.data[idx,])
}

# return position where peptide sequence and protein sequence match
matchingPosition <- function(peptide.sequence,protein.sequence) {
  pos <- regexpr(peptide.sequence,protein.sequence)
  return(pos)
}

# remove protein N-termini (i.e. peptides matched to position 1 to 6 in the protein)
# and add (non)prime regions
addSequenceVicinity <- function(peptide.data) {
  peptide.data$matching.position <- mapply(matchingPosition, peptide.data$sequence, peptide.data$SEQUENCE)
  idx <- which(peptide.data$matching.position > 6)
  peptide.data <- peptide.data[idx,]
  peptide.data$nonprime.sequence <- mapply(substr, peptide.data$SEQUENCE, peptide.data$matching.position-6, peptide.data$matching.position-1)
  peptide.data$prime.sequence <- mapply(substr, peptide.data$SEQUENCE, peptide.data$matching.position, peptide.data$matching.position+5)
  return(peptide.data)
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
  abline(v=0, col = "gray")
  dev.off()
}

# return number of positions where peptide sequence and protein sequence match
numberOfPositions <- function(peptide.sequence,protein.sequence) {
  n <- length(gregexpr(peptide.sequence,protein.sequence))
  return(n)
}

# generate frequency scaling factors for each amino acid
generateGlobalAminoAcidFrequency <- function(proteins) {
  x <- paste(unique(proteins), collapse="")
  x <- substring(x, seq(1,nchar(x), 1), seq(1,nchar(x),1))
  t <- as.data.frame(table(x))

  # remove non-standard aminoacids Selenocysteine (U) and Pyrrolysine (O)
  t <- t[which(t[,1]!="U"),]
  t <- t[which(t[,1]!="O"),]

  freq <- t$Freq/sum(t$Freq)
  return(freq)
}

# generate frequency matrix
generateFrequencyMatrix <- function(protein.sequence, position, scaling, pdf.file) {
  
  table <- as.data.frame(all.amino.acids)
  for (i in -6:5)
  {
    aa <- mapply(substr, protein.sequence, position+i, position+i)
    df <- as.data.frame(table(aa))
    count <- rep(0,20)
    idx <- match(df[,1],all.amino.acids)
    count[idx] <- df[,2]
    count <- scaling * count
    table <- cbind(table,count)
  }
  
  rownames(table) <- table$all.amino.acids
  table <- table[,-c(1)]
  colnames(table) <- c("P6","P5","P4","P3","P2","P1","P1'","P2'","P3'","P4'","P5'","P6'")
  
  table <- table[c(20:1),]
  table_matrix <- data.matrix(table)

  pdf(file=pdf.file)
  #my.colours <- colorRampPalette(c("white", "grey", "red"))(n = 299)
  #heatmap(table_matrix, Rowv=NA, Colv=NA, col=my.colours, scale="none", margins=c(5,10), add.expr = abline(h=0, v=6.5))
  
  my.theme <- BuRdTheme()
  my.min <- min(table_matrix)
  my.max <- max(table_matrix)
  my.at <- seq(my.min, my.max, length.out=length(my.theme$regions$col)-1)
  my.ckey <- list(at=my.at, col=my.theme$regions$col)
  
  table_matrix <- t(table_matrix)
  p <- levelplot(table_matrix, par.settings=my.theme, at=my.at, colorkey=my.ckey, xlab="", ylab="",
    scales = list(tck = c(0,0)),    # removes ticks from axes
    panel = function(...){    # add vertical line between P1 and P1'
    panel.levelplot(...)
    panel.abline(v = 6.5)
  })
  print(p)
  
  dev.off()
}

# check species
checkSpecies <- function(species, proteins) {
  # Do all protein accessions have the substring <species> in their name?
  return(length(proteins) == length(proteins[grepl(species,proteins)]))
}









peptide.data <- readMzTabPEP(input.file)

# remove null accessions (TODO: check why there are nulls, FYVPGVAPINFHQND -> should be Q92544 in BM2321 and BM2322)
peptide.data <- peptide.data[which(peptide.data$accession!="null"),]

# remove decoy and contaminant hits
peptide.data <- peptide.data[which(substr(peptide.data$accession,1,4)!="dec_"),]
peptide.data <- peptide.data[which(substr(peptide.data$accession,1,4)!="CON_"),]

# total number of quantified peptides
n.total <- dim(peptide.data)[1]

# try to determine the species from the protein accessions
if (checkSpecies("MOUSE",peptide.data$accession))
{
  species <- 10090
}
if (checkSpecies("HUMAN",peptide.data$accession))
{
  species <- 9606
}

# load UniProt database for this species
up <- UniProt.ws(taxId=species)

peptide.data <- splitAccession(peptide.data)
peptide.data <- annotateAccession(peptide.data)

# remove unannotated peptides
peptide.data <- peptide.data[which(peptide.data$SEQUENCE!="NA"),]

# simplify table
#peptide.data <- peptide.data[,c("sequence", "best_search_engine_score[1]", "opt_global_XTandem_score", "accession", "peptide_abundance_study_variable[1]", "peptide_abundance_study_variable[2]", "peptide_abundance_study_variable[3]")]

# calculate intensities and fold changes
peptide.data <- calculateIntensityFC(peptide.data)

# pick the quants with highest intensity
peptide.data <- makeSequencesUnique(peptide.data)

# number of unique, quantified peptides
n.unique <- dim(peptide.data)[1]

# add N-term acetylation column
peptide.data$n.term.acetylation <- (grepl("0-UniMod:1,", peptide.data$modifications) | (peptide.data$modifications=="0-UniMod:1"))

# add (non)prime regions
peptide.data <- addSequenceVicinity(peptide.data)

# add column with amino acid at position P1 and P1'
peptide.data$P1 <- substr(peptide.data$nonprime.sequence, nchar(peptide.data$nonprime.sequence), nchar(peptide.data$nonprime.sequence))
peptide.data$P1prime <- substr(peptide.data$prime.sequence, 1, 1)

# number of unique, quantified peptides
n.after.vicinity.mapping <- dim(peptide.data)[1]

# plot fold change vs log intensity
plotFcLogIntensity(peptide.data$fc.H.L, peptide.data$intensity, "fc (H:L)", plot.HL)

# global amino acid frequency
frequency <- generateGlobalAminoAcidFrequency(peptide.data$SEQUENCE)
frequency.scaling <- 1/frequency
frequency.scaling <- frequency.scaling/sum(frequency.scaling)


idx.HL.mInf <- which(peptide.data$fc.H.L == -FcCutoff)
idx.HL.zero <- which((peptide.data$fc.H.L > -FcCutoff) & (peptide.data$fc.H.L < FcCutoff))
idx.HL.pInf <- which(peptide.data$fc.H.L == FcCutoff)

p <- peptide.data[idx.HL.mInf,]
generateFrequencyMatrix(p$SEQUENCE, p$matching.position, frequency.scaling, frequency.heatmap.HL.mInf)
p <- peptide.data[idx.HL.zero,]
generateFrequencyMatrix(p$SEQUENCE, p$matching.position, frequency.scaling, frequency.heatmap.HL.zero)
p <- peptide.data[idx.HL.pInf,]
generateFrequencyMatrix(p$SEQUENCE, p$matching.position, frequency.scaling, frequency.heatmap.HL.pInf)

write.table(peptide.data, output.file, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
