# helper functions for preprocessing and plotting of proteomics data in mzTab format
# J. Griss et al. (2014). The mzTab Data Exchange Format: Communicating Mass-spectrometry-based Proteomics and Metabolomics Experimental Results to a Wider Audience. Molecular & Cellular Proteomics, 13(10), 2765–2775. http://doi.org/10.1074/mcp.O113.036681

# install.packages("corrplot")     # for correlation of peptide intensities
# install.packages("xtable")       # for peptides/proteins of interest tables
# install.packages("rasterVis")
# install.packages("reshape2")

#' Simple Moving Average without leading NA.
#'
#' @param file path to mzTab file
#' @param section.identifier identifier at the start of the section, either 'PEH', 'PRH' or 'PSH'
#'
#' @return first row of the section
#'
#' @export
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

#' Check if the vector/column is empty.
#'
#' i.e. all entries are NA or "" etc. or the the vector is of length 0
#'
#' @param column column to be checked
#'
#' @return Boolean. Is the entire column empty?
#'
#' @export
isEmpty <- function(column)
{
  column <- column[!is.na(column)]
  column <- column[column != c("")]
  return(length(column) == 0)
}

#' Count the occurences of character char in string s.
#'
#' @param char character
#' @param s string
#'
#' @return number of occurences
#'
#' @export
countOccurrences <- function(char, s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2))
}

#' Check that the protein accession is of the format ..|..|..
#'
#' Note that NA returns TRUE.
#'
#' @param accession protein accession
#'
#' @return Boolean. Is the accession of the format ..|..|.. ?
#'
#' @export
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

#' Extracts the second entry from a string of the form ..|..|..
#'
#' @param string string of the form x|y|z
#'
#' @return y
#'
#' @export
getAccession <- function(string) {
  if (is.na(string)) {
    return (NA)
  }
  return (unlist(strsplit(string, "[|]"))[2])
}

#' Extracts the third entry from a string of the form ..|..|..
#'
#' @param string string of the form x|y|z
#'
#' @return z
#'
#' @export
getGene <- function(string) {
  if (is.na(string)) {
    return (NA)
  }
  return (unlist(strsplit(string, "[|]"))[3])
}

#' Read the PEP section of an mzTab file.
#'
#' @param file path to mzTab file
#'
#' @return dataframe of PEP section
#'
#' @export
readMzTabPEP <- function(file) {

  # find start of the PEP section
  first.row <- startSection(file, "PEH")

  # read a single row and determine class types
  single.row <- read.table(file, sep="\t", skip=first.row-1, nrow=1, fill=TRUE, header=TRUE, quote="", na.strings=c("null","NA"), stringsAsFactors=FALSE, check.names=FALSE)
  class.types <- sapply(single.row, class)

  # read entire mzTab
  data <- read.table(file, sep="\t", skip=first.row-1, fill=TRUE, header=TRUE, quote="", na.strings=c("null","NA"), stringsAsFactors=FALSE, check.names=FALSE)

  # extract PEP data
  peptide.data <- data[which(data[,1]=="PEP"),]
  for (i in 1:length(class.types))
  {
    if (class.types[i]=="character")
    {
      peptide.data[,i] <- as.character(peptide.data[,i])
    }
    if (class.types[i]=="numeric")
    {
      peptide.data[,i] <- as.numeric(peptide.data[,i])
    }
    if (class.types[i]=="logical")
    {
      peptide.data[,i] <- as.logical(peptide.data[,i])
    }
  }
  # set the class types for 'unique' and 'charge' explicitly
  peptide.data$unique <- as.numeric(peptide.data$unique)
  peptide.data$charge <- as.numeric(peptide.data$charge)
  peptide.data$PEH <- NULL

  # In case the accession column is of the format *|*|*, we split this column into an accession and a gene column.
  if (all(sapply(peptide.data$accession, checkAccessionFormat))) {
    peptide.data$gene <- sapply(peptide.data$accession, getGene)
    peptide.data$accession <- sapply(peptide.data$accession, getAccession)
  }

  return (peptide.data)
}

#' Splits fasta protein accession into UniProt accession and gene name.
#'
#' @param peptide.data dataframe with <accession> column
#'
#' @return dataframe with UniProt accession and gene name
#'
#' @export
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

#' Check if a specific "peptide_abundance_study_variable[n]" column exists.
#'
#' @param data dataframe
#' @param n index of study variable
#'
#' @return Boolean. Does the column "peptide_abundance_study_variable[n]" exist?
#'
#' @export
studyVariableExists <- function(data, n)
{
  column <- paste("peptide_abundance_study_variable[",as.character(n),"]",sep="")
  return (column %in% colnames(data))
}

#' Returns the number of quantification channels i.e. the number of "peptide_abundance_study_variable[*]" columns.
#'
#' @param data dataframe
#'
#' @return number of quantification channels
#'
#' @export
numberOfStudyVariables <- function(data)
{
  columns <- colnames(data)
  return(length(which(grepl("peptide_abundance_study_variable", columns))))
}

#' Determine fold changes and map to finite numbers.
#'
#' fc = log2(abundances1/abundances2)
#' Should be NaN-save, i.e. NaN in any of the two abundance vectors results in NaN in the fold change vector.
#' NaN = not quantifiable, hence fc = NaN = not quantifiable
#'
#' @param abundances1 first abundance vector
#' @param abundances2 second abundance vector
#'
#' @return fold change
#'
#' @export
calculateFoldChange <- function(abundances1, abundances2)
{
  offset <- 1e-10       # avoids devisions by zero
  max.fc <- FcCutoff    # map knock-out fold changes to finite values

  # Alternatively, NaNs can be replaced by zeros.
  #abundances1[is.nan(abundances1)] <- 0
  #abundances2[is.nan(abundances2)] <- 0

  # calculate fold changes
  abundances1 <- abundances1 + offset
  abundances2 <- abundances2 + offset
  fc <- log2(abundances1/abundances2)

  # map fc to finite values
  fc[which(fc < -FcCutoff)] <- -FcCutoff
  fc[which(fc > +FcCutoff)] <- +FcCutoff

  return (fc)
}

#' Returns a unique colour for each string.
#'
#' @param string.vector vector of strings
#'
#' @return colours
#'
#' @export
uniqueColors <- function(string.vector)
{
  unique <- unique(string.vector)
  palette <- rainbow(length(unique))
  names(palette) <- unique

  return(palette[string.vector])
}

#' Plot the distribution of peptide elution times.
#'
#' Each peptide reports a minimum/maximum retention time in the retention_time_window column.
#'
#' @param data dataframe with retention_time_window column
#' @param pdf.file path to output pdf file
#'
#' @export
plotElutionTimeDistribution <- function(data, pdf.file)
{
  # split the input string at '|' and return the numeric element at a certain position.
  splitAtPipe <- function(string, position=1)
  {
    return(as.numeric(unlist(strsplit(string,split='\\|')))[position])
  }

  rt.window <- data$retention_time_window
  rt.min <- unlist(lapply(rt.window, splitAtPipe, position=1))
  rt.max <- unlist(lapply(rt.window, splitAtPipe, position=2))
  rt.width <- rt.max-rt.min

  plotDistribution(log10(rt.width), expression('log'[10]*' of peptide elution time [sec]'), "plot_DistributionElutionTime.pdf")
}

#' Plot fold change vs log intensity.
#'
#' @param fc.vector fold change vector (x-axis)
#' @param intensity.vector intensity vector (y-axis)
#' @param fc.label label for x-axis
#' @param pdf.file path to output pdf file
#'
#' @export
plotFcLogIntensity <- function(fc.vector, intensity.vector, fc.label, pdf.file)
{
  pdf(file=pdf.file)
  idx <- complete.cases(fc.vector) & complete.cases(intensity.vector)
  x <- fc.vector[idx]
  y <- log10(intensity.vector[idx])
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

#' Plot the difference between 2+ and 3+ fold changes of the same (modified) sequence against the logarithm of the average peptide intensity.
#'
#' @param data data frame
#' @param sample.1 number of the sample i.e. study variable
#' @param sample.2 number of the sample i.e. study variable
#' @param pdf.file path to output pdf file
#'
#' @export
plotDeltaFcLogIntensity <- function(data, sample.1, sample.2, pdf.file)
{
  column.1 <- paste("peptide_abundance_study_variable[", as.character(sample.1), "]", sep="")
  column.2 <- paste("peptide_abundance_study_variable[", as.character(sample.2), "]", sep="")

  # filter for relevant data (for run time improvement)
  data <- data[c("opt_global_modified_sequence", "charge", column.1, column.2)]

  # make (sequence, charge) unique
  # Only then can we match 2+ and 3+ quantifications unambiguously
  data <- makeModifiedSequenceChargeUnique(data)

  # calculate all (finite) fold changes for 2+ quantifications
  data.2 <- data[which(data$charge == 2),]
  data.2 <- data.2[which(data.2[[column.1]] > 0),]
  data.2 <- data.2[which(data.2[[column.2]] > 0),]
  data.2$fc_2 <- log2(data.2[[column.1]]/data.2[[column.2]])
  data.2$intensity_2 <- getAverageIntensity(data.2)
  data.2 <- data.2[c("opt_global_modified_sequence","fc_2","intensity_2")]

  # calculate all (finite) fold changes for 3+ quantifications
  data.3 <- data[which(data$charge == 3),]
  data.3 <- data.3[which(data.3[[column.1]] > 0),]
  data.3 <- data.3[which(data.3[[column.2]] > 0),]
  data.3$fc_3 <- log2(data.3[[column.1]]/data.3[[column.2]])
  data.3$intensity_3 <- getAverageIntensity(data.3)
  data.3 <- data.3[c("opt_global_modified_sequence","fc_3","intensity_3")]

  # find peptides quantified as both 2+ and 3+
  data <- merge(data.2, data.3, by="opt_global_modified_sequence")

  # determine simple statistics
  n <- dim(data)[1]
  sd.delta.fc <- sd(data$fc_2 - data$fc_3)

  if (n>3)
  {
    pdf(file=pdf.file)
    x <- data$fc_2 - data$fc_3
    y <- (data$intensity_2 + data$intensity_3)/2
    df <- data.frame(x,y)
    x <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
    df$dens <- col2rgb(x)[1,] + 1L
    cols <-  colorRampPalette(c("#2166AC", "#3F8EC0", "#80B9D8", "#BCDAEA", "#E6EFF3", "#F9EAE1", "#FAC8AF", "#ED9576", "#D25749", "#B2182B"))(256)
    df$col <- cols[df$dens]
    plot(y~x, data=df[order(df$dens),], pch=20, col=col, log="y", ylab="intensity", xlab="fc 2+ - fc 3+", main=paste("n = ",as.character(n),"    sd(fc 2+ - fc 3+) = ",as.character(round(sd.delta.fc, digits=4)),sep=""))
    dev.off()
  }
}

#' Plot peptide fold change vs log intensity for a single specific protein.
#'
#' (same peptide sequence -> same colour)
#'
#' @param data data frame
#' @param protein protein
#' @param sample.1 number of the sample i.e. study variable
#' @param sample.2 number of the sample i.e. study variable
#' @param pdf.file path to output pdf file
#'
#' @export
plotFcLogIntensitySingleProtein <- function(data, protein, sample.1, sample.2, pdf.file)
{
  idx <- which(data$accession == protein)
  data <- data[idx,]

  # make (modified sequence, charge) pair unique i.e. use best quants only
  data <- makeModifiedSequenceChargeUnique(data)

  # calculate fc and average intensity
  column.1 <- paste("peptide_abundance_study_variable[", as.character(sample.1), "]", sep="")
  column.2 <- paste("peptide_abundance_study_variable[", as.character(sample.2), "]", sep="")
  data$fc <- calculateFoldChange(data[[column.1]], data[[column.2]])
  data$intensity <- getAverageIntensity(data)
  idx <- complete.cases(data[,c("fc","intensity")])
  data <- data[idx,]

  # define colours
  colours <- c("grey", "red")

  pdf(file=pdf.file)
  plot(data$fc, data$intensity, pch=20, col=colours[data$unique+1], xlab=paste("fold change (sample ", sample.1, " vs ", sample.2, ")", sep=""), ylab="intensity", log="y", main=protein)
  abline(v=0, col = "gray", lty=1)
  abline(v=median(data$fc, na.rm=TRUE), col = "gray", lty=2)
  legend("topright", legend=c("unique", "non-unique"), fill=c("red", "grey"), box.lty=0)
  dev.off()
}

