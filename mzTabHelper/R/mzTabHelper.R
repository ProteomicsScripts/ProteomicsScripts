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



