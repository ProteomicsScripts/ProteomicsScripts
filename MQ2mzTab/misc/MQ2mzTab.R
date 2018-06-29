## This is an R script for the conversion of mzTab to a better readable tsv format
## To install dependencies, run in R:
## install.packages(devtools)
## install.packages(tidyr)

library("tidyr")
library("dplyr")

## Questions ##
# 1. Example folder with maxquant output files?
# 2. are 

# clear entire workspace
rm(list = ls())

# options and parameters
options(digits=10)

input.folder <- 'misc/maxquant_example'

# generate a PEP section from files in `maxquant_folder`.
# The following columns are considered necessary:
## sequence = peptide sequence
## accession 
## peptide abundance_[1 -- n] abundances for each of 1..n study variables
## RT
## mz
## charge

# Each row is one (in mzTab potentially non-unique) detected peptide.
# Q: From which files can we obtain this information in "maxquant_folder"
generatePEP<- function(allPeptidesFile) {
  t = read.table(allPeptidesFile, sep="\t", header=TRUE)
  print(colnames(t))
  print(dim(t))
  mass = t["Mass"]
  # num_study_variables = length(strsplit(as.character(t["Intensities"][[1]]), ";"))
  # study_variables = sprintf("peptide_abundance_study_variable[%s]", seq(1: num_study_variables))
  # split intensities 
  # (they are ";"-seperated in max quant, but mzTab expects different study variables there)
  # TODO: Pre-generate list of "peptide_abundance_study_variable[i]" things

  # maxquant lists matching proteins in a ";" seperated string for each peptide
  # => mzTab expects a single protein accession per listed peptide, 
  # so we need to seperate each row whose "Proteins" column lists multiple proteins
  # into many rows with one protein each.
  t = separate_rows(t, col="Proteins", sep=";", convert=TRUE)
  # seperate(t, col="Intensities", into=c("peptide_abundance_study_variable_1", "peptide_abundance_study_variable_2") sep=";", convert=TRUE)
  print(t["Proteins"])


  sequence = c(0, 0)
  accession = c(0, 0)
  peptide_abundance_study_variable_1 = c(0, 0)
  peptide_abundance_study_variable_2 = c(0, 0)
  rt = c(0, 0)
  mz = c(0, 0)
  charge = c(0, 0)

  # mz vs uncalibrated mz?
  # Proteins = is this comma seperated? mzTab spec says "duplicate rows in case of multiple protein accessions for a peptide"
  frame = data.frame(t["Sequence"], t["Proteins"], t["m.z"], t["Retention.time"], t["Charge"])

  colnames(frame) = c("sequence", "accession", "mz", "rt", "charge")
  return (frame)
}


f = file.path(input.folder, "allPeptides.txt")
pep <- generatePEP(f)
