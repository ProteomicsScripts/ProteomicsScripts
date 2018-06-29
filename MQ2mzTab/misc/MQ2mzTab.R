## This is an R script for the conversion of mzTab to a better readable tsv format
## To install dependencies, run in R:
## install.packages(devtools)
## library(devtools)
## install_github("corrplot", username = "MFreidank")

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
  print(t)
  sequence = c(0, 0)
  accession = c(0, 0)
  peptide_abundance_study_variable_1 = c(0, 0)
  peptide_abundance_study_variable_2 = c(0, 0)
  rt = c(0, 0)
  mz = c(0, 0)
  charge = c(0, 0)

  frame = data.frame(sequence, accession, peptide_abundance_study_variable_1, peptide_abundance_study_variable_2, rt, mz, charge)
  colnames(frame) = c("sequence", "accession", "peptide_abundance_study_variable[1]", "peptide_abundance_study_variable[2]", "rt", "mz", "charge")
  return (frame)
}


files <- list.files(path=input.folder, pattern = "allPeptides.txt")
print(files)
pep <- generatePEP(allPeptidesFile=list.files(path=input.folder, pattern="allPeptides.txt"))
print(pep)
