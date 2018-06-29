## This is an R script for the conversion of mzTab to a better readable tsv format
## To install dependencies, run in R:
## install.packages(devtools)
## install.packages(tidyr)

library("tidyr")
library("dplyr")

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
# NOTE: Check how modifications are represented in maxquant => 
generatePEP<- function(allPeptidesFile) {
  t = read.table(allPeptidesFile, sep="\t", header=TRUE)
  column_names = sort(colnames(t))
  t = separate_rows(t, col="Proteins", sep=";", convert=TRUE)
  print(t["Raw.file"])

  # Check which type of analysis this is.
  is_tmt = any(grepl("Reporter.intensity", column_names))
  is_labeled = any(grepl("Intensity.L", column_names)) && any(grepl("Intensity.H", column_names))

  if (is_labeled) {
      mztab_column_names = c("sequence", "accession", "mass", 
                             "best_search_engine_score[1]",
                             "mz", "rt", "charge",
                             "modifications",
                             "peptide_abundance_study_variables[1]", 
                             "peptide_abundance_study_variables[2]")

      if (any(grepl("Intensity.M", column_names))) {
          mztab_column_names = c(mztab_column_names, 
                                 "peptide_abundance_study_variables[3]")
          df = data.frame(t["Sequence"], t["Proteins"], t["Mass"], 
                          t["Score"], t["m.z"], t["Retention.time"], t["Charge"], 
                          t["Modifications"], 
                          t["Intensity.L"], t["Intensity.H"], t["Intensity.M"])
      } else {
          df = data.frame(t["Sequence"], t["Proteins"], t["Mass"], 
                          t["Score"], t["m.z"], t["Retention.time"], t["Charge"], 
                          t["Modifications"], t["Intensity.L"], t["Intensity.H"])
      }
  } else if (is_tmt) {
      # XXX: handle tmt, reporter.intensity seems to be how peptide_abundance_study_variables are represented here.
      stop("It seems as if your input is a TMT analysis. Currently only labeled analyses are supported!")
  } else {
      stop("Unsupported type of analysis.\nExpected TMT analysis (with at least one column starting with 'Reporter.intensity') or labeled analysis (with at least columns 'Intensity H', 'Intensity L').")
  }

  print(is_tmt)
  print(is_labeled)
  # maxquant lists matching proteins in a ";" seperated string for each peptide
  # => mzTab expects a single protein accession per listed peptide, 
  # so we need to seperate each row whose "Proteins" column lists multiple proteins
  # into many rows with one protein each.
  # NOTE: mz vs uncalibrated mz?

  df <- cbind(rep("PEP", nrow(df)), df)
  colnames(df)[1] <- "PEH"
  return (df)
}


f = file.path(input.folder, "allPeptides.txt")
pep <- generatePEP(f)

output_file <- file("misc/maxquant_example.mzTab", open="wt")
on.exit(close(output_file))

header = c("MTD\tmzTab-version\t1.0.0", "MTD\tmzTab-mode\tSummary",
           "MTD\tmzTab-type\tQuantification", "MTD\tdescription\tGenerated automatically from MaxQuant Output",
           "MTD\tpeptide_search_engine_score[1]\tnull",
           "MTD\tpsm_search_engine_score[1]\tnull",
           paste("MTD\turi[1]\t", normalizePath(input.folder, "allPeptides.txt")),
           "\n")
           # NOTE: modifications missing and maybe there is an equivalent of `ms_run[1]-location`?
    

# Write mzTab header to file
cat(header, sep="\n", file=output_file)

write.table(pep, file=output_file, quote=FALSE, sep="\t", append=T, row.names=FALSE)
