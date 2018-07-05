## This is an R script for the conversion of mzTab to a better readable tsv format
## To install dependencies, run in R:
## install.packages(devtools)
## install.packages(tidyr)

library("tidyr")

# clear entire workspace
rm(list = ls())

# options and parameters
options(digits=10)


# ANSI Escape sequences to print in (red) color to terminal screens and 
# to turn color printing off again.
# https://en.wikipedia.org/wiki/ANSI_escape_code
FAIL_COLOR_ANSI = "\033[91m"
END_COLOR_ANSI = "\033[0m"


#' Assert that a given folder exists and contains all necessary MaxQuant result files.

#' Required files for processing are: allPeptides.txt 

#' @param folder Character, path to a folder on disk that should contain MaxQuant results.
#' @param relevant_files List of files that are necessary for conversion to \code{mzTab}.
checkMaxQuantFolder <- function(folder, relevant_files=c("allPeptides.txt")) {
    if (!dir.exists(file.path(folder))) {
        message = sprintf("Folder at path '%s' not found.", input.folder)
        stop(sprintf("%s%s%s", FAIL_COLOR_ANSI, message, END_COLOR_ANSI))    
    }
    for (file in relevant_files) {
        if (!file.exists(file.path(folder, file))) {
            message = sprintf("Found folder at path '%s' but it does not contain necessary file '%s'.", folder, file
            )
            stop(sprintf("%s%s%s", FAIL_COLOR_ANSI, message, END_COLOR_ANSI))
        }
    }
}

separateRows <- function(dataframe, col, sep) {
    rest_of_dataframe = dataframe[, -which(names(dataframe) == col)]
    # print(rep(rest_of_dataframe, 2))
    return (dataframe)
}

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
generatePEP<- function(max_quant_peptides) {
  column_names = sort(colnames(max_quant_peptides))

  # max_quant_peptides = separate_rows(max_quant_peptides, 
  #                                    col="Proteins", 
  #                                    sep=";", 
  #                                    convert=TRUE)

  # rest_of_dataframe = max_quant_peptides[, -which(names(max_quant_peptides) == "Proteins")]
  # proteins_split = strsplit(max_quant_peptides$Proteins, ";")
  # duplicated = rep(rest_of_dataframe$, sapply(proteins_split, length))
  # max_quant_peptides = data.frame(unlist(proteins_split), rest_of_dataframe)
  max_quant_peptides = separateRows(max_quant_peptides, col="Proteins", sep=";")

  # Check which type of analysis this is.
  is_tmt = any(grepl("Reporter.intensity", column_names))
  is_labeled = any(grepl("Intensity.L", column_names)) && any(grepl("Intensity.H", column_names))
  mztab_column_names = c(
    "PEH", "sequence", "accession", "unique", "database",	
    "database_version", "search_engine", "best_search_engine_score[1]", 
    "search_engine_score[1]_ms_run[1]", "modifications", "retention_time",
    "retention_time_window",	
    "charge",	"mass_to_charge	spectra_ref")


  if (is_labeled) {
      mztab_column_names = c(mztab_column_names, 
                             "peptide_abundance_study_variable[1]", 
                             "peptide_abundance_stdev_study_variable[1]", 
                             "peptide_abundance_std_error_study_variable[1]",
                             "peptide_abundance_study_variable[2]", 
                             "peptide_abundance_stdev_study_variable[2]", 
                             "peptide_abundance_std_error_study_variable[2]")

      nulls = rep("null", nrow(max_quant_peptides))  # a column with only null values
      if (any(grepl("Intensity.M", column_names))) {
          mztab_column_names = c(mztab_column_names, 
                             "peptide_abundance_study_variable[3]", 
                             "peptide_abundance_stdev_study_variable[3]", 
                             "peptide_abundance_std_error_study_variable[3]")
          df = data.frame(rep("PEP", nrow(max_quant_peptides)), 
                          max_quant_peptides["Sequence"], 
                          max_quant_peptides["Proteins"],
                          nulls, nulls, nulls, nulls, 
                          max_quant_peptides["Score"], nulls, 
                          max_quant_peptides["Modifications"], 
                          max_quant_peptides["Retention.time"],
                          nulls, max_quant_peptides["Charge"], nulls, 
                          max_quant_peptides["Intensity.L"], nulls, nulls,
                          max_quant_peptides["Intensity.H"], nulls, nulls,
                          max_quant_peptides["Intensity.M"], nulls, nulls)
      } else {
          df = data.frame(rep("PEP", nrow(max_quant_peptides)), 
                          max_quant_peptides["Sequence"], 
                          max_quant_peptides["Proteins"],
                          nulls, nulls, nulls, nulls, 
                          max_quant_peptides["Score"],
                          nulls, 
                          max_quant_peptides["Modifications"], 
                          max_quant_peptides["Retention.time"],
                          nulls, 
                          max_quant_peptides["Charge"], 
                          nulls, 
                          max_quant_peptides["Intensity.L"], nulls, nulls,
                          max_quant_peptides["Intensity.H"], nulls, nulls)
      }
  } else if (is_tmt) {
      # XXX: handle tmt, reporter.intensity seems to be how peptide_abundance_study_variables are represented here.
      stop("It seems as if your input is a TMT analysis. Currently only labeled analyses are supported!")
  } else {
      stop("Unsupported type of analysis.\nExpected TMT analysis (with at least one column starting with 'Reporter intensity') or labeled analysis (with at least columns 'Intensity H', 'Intensity L').")
  }

  print(is_tmt)
  print(is_labeled)
  # maxquant lists matching proteins in a ";" seperated string for each peptide
  # => mzTab expects a single protein accession per listed peptide, 
  # so we need to seperate each row whose "Proteins" column lists multiple proteins
  # into many rows with one protein each.
  # NOTE: mz vs uncalibrated mz?
  colnames(df) <- mztab_column_names

  # replace NA with "null" (as this is the 'missing information marker' for mzTab)
  df[is.na(df)] <- "null"
  return (df)
}

outputFilename<- function(input_files) {
    single_input_file = length(unique(input_files)) == 1
    first_input_file = input_files[[1]]
    if (single_input_file) {
        return (paste(first_input_file, ".mzTab", sep=""))
    }
    return (paste(first_input_file, "_etal.mzTab", sep=""))
}

mzTabHeader <- function(uri, maintainer="Moritz Freidank, freidankm@gmail.com") {
    return (c("MTD\tmzTab-version\t1.0.0", "MTD\tmzTab-mode\tSummary",
               "MTD\tmzTab-type\tQuantification", "MTD\tdescription\tGenerated using MQ2mzTab.R from MaxQuant Output",
               "MTD\tmaintainer\tMoritz Freidank, freidankm@gmail.com", 
               "MTD\tpeptide_search_engine_score[1]\tnull",
               "MTD\tpsm_search_engine_score[1]\tnull",
               paste("MTD\turi[1]\t", uri)))

}


command_line_arguments <- commandArgs(trailingOnly = TRUE)

if (length(command_line_arguments) != 1) {
    message = "Invalid amount of arguments!\nUsage: MQ2mzTab.R MAXQUANT_OUTPUT_FOLDER"
    stop(sprintf("%s%s%s", FAIL_COLOR_ANSI, message, END_COLOR_ANSI))    
}

input.folder <- commandArgs(trailingOnly = TRUE)[[1]]
checkMaxQuantFolder(input.folder)


f = file.path(input.folder, "allPeptides.txt")
max_quant_peptides = read.table(f, sep="\t", header=TRUE, stringsAsFactors=FALSE)

output_filename <- outputFilename(input_files=max_quant_peptides$Raw.file)

pep_section <- generatePEP(max_quant_peptides)

output_file <- file(output_filename, open="wt")
on.exit(close(output_file))

           # NOTE: modifications missing and maybe there is an equivalent of `ms_run[1]-location`?
    
#  Write mzTab header to file {{{ # 

# uri is the source used to generate a given mzTab file. 
# In our context here it is 'MAXQUANT_RESULTS_FOLDER/allPeptides.txt'
uri = normalizePath(input.folder, "allPeptides.txt")

cat(mzTabHeader(uri=uri), sep="\n", file=output_file)
cat("\n", file=output_file)
#  }}} Write mzTab header to file # 

# Append PEP section.
write.table(pep_section, 
            file=output_file, 
            quote=FALSE, 
            sep="\t", 
            row.names=FALSE)
