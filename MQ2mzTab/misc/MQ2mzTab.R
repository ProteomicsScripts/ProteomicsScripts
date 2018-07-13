#!/usr/bin/env Rscript
## TODO: Change docu here.
## This is an R script for the conversion of mzTab to a better readable tsv format

library("methods")  # required to avoid a warning from R
library("tidyr")  # allows simple and efficient separation of semi-colon separated "Proteins" column into multiple rows

# clear entire workspace
rm(list = ls())

# options and parameters
options(digits=10)


# ANSI Escape sequences to print in (red) color to terminal screens and 
# to turn color printing off again.
# https://en.wikipedia.org/wiki/ANSI_escape_code
FAIL_COLOR_ANSI = "\033[91m"
SUCCESS_COLOR_ANSI = '\033[32m'
END_COLOR_ANSI = "\033[0m"


#' Assert that a given folder exists and contains all necessary MaxQuant result files.

#' Required files for processing are: allPeptides.txt 

#' @param folder Character, path to a folder on disk that should contain MaxQuant results.
#' @param relevant_files Character List, Files that are necessary for conversion to \code{mzTab}.
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


#' Determine type of analysis used to obtain our MaxQuant input files.
#' 
#' Supported types of analyses: \code{"TMT"}, \code{"Labeled"}
#' A TMT analysis is expected to have a column with name \code{"Reporter.intensity"}.
#' A Labeled analysis is expected to have a column with name \code{"Intensity.H"} 
#' and a column with name \code{Intensity.L}.
#' 
#' @param column_names List of Characters, Names of all columns of MaxQuant's \code{allPeptides.txt}.

analysisType <- function(column_names) {
    if (any(grepl("Reporter.intensity", column_names))) {
        return("TMT")
    } else if (any(grepl("Intensity.L", column_names)) && 
               any(grepl("Intensity.H", column_names)))
    {
        return("Labeled")
    } else if (any(grepl("Intensities", column_names))) {
        return ("Labelfree")

    }
    message = "Unknown type of analysis.\nTMT analyses are expected to have a column named 'Reporter intensity'.\nLabeled analyses are expected to have at least a column named 'Intensity.L' and a column named 'Intensity.H'."
    stop(sprintf("%s%s%s", FAIL_COLOR_ANSI, message, END_COLOR_ANSI))
}

# Each row is one (in mzTab potentially non-unique) detected peptide.
#' Generate a \code{mzTab} PEP section from the given \code{max_quant_peptides}. 
#'
#' \code{max_quant_peptides} are expected to have a column-format like in MaxQuant's 
#' \code{allPeptides.txt}. 
#'
#' @param max_quant_peptides data.frame, Peptides parsed from MaxQuant's \code{allPeptides.txt}, possibly with semicolon separated \code{"Proteins"} column for each peptide. 
generatePEP<- function(max_quant_peptides) {

  # Determine type of analysis used.
  column_names = sort(colnames(max_quant_peptides))
  analysis = analysisType(column_names=column_names)

  # use tidyr magic to split a row with semicolon-separated column "Proteins" into 
  # multiple rows that each have a single protein 
  # (and the same values in all remaining columns).
  max_quant_peptides = separate_rows(max_quant_peptides, col="Proteins", sep=";")

  # NOTE: The code below this is quite lengthy and can probably be improved.
  mztab_column_names = c("PEH", "sequence", "accession", "unique", "database",	
    "database_version", "search_engine", "best_search_engine_score[1]", 
    "search_engine_score[1]_ms_run[1]", "modifications", "retention_time",
    "retention_time_window",	
    "charge",	"mass_to_charge", "spectra_ref"
  )
  print(colnames(max_quant_peptides))


  # NOTE: modifications missing 
  if (analysis == "Labeled") {
      # Results stem from TMT analysis.
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
                          max_quant_peptides["Retention.time"] * 60,  # minutes to seconds
                          max_quant_peptides["Retention.length"] * 60,  # minutes to seconds
                          max_quant_peptides["Charge"], 
                          max_quant_peptides["Uncalibrated.m.z"],
                          nulls,
                          max_quant_peptides["Intensity.L"], nulls, nulls,
                          max_quant_peptides["Intensity.H"], nulls, nulls,
                          max_quant_peptides["Intensity.M"], nulls, nulls,)
      } else {
          df = data.frame(rep("PEP", nrow(max_quant_peptides)), 
                          max_quant_peptides["Sequence"], 
                          max_quant_peptides["Proteins"],
                          nulls, nulls, nulls, nulls, 
                          max_quant_peptides["Score"], nulls, 
                          max_quant_peptides["Modifications"], 
                          max_quant_peptides["Retention.time"] * 60,
                          max_quant_peptides["Retention.length"] * 60,
                          max_quant_peptides["Charge"], 
                          max_quant_peptides["Uncalibrated.m.z"],
                          nulls,
                          max_quant_peptides["Intensity.L"], nulls, nulls,
                          max_quant_peptides["Intensity.H"], nulls, nulls)
      }
  }
  else if (analysis == "TMT" || analysis == "Labelfree") {
      # XXX: handle tmt, reporter.intensity seems to be how peptide_abundance_study_variables are represented here.
          mztab_column_names = c(mztab_column_names, 
                             "peptide_abundance_study_variable[1]", 
                             "peptide_abundance_stdev_study_variable[1]", 
                             "peptide_abundance_std_error_study_variable[1]")
          df = data.frame(rep("PEP", nrow(max_quant_peptides)), 
                          max_quant_peptides["Sequence"], 
                          max_quant_peptides["Proteins"],
                          nulls, nulls, nulls, nulls, 
                          max_quant_peptides["Score"], nulls, 
                          max_quant_peptides["Modifications"], 
                          max_quant_peptides["Retention.time"] * 60,  # minutes to seconds
                          max_quant_peptides["Retention.length"] * 60,  # minutes to seconds
                          max_quant_peptides["Charge"],  
                          max_quant_peptides["Uncalibrated.m.z"],
                          nulls, null, nulls, nulls,
                          )
  }

  # Overwrite column names with correct ones.
  colnames(df) <- mztab_column_names

  # replace NA with "null" (as this is the 'missing information marker' for mzTab)
  df[is.na(df)] <- "null"
  return (df)
}

#' Determine an output filename from the `Raw file` column of the MaxQuant inputs.
#'
#' If there is a single unique filename in the `Raw file` column of the MaxQuant inputs, 
#' the output filename is exactly that unique filename with suffix \code{".mzTab"}. 
#' If there are multiple filenames in the `Raw file` column of the MaxQuant inputs, 
#' the output filename is \code{RAWFILE_etal.mzTab} where \code{RAWFILE} is 
#' the name of the first filename in the `Raw file` column of the MaxQuant inputs.
#'
#' @param input_files Character List, Filenames present in the `Raw file` column of \code{allPeptides.txt} of MaxQuant.
outputFilename<- function(input_files) {
    single_input_file = length(unique(input_files)) == 1
    first_input_file = input_files[[1]]
    if (single_input_file) {
        return (paste(first_input_file, ".mzTab", sep=""))
    }
    return (paste(first_input_file, "_etal.mzTab", sep=""))
}

#' Write a header for a resulting \code{mzTab} file.
#'
#' Header has information about information source (\code{uri}) and 
#' \code{mzTab} metadata.
#'
#' @param uri Character, Name of input folder containing MaxQuant results.
#' @param maintainer Character, Description of current maintainer of this tool.
mzTabHeader <- function(uri, maintainer="Moritz Freidank, freidankm@gmail.com") {
    return (c("MTD\tmzTab-version\t1.0.0", "MTD\tmzTab-mode\tSummary",
               "MTD\tmzTab-type\tQuantification", "MTD\tdescription\tGenerated using MQ2mzTab.R from MaxQuant Output",
               "MTD\tmaintainer\tMoritz Freidank, freidankm@gmail.com", 
               "MTD\tpeptide_search_engine_score[1]\tnull",
               "MTD\tpsm_search_engine_score[1]\tnull",
               paste("MTD\turi[1]\t", uri)))
}


#  Parse folder name from command line and validate it. {{{ # 
command_line_arguments <- commandArgs(trailingOnly = TRUE)

if (length(command_line_arguments) != 1) {
    message = "Invalid amount of arguments!\nUsage: MQ2mzTab.R MAXQUANT_OUTPUT_FOLDER"
    stop(sprintf("%s%s%s", FAIL_COLOR_ANSI, message, END_COLOR_ANSI))    
}

input.folder <- commandArgs(trailingOnly = TRUE)[[1]]
checkMaxQuantFolder(input.folder)
#  }}} Parse folder name from command line # 


#  Generate PEP section {{{ # 

f = file.path(input.folder, "allPeptides.txt")
max_quant_peptides = read.table(f, sep="\t", header=TRUE, stringsAsFactors=FALSE, na.strings=c("", "NA", " ", "  "))

output_filename <- outputFilename(input_files=max_quant_peptides$Raw.file)

pep_section <- generatePEP(max_quant_peptides)
#  }}} Generate PEP section # 


    
#  Write mzTab header to file {{{ # 
output_file <- file(output_filename, open="wt")
on.exit(close(output_file))

# uri is the source used to generate a given mzTab file. 
# In our context here it is 'MAXQUANT_RESULTS_FOLDER/allPeptides.txt'
uri = normalizePath(input.folder, "allPeptides.txt")

cat(mzTabHeader(uri=uri), sep="\n", file=output_file)
cat("\n", file=output_file)
#  }}} Write mzTab header to file # 

#  Append PEP section. {{{ # 

write.table(pep_section, 
            file=output_file, 
            quote=FALSE, 
            sep="\t", 
            row.names=FALSE)
#  }}} Append PEP section. # 

cat(sprintf("%sDone. Your mzTab file is available at: %s%s\n", SUCCESS_COLOR_ANSI, output_filename, END_COLOR_ANSI))
