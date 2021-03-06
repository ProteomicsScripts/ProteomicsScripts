% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mzTabHelper.R
\name{plotMultiplicityFrequency}
\alias{plotMultiplicityFrequency}
\title{Plot (modified sequence, charge) pair multiplicity vs frequency plot.}
\usage{
plotMultiplicityFrequency(data, pdf.file)
}
\arguments{
\item{data}{peptide dataframe}

\item{pdf.file}{path to output pdf file}
}
\description{
Each peptide feature (characterised by a (possibly) modified peptide sequence and a charge state) should ideally occur only once in the analysis.
In other words, peptides of multiplicity 1 should have a very high frequency. The plot below should show a significant spike on the left and can be used as QC of the analysis.
}
