## This R script generates a set of basic plots and tables which summarise the PEP section of mzTab files.
##
## It is divided into two parts:
## (1) definition of global options and parameters
## (2) main part i.e. generation of plots and tables

# package with all helper functions
# https://github.com/ProteomicsScripts/ProteomicsScripts/mzTabHelper
library(mzTabHelper)

# clear entire workspace
rm(list = ls())

####
## (1) definition of global options and parameters
####

#input.file <- 'analysis.mzTab'
input.file <- 'example_7.mzTab'

# maximum number of digits
options(digits=10)

# switch off warnings
options(warn=-1)

# fold change cutoff, i.e. infinite fc values are mapped to +/-FcCutoff
FcCutoff <- 8

# study variable labels
# The PEP section of the mzTab file contains columns peptide_abundance_study_variable[*]. Each column can be labelled with a string such as "control" or "treated".
# If the number of study variable columns is not equal to the length of the labels vector, then the labels vector is ignored.
labels.of.study.variables <- rep(c("H", "BL", "fu48"), times = 18)

# Biognosys iRT spike-in peptides
#peptides.of.interest <- c("LGGNEQVTR", "GAGSSEPVTGLDAK", "VEATFGVDESNAK", "YILAGVENSK", "TPVISGGPYEYR", "TPVITGAPYEYR", "DGLDAASYYAPVR", "ADVTPADFSEWSK", "GTFIIDPGGVIR", "GTFIIDPAAVIR", "LFLQFGAQGSPFLK")

# Pierce spike-in peptides
peptides.of.interest <- c("SSAAPPPPPR", "GISNEGQNASIK", "HVLTSIGEK", "DIPVPKPK", "IGDYAGIK", "TASEFDSAIAQDK", "SAAGAFGPELSR", "ELGQSGVDTYLQTK", "GLILVGGYGTR", "GILFVGSGVSGGEEGAR", "SFANQPLEVVYSK", "LTILEELR", "NGFILDGFPR", "ELASGLSFPVGFK", "LSSEAPALFQFDLK")

# some random human peptides
#peptides.of.interest <- c("LSLMYAR", "EQCCYNCGKPGHLAR", "LSAIYGGTYMLNKPVDDIIMENGKVVGVK", "MVQEAEKYKAEDEKQR", "TVPFCSTFAAFFTR", "GNFGGSFAGSFGGAGGHAPGVAR", "LGWDPKPGEGHLDALLR")

# proteins of interest
proteins.of.interest <- c("P46783", "P12270")





####
## (2) main part
####

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

# create quant summary statistics
# How many peptides have finite, zero and NaN abundance in each sample?
#stats.quants <- getQuantSummary(peptide.data)

# create mod summary statistics
stats.mods <- getModsSummary(peptide.data)

# total number of quantified and identified peptides
n.peptides <- dim(peptide.data)[1]
n.peptides.any.nan <- length(which(apply(getPeptideQuants(peptide.data), 1, function(row) any(is.nan(unlist(row))))))
n.peptides.any.zero <- length(which(apply(getPeptideQuants(peptide.data), 1, function(row) any(unlist(row) == 0))))
peptide.data.identified <- peptide.data[which(!is.na(peptide.data$sequence)),]
n.peptides.identified <- dim(peptide.data.identified)[1]
n.peptides.identified.modified.unique <- length(unique(peptide.data.identified$opt_global_modified_sequence))
n.peptides.identified.stripped.unique <- length(unique(peptide.data.identified$sequence))

# plot elution time distribution
if (!isEmpty(peptide.data$retention_time_window))
{
  plotElutionTimeDistribution(peptide.data, "plot_DistributionElutionTime.pdf")
}

# plot frequency of peptide quants
if (numberOfStudyVariables(peptide.data) >= 2) {
  plotQuantFrequency(getPeptideQuants(peptide.data), "plot_QuantFrequency.pdf")
}

# plot charge distribution
if (!isEmpty(peptide.data$charge))
{
  plotChargeDistribution(peptide.data, "plot_ChargeDistribution.pdf")
}

# plot frequency of multiply quantified sequences
if (!isEmpty(peptide.data$opt_global_modified_sequence) && !isEmpty(peptide.data$charge))
{
  plotMultiplicityFrequency(peptide.data, "plot_MultiplicityFrequency.pdf")
}

# plot frequency of peptides per protein
if (!isEmpty(peptide.data$accession) && !isEmpty(peptide.data$sequence) && !isEmpty(peptide.data$opt_global_modified_sequence) && !isEmpty(peptide.data$charge))
{
  plotPeptidesPerProtein(makeModifiedSequenceChargeUnique(peptide.data), "plot_PeptidesPerProteinFrequency.pdf")
}

# extract peptides and proteins of interest
interest.peptides.matches <- findPeptidesOfInterest(peptide.data)
interest.proteins.matches <- findProteinsOfInterest(peptide.data)

# plot abundances of peptides and proteins of intrest
plotPeptidesOfInterest(peptide.data, "plot_PeptidesOfInterest.pdf")
plotProteinsOfInterest(peptide.data, "plot_ProteinsOfInterest.pdf")

# pre-define variables which will be called from LaTeX
# Needed even with IfFileExists()

median.abundance.1 <- 1
median.abundance.2 <- 1
median.abundance.3 <- 1

median.fc.12 <- 0
median.fc.13 <- 0
median.fc.23 <- 0

sd.fc.12 <- 0
sd.fc.13 <- 0
sd.fc.23 <- 0

# Kendrick plot
plotKendrick((peptide.data$mass_to_charge - 1.00784) * peptide.data$charge, "plot_Kendrick.pdf")

# plot peptide abundance distributions
if (studyVariableExists(peptide.data,1)) {
  abundances <- peptide.data$"peptide_abundance_study_variable[1]"
  abundances <- abundances[complete.cases(abundances)]
  median.abundance.1 <- median(abundances, na.rm=TRUE)
  plotDistribution(log10(abundances), expression('log'[10]*' intensity'), "plot_DistributionIntensity_1.pdf")
}
if (studyVariableExists(peptide.data,2)) {
  abundances <- peptide.data$"peptide_abundance_study_variable[2]"
  abundances <- abundances[complete.cases(abundances)]
  median.abundance.2 <- median(abundances, na.rm=TRUE)
  plotDistribution(log10(abundances), expression('log'[10]*' intensity'), "plot_DistributionIntensity_2.pdf")
}
if (studyVariableExists(peptide.data,3)) {
  abundances <- peptide.data$"peptide_abundance_study_variable[3]"
  abundances <- abundances[complete.cases(abundances)]
  median.abundance.3 <- median(abundances, na.rm=TRUE)
  plotDistribution(log10(abundances), expression('log'[10]*' intensity'), "plot_DistributionIntensity_3.pdf")
}

# plot fold change distributions and scatter plots
if (studyVariableExists(peptide.data,1) && studyVariableExists(peptide.data,2)) {
  a <- peptide.data$"peptide_abundance_study_variable[1]"
  b <- peptide.data$"peptide_abundance_study_variable[2]"
  fc <- calculateFoldChange(a, b)
  a[is.na(a)] <- 0
  b[is.na(b)] <- 0
  intensity <- (a + b)/2
  median.fc.12 <- median(fc, na.rm=TRUE)
  sd.fc.12 <- sd(fc, na.rm=TRUE)
  plotFcLogIntensity(fc, intensity, "fold change", "plot_FoldChangeLogIntensity_12.pdf")
  plotDistribution(fc, "fold change", "plot_DistributionFoldChange_12.pdf")

}
if (studyVariableExists(peptide.data,1) && studyVariableExists(peptide.data,2) && !(isEmpty(peptide.data$opt_global_modified_sequence)))
{
  plotDeltaFcLogIntensity(peptide.data, 1, 2, "plot_DeltaFoldChangeLogIntensity_12.pdf")
}
if (studyVariableExists(peptide.data,1) && studyVariableExists(peptide.data,3)) {
  a <- peptide.data$"peptide_abundance_study_variable[1]"
  b <- peptide.data$"peptide_abundance_study_variable[3]"
  fc <- calculateFoldChange(a, b)
  a[is.na(a)] <- 0
  b[is.na(b)] <- 0
  intensity <- (a + b)/2
  median.fc.13 <- median(fc, na.rm=TRUE)
  sd.fc.13 <- sd(fc, na.rm=TRUE)
  plotFcLogIntensity(fc, intensity, "fold change", "plot_FoldChangeLogIntensity_13.pdf")
  plotDistribution(fc, "fold change", "plot_DistributionFoldChange_13.pdf")
}
if (studyVariableExists(peptide.data,1) && studyVariableExists(peptide.data,3) && !(isEmpty(peptide.data$opt_global_modified_sequence)))
{
  plotDeltaFcLogIntensity(peptide.data, 1, 3, "plot_DeltaFoldChangeLogIntensity_13.pdf")
}
if (studyVariableExists(peptide.data,2) && studyVariableExists(peptide.data,3)) {
  a <- peptide.data$"peptide_abundance_study_variable[2]"
  b <- peptide.data$"peptide_abundance_study_variable[3]"
  fc <- calculateFoldChange(a, b)
  a[is.na(a)] <- 0
  b[is.na(b)] <- 0
  intensity <- (a + b)/2
  median.fc.23 <- median(fc, na.rm=TRUE)
  sd.fc.23 <- sd(fc, na.rm=TRUE)
  plotFcLogIntensity(fc, intensity, "fold change", "plot_FoldChangeLogIntensity_23.pdf")
  plotDistribution(fc, "fold change", "plot_DistributionFoldChange_23.pdf")
}
if (studyVariableExists(peptide.data,2) && studyVariableExists(peptide.data,3) && !(isEmpty(peptide.data$opt_global_modified_sequence)))
{
  plotDeltaFcLogIntensity(peptide.data, 2, 3, "plot_DeltaFoldChangeLogIntensity_23.pdf")
}

# plot correlation matrix of peptide abundances
corr.min <- 1
corr.median <- 1
corr.max <- 1
if (numberOfStudyVariables(peptide.data) >= 3) {
  corr <- plotCorrelations(data = peptide.data, pdf.file = "plot_Correlations.pdf")
  corr.min <- min(corr)
  corr.median <- median(corr)
  corr.max <- max(corr)
}

# plot boxplot of peptide abundances
if (numberOfStudyVariables(peptide.data) >= 3) {
  plotBoxplot(peptide.data, "plot_Boxplot.pdf")
}

# start of Principal Component Analysis
# (Even if we do not run the PCA, we generate these three non-empty tables in order to prevent LaTeX from crashing.)
important.peptides.principal.component.1 <- data.frame(c(42))
important.peptides.principal.component.2 <- data.frame(c(42))
important.peptides.principal.component.3 <- data.frame(c(42))
if (numberOfStudyVariables(peptide.data) >= 3) {

  ## simple ggplot2 version of PCA plot
  #plotPCA(peptide.data, "plot_PCA.pdf")

  pca <- getPCA(peptide.data)

  plotPCAscatter(pca, "plot_PCA_scatter.pdf")

  plotPCAcomponents(pca, "plot_PCA_components.pdf")

  plotPCAeigenvector(pca, peptide.data, 1, "plot_PCA_eigenvector1st.pdf")
  plotPCAeigenvector(pca, peptide.data, 2, "plot_PCA_eigenvector2nd.pdf")
  plotPCAeigenvector(pca, peptide.data, 3, "plot_PCA_eigenvector3rd.pdf")


  # Note that getPCAeigenvector() returns the row indices with respect to the complete cases.
  # Since the complete cases appear first in the PEP section, the row indices are the same as for the entire peptide data.
  # But let's play it save.
  idx.complete <- which(complete.cases(getPeptideQuants(peptide.data)))

  idx.1 <- idx.complete[getPCAeigenvector(pca, 1)]
  idx.2 <- idx.complete[getPCAeigenvector(pca, 2)]
  idx.3 <- idx.complete[getPCAeigenvector(pca, 3)]

  # add column with row index
  peptide.data$'row index' <- rownames(peptide.data)

  retain.columns=c("row index", "opt_global_modified_sequence", "accession", "charge", "retention_time", "mass_to_charge")
  new.column.names=c("row index", "modified sequence", "accession", "charge", "retention time", "m/z")

  important.peptides.principal.component.1 <- peptide.data[idx.1, retain.columns]
  important.peptides.principal.component.2 <- peptide.data[idx.2, retain.columns]
  important.peptides.principal.component.3 <- peptide.data[idx.3, retain.columns]

  colnames(important.peptides.principal.component.1) <- new.column.names
  colnames(important.peptides.principal.component.2) <- new.column.names
  colnames(important.peptides.principal.component.3) <- new.column.names

  # reduce sequence length
  important.peptides.principal.component.1$'modified sequence' <- unlist(lapply(important.peptides.principal.component.1$'modified sequence', cutSequence))
  important.peptides.principal.component.2$'modified sequence' <- unlist(lapply(important.peptides.principal.component.2$'modified sequence', cutSequence))
  important.peptides.principal.component.3$'modified sequence' <- unlist(lapply(important.peptides.principal.component.3$'modified sequence', cutSequence))

}
# end of Principal Component Analysis


# plot fc vs log intensity for all proteins of interest
if (!isEmpty(peptide.data$accession) && !isEmpty(peptide.data$unique) && !isEmpty(peptide.data$opt_global_modified_sequence) && !isEmpty(peptide.data$charge))
{
  for (p in 1:length(proteins.of.interest))
  {
    if ((length(which(peptide.data$accession == proteins.of.interest[p])) > 0) && (numberOfStudyVariables(peptide.data) > 1))
    {
      pdf.file <- paste("plot_ProteinsOfInterest_", as.character(p), ".pdf", sep="")
      plotFcLogIntensitySingleProtein(peptide.data, proteins.of.interest[p], 1, 2, pdf.file)
    }
  }
}

# plot fc vs log intensity for the proteins with the most quantified peptides
if (!isEmpty(peptide.data$accession) && !isEmpty(peptide.data$unique) && !isEmpty(peptide.data$opt_global_modified_sequence) && !isEmpty(peptide.data$charge))
{
  # remove duplicate, low-intensity peptide quantifications
  data <- makeModifiedSequenceChargeUnique(peptide.data)

  # count number of peptides per protein
  frequency.table <- data.frame(table(data$accession))
  colnames(frequency.table) <- c("accession","frequency")
  data <- merge(data, frequency.table, by="accession")

  # order by frequency and accession
  data <- data[order(data$accession),]
  data <- data[order(data$frequency, decreasing=TRUE),]

  accessions <- unique(data$accession)

  # plot fc vs log intensity for the 9 best quantified proteins
  for (p in 1:9)
  {
    if ((length(which(peptide.data$accession == accessions[p])) > 0) && (numberOfStudyVariables(peptide.data) > 1))
    {
      pdf.file <- paste("plot_ProteinsManyPeptides_", as.character(p), ".pdf", sep="")
      plotFcLogIntensitySingleProtein(peptide.data, accessions[p], 1, 2, pdf.file)
    }
  }
}

# plot retention time shift distribution
plotRetentionTimeShiftDistribution(peptide.data, "plot_RetentionTimeShift.pdf")
