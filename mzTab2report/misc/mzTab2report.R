## This is an R script for the conversion of mzTab to a better readable tsv format
## To install dependencies, run in R:
## install.packages("devtools")
## library("devtools")
## install_github("taiyun/corrplot")

library("corrplot")

# clear entire workspace
rm(list = ls())

# options and parameters
options(digits=10)
FcCutoff <- 8    # fold change cutoff, i.e. infinite fc values are mapped to +/-FcCutoff

# input.file <- 'example_3.mzTab'
input.file <- 'misc/Galaxy150_copy2.mzTab'

# find start of the section
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

# count the occurences of character c in string s
countOccurrences <- function(char,s) {
  s2 <- gsub(char,"",s)
  return (nchar(s) - nchar(s2))
}

# check that the protein accession is of the format *|*|*
# Note that NA returns TRUE.
checkAccessionFormat <- function(accession) {
  if (is.na(accession)) {
    return (TRUE)
  }
  n <- length(accession)
  count <- countOccurrences("[|]",accession)
  m <- length(which(count==2))
  return (n==m)
}

# Extracts the second entry from a string of the form *|*|*.
getAccession <- function(string) {
  if (is.na(string)) {
    return (NA)
  }
  return (unlist(strsplit(string, "[|]"))[2])
}

# Extracts the third entry from a string of the form *|*|*.
getGene <- function(string) {
  if (is.na(string)) {
    return (NA)
  }
  return (unlist(strsplit(string, "[|]"))[3])
}

# read the PEP section of an mzTab file
readMzTabPEP <- function(file) {
  
  # find start of the PEP section
  first.row <- startSection(file, "PEH")
  
  # read entire mzTab
  data <- read.table(file, sep="\t", skip=first.row-1, fill=TRUE, header=TRUE, quote="", na.strings=c("null","NA"), stringsAsFactors=FALSE, check.names=FALSE)
  
  # extract PEP data
  peptide.data <- data[which(data[,1]=="PEP"),]
  peptide.data$PEH <- NULL
  
  # In case the accession column is of the format *|*|*, we split this column into an accession and a gene column.
  if (all(sapply(peptide.data$accession, checkAccessionFormat))) {
    peptide.data$gene <- sapply(peptide.data$accession, getGene)
    peptide.data$accession <- sapply(peptide.data$accession, getAccession)
  }
  
  return (peptide.data)
}

# check if a specific "peptide_abundance_study_variable[n]" column exists
abundance.exists <- function(data, n) {
  column <- paste("peptide_abundance_study_variable[",as.character(n),"]",sep="")
  return (column %in% colnames(data))
}

# determine fold changes and map to finite numbers
# fc = log2(abundances1/abundances2)
calculateFoldChange <- function(abundances1, abundances2) {
  offset <- 1e-10       # avoids devisions by zero
  max.fc <- FcCutoff    # map knock-out fold changes to finite values 
  
  # calculate fold changes
  abundances1 <- abundances1 + offset
  abundances2 <- abundances2 + offset
  fc <- log2(abundances1/abundances2)
  
  # map fc to finite values
  fc[which(fc < -FcCutoff)] <- -FcCutoff
  fc[which(fc > +FcCutoff)] <- +FcCutoff
  
  return (fc)
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
  abline(v=0, col = "gray", lty=1)
  abline(v=median(fc.vector, na.rm=TRUE), col = "gray", lty=2)
  dev.off()
}

# plot distribution
plotDistribution <- function(vector, label, pdf.file) {
  pdf(file=pdf.file, height=4)
  density <- density(vector, na.rm=TRUE, bw="nrd0")
  plot(density$x, density$y, xlab=label, type="n", ylab="density", main="", yaxt='n')
  lines(density$x, density$y, col="gray", lwd=2)
  abline(v=0, col = "gray", lty=1)
  abline(v=median(vector, na.rm=TRUE), col = "gray", lty=2)
  dev.off()
}

# Kendrick nominal fractional mass plot
plotKendrick <- function(mass, pdf.file) {
  nominal = floor(mass)
  fractional = mass - nominal
  pdf(file=pdf.file)
  x <- nominal
  y <- fractional
  df <- data.frame(x,y)
  x <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  #cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
  cols <-  colorRampPalette(c("#2166AC", "#3F8EC0", "#80B9D8", "#BCDAEA", "#E6EFF3", "#F9EAE1", "#FAC8AF", "#ED9576", "#D25749", "#B2182B"))(256)
  df$col <- cols[df$dens]
  plot(y~x, data=df[order(df$dens),], pch=20, col=col, xlab='nominal mass [Da]', ylab='fractional mass [Da]')
  dev.off()
}

plotCorrelations <- function(data, pdf.file) {
    # extract study variables
    study_variables.index <- grepl("peptide_abundance_study_variable", colnames(peptide.data))
    study_variables.n <- sum(study_variables.index, na.rm=TRUE)
    study_variables.data = peptide.data[, study_variables.index]

    corr = cor(study_variables.data[complete.cases(study_variables.data),])

    # rename columns and rows
    colnames(corr) <- 1: study_variables.n
    rownames(corr) <- 1: study_variables.n
    cols <- colorRampPalette(c("#2166AC", "#3F8EC0", "#80B9D8", "#BCDAEA", "#E6EFF3", "#F9EAE1", "#FAC8AF", "#ED9576", "#D25749", "#B2182B"))(256)
    pdf(file=pdf.file)
    corrplot(corr, cl.lim=c(min(corr),max(corr)), col = cols, is.corr=FALSE)

    dev.off()

}

# Table with modifications
ModificationTable <- function(data) {
    modifications = data[!is.na(data)]
    df <- data.frame("Mod"=character(0), 
                     "City"=character(0), 
                     "Total"=double(0),
                     "Per Peptide"=character(0)
                    )
    df <- data.frame("allMods"=modifications)

    positions = NULL
    accessions = NULL
    peptide_ids = NULL

    id = 0
    for (modified_peptide in modifications) {
        id <- id + 1
        peptide_modifications <- strsplit(modified_peptide, ",")[[1]]
        for (modification in peptide_modifications) {
            position <- substr(modification, 1, 1)[[1]]
            accession <- (strsplit(modification, ":")[[1]][[2]])

            # XXX: Hardcode and recover mod + city combinations

            positions <- c(positions, position)
            accessions <- c(accessions, accession)
            
            peptide_ids <- c(peptide_ids, id)
        }
        # XXX: Per peptide?
    }
    df <- data.frame("ID"=peptide_ids, "Accession"=accessions, "Position"=positions)

    df2 <- unique(df[,c('Accession','Position')])
    df2["Total"] = 0

    # Increment "total counter" in unique dataframe for each row in non-unique dataframe.
    apply(df, 1, FUN=function(row) {
      row_accession = row["Accession"]
      row_position = row["Position"]
      unique_row <- df2[df2["Accession"] == row_accession & df2["Position"] == row_position,]
      df2[df2["Accession"] == row_accession & df2["Position"] == row_position,]["Total"] <<- unique_row["Total"] + 1
    })
    colnames(df2) <- c("Mod", "Site", "Total")

    ToModification <- function(accession) {
      if (accession == "1") {
        return ("Acetyl")
      } else if (accession == "4") {
        return ("Carbamidomethyl")
      } else if (accession == "737") {
        return ("TMT6plex")
      } else if (accession == "738") {
        return ("TMT2plex")
      }  else if (accession == "35") {
        return ("Oxidation")
      }
      else {
        stop(sprintf("Unsupported accession %s", accession))
      }
    }

    ToSite <- function(position) {
      return (position);
    }
    df2["Mod"] <- apply(FUN=ToModification, MARGIN=1, X=df2["Mod"])
    df2["Site"] <- apply(FUN=ToSite, MARGIN=1, X=df2["Site"])
    df2 <- df2[order(df2["Mod"], -df2["Total"]),]
    rownames(df2) <- c()
    return (df2)
}

# read mzTab data
peptide.data <- readMzTabPEP(input.file)

n.peptides = dim(peptide.data)[1]

# pre-define variables which will be called from LaTeX
# Needed even with IfFileExists()

median.abundance.1 <- 1
median.abundance.2 <- 1
median.abundance.3 <- 1

median.fc.12 <-0
median.fc.13 <-0
median.fc.23 <-0

sd.fc.12 <-0
sd.fc.13 <-0
sd.fc.23 <-0

# Kendrick plot
# plotKendrick((peptide.data$mass_to_charge - 1.00784) * peptide.data$charge, "plot_Kendrick.pdf")

# Modification Tabular
df <- ModificationTable(data=peptide.data$modifications)
write.table(df, "./test.csv", row.names=FALSE)


# plot peptide abundance distributions
# if (abundance.exists(peptide.data,1)) {
#   median.abundance.1 <- median(peptide.data$"peptide_abundance_study_variable[1]", na.rm=TRUE)
#   plotDistribution(log10(peptide.data$"peptide_abundance_study_variable[1]"), expression('log'[10]*' intensity'), "plot_DistributionIntensity_1.pdf")
# }
# if (abundance.exists(peptide.data,2)) {
#   median.abundance.2 <- median(peptide.data$"peptide_abundance_study_variable[2]", na.rm=TRUE)
#   plotDistribution(log10(peptide.data$"peptide_abundance_study_variable[2]"), expression('log'[10]*' intensity'), "plot_DistributionIntensity_2.pdf")
# }
# if (abundance.exists(peptide.data,3)) {
#   median.abundance.3 <- median(peptide.data$"peptide_abundance_study_variable[3]", na.rm=TRUE)
#   plotDistribution(log10(peptide.data$"peptide_abundance_study_variable[3]"), expression('log'[10]*' intensity'), "plot_DistributionIntensity_3.pdf")
# }
# 
# # plot fold change distributions and scatter plots
# if (abundance.exists(peptide.data,1) && abundance.exists(peptide.data,2)) {
#   a <- peptide.data$"peptide_abundance_study_variable[1]"
#   b <- peptide.data$"peptide_abundance_study_variable[2]"
#   fc <- calculateFoldChange(a, b)
#   intensity <- a + median(a/b, na.rm=TRUE)*b
#   median.fc.12 <- median(fc, na.rm=TRUE)
#   sd.fc.12 <- sd(fc, na.rm=TRUE)
#   plotFcLogIntensity(fc, intensity, "fold change", "plot_FoldChangeLogIntensity_12.pdf")
#   plotDistribution(fc, "fold change", "plot_DistributionFoldChange_12.pdf")
# }
# if (abundance.exists(peptide.data,1) && abundance.exists(peptide.data,3)) {
#   a <- peptide.data$"peptide_abundance_study_variable[1]"
#   b <- peptide.data$"peptide_abundance_study_variable[3]"
#   fc <- calculateFoldChange(a, b)
#   intensity <- a + median(a/b, na.rm=TRUE)*b
#   median.fc.13 <- median(fc, na.rm=TRUE)
#   sd.fc.13 <- sd(fc, na.rm=TRUE)
#   plotFcLogIntensity(fc, intensity, "fold change", "plot_FoldChangeLogIntensity_13.pdf")
#   plotDistribution(fc, "fold change", "plot_DistributionFoldChange_13.pdf")
# }
# if (abundance.exists(peptide.data,2) && abundance.exists(peptide.data,3)) {
#   a <- peptide.data$"peptide_abundance_study_variable[2]"
#   b <- peptide.data$"peptide_abundance_study_variable[3]"
#   fc <- calculateFoldChange(a, b)
#   intensity <- a + median(a/b, na.rm=TRUE)*b
#   median.fc.23 <- median(fc, na.rm=TRUE)
#   sd.fc.23 <- sd(fc, na.rm=TRUE)
#   plotFcLogIntensity(fc, intensity, "fold change", "plot_FoldChangeLogIntensity_23.pdf")
#   plotDistribution(fc, "fold change", "plot_DistributionFoldChange_23.pdf")
# }
# 
# # Determine number of study variables.
# study_variables.index <- grepl("peptide_abundance_study_variable", colnames(peptide.data))
# study_variables.n <- sum(study_variables.index, na.rm=TRUE)
# 
# # plot correlation matrix of peptide abundances
# if (study_variables.n >= 3) {
#     plotCorrelations(data = peptide.data, pdf.file = "correlations.pdf")
# }
