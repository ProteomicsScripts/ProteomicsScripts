df <- data.frame("DT51", "null", 1, 1, "null")
colnames(df) <- c('Raw.file', 'Condition', 'BioReplicate', 'Run', 'IsotopeLabelType')
write.csv(df, file="./annotation.csv")
