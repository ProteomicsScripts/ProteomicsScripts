
# load analysis
analysis <- read.delim("analysis.csv", header=TRUE, sep=",", fill=TRUE, skip=5, stringsAsFactors=FALSE)
analysis <- data.frame(analysis[which(analysis[,1]=="CONSENSUS"),])
for(i in 2:dim(analysis)[2]) {
	analysis[,i] <- as.numeric(analysis[,i])
}

## sort by intensity
analysis <- analysis[order(analysis$intensity_cf, decreasing=TRUE),] 


pdf(file="RatioIntensity.pdf")
medianRatio <- median(analysis$intensity_1/analysis$intensity_0)
x <- log2(analysis$intensity_1/analysis$intensity_0)
y <- log10(analysis$intensity_1 + medianRatio * analysis$intensity_0)
df <- data.frame(x,y)
x <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", "#FCFF00", "#FF9400", "#FF3100"))(256)
df$col <- cols[df$dens]
plot(y~x, data=df[order(df$dens),], pch=20, col=col, xlab="fc (H:L)", ylab=expression('log'[10]*' intensity'))
abline(v=log2(medianRatio), col = "gray")
dev.off()

pdf(file="densityRatios.pdf", height=4)
density <- density(log2(analysis$intensity_1/analysis$intensity_0), na.rm=TRUE, bw="nrd0")
plot(density$x, density$y, xlab="fc (H:L)", type="n", ylab="density", main="", yaxt='n')
lines(density$x, density$y, col="#000099", lwd=2)
abline(v=log2(medianRatio), col = "gray")
dev.off()

pdf(file="densityIntensities.pdf", height=4)
density <- density(log10(analysis$intensity_1 + medianRatio * analysis$intensity_0), na.rm=TRUE, bw="nrd0")
plot(density$x, density$y, xlab=expression('log'[10]*' intensity'), type="n", ylab="density", main="", yaxt='n')
lines(density$x, density$y, col="#000099", lwd=2)
dev.off()

numberPairs <- dim(analysis)[1]
