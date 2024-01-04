library(rtracklayer)
library(GenomicRanges)
library(Rsamtools)
library(dplyr)

gff <- import.gff("Mygtf.gtf")
gffgene<-getFeat2(x=gff,range_types=c("gene"),format="gtf")
gffgene <- getFeat2b(x=gff,range_types=c("gene"),format="gtf")


data <- read.csv(file="intergenic_dist.csv")
head(data)
summary(data)
data$X5_intergenic <- abs(data$X5_intergenic)
data$X3_intergenic <- abs(data$X3_intergenic)

data <- data %>% select(gene,strand,X5_intergenic,X3_intergenic)
colnames(data) <- c("geneid","strand","fiveprime","threeprime")

NumBins=40
if ((max(data$fiveprime,na.rm=TRUE)> max(data$threeprime,na.rm=TRUE)) == TRUE){
  FIR2Bin <- data$fiveprime
}else{FIR2Bin<- data$threeprime}

BinSteps <- round(length(FIR2Bin)/(NumBins-1),digits = 0)

FIR2BinOrd <- sort(FIR2Bin)

TempBinLimits <- FIR2BinOrd[seq(FIR2BinOrd[2*BinSteps],length(FIR2BinOrd),BinSteps)]

TempBinLimits[length(TempBinLimits)+1]<- max(FIR2Bin,na.rm = TRUE)

x <- seq(length(TempBinLimits))
fit <- nls(log(TempBinLimits) ~ a*x+b,start = c(a=0,b=0),algorithm  ='port',weights=((x-0.5*NumBins)^2))

pred=predict(fit,x)
BinLimits=c(1,round(exp(pred),0),max(FIR2Bin))
xbin=cut(data$fiveprime,breaks = c(BinLimits))
ybin=cut(data$threeprime,breaks=c(BinLimits))

FIRdata <- cbind(data,xbin,ybin,genevalue=rep(1,length(data$fiveprime)))

GeneValMatrix <- with(FIRdata,tapply(genevalue,list(xbin,ybin),sum))

x <- 1:ncol(GeneValMatrix)
y <- 1:nrow(GeneValMatrix)
zlim=range(as.numeric(unlist(GeneValMatrix)),finite=TRUE)

mypalette <- colorRampPalette(c("white","darkblue","forestgreen","goldenrod1","orangered","red3","darkred"),space="rgb")

mycol=mypalette(2*max(GeneValMatrix,na.rm = TRUE))

mylabels<- paste(BinLimits[1:length(BinLimits)-1],BinLimits[2:length(BinLimits)],sep="-",collapse = NULL)
mylabels <- head(mylabels,-1)
filled.contour(x,y,z=GeneValMatrix,
               plot.title = title(main="Phytophthora meadii genome",xlab="five prime intergenic regions",ylab="three prime intergenic regions",
                                  cex.main=0.8,cex.lab=0.5),key.title = title(main="Number of genes in bins",cex.main=0.5,line=1),
               col=mycol,levels = pretty(zlim,2*max(GeneValMatrix,na.rm = TRUE)),plot.axes={axis(1,at=x,labels=mylabels,las=2,cex.axis=0.5);
                 axis(2,at=y,labels=mylabels,cex.axis=0.5)})

