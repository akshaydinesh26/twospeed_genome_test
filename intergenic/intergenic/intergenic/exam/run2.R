# load filled contour 3
install.packages("png")
library(png)
library(gridExtra)
library(ggplot2)

png(filename = paste("afs",".png",sep=""))
filled.contour3(x,y,z=GeneValMatrix,col=mycol,levels=pretty(zlim,2*max(GeneValMatrix,na.rm=TRUE)),
                frame.plot = FALSE,axes=FALSE)
dev.off()
getwd()
