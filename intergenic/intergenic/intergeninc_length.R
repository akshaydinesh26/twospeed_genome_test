library(dplyr)
library(readxl)

data <- read_excel("Akshay/4_ongoing/fungal_assembly/intergenic/intergenic.xlsx")

head(data)
dim(data)

int <- data %>% group_by(scaffold) %>% mutate(inter=start-lag(stop,default=first(stop))) %>% dplyr::mutate(rank = 1:length(scaffold)) %>% filter(rank < max(rank)) %>%
 filter(rank> min(rank)) %>% dplyr::mutate(rank = NULL)
View(int)

head(int)

data2 <- int %>% select(scaffold,strand,inter,mRNA,gene)
head(data2)

three <- c()
five <- c()

for(a in 1:nrow(data2)){
  if(data2[a,"strand"]=="+"){
    five <- c(five,data2[a,"inter"])
    three <- c(three,data2[a+1,"inter"])
  }
  else if(data2[a,"strand"]=="-") {
    three <- c(three,data2[a,"inter"])
    five <- c(five,data2[a+1,"inter"])
  }
}

data3 <- data.frame(data2,unlist(five),unlist(three))
head(data3)

final <- data3 %>% select(scaffold,`unlist.five.`,`unlist.three.`)
final <- data.frame(final,data2$mRNA,data2$gene)
colnames(final) <- c("scaffold","five_intergenic","three_intergenic","mRNA","gene")
head(final)
write.csv(final,file="intergenic_dist.csv",row.names = F)
final$`5_intergenic` <- abs(final$`5_intergenic`)
final$`3_intergenic` <- abs(final$`3_intergenic`)
summary(final)

head(final)

library('epiDisplay')

a <- tab1(final$five_intergenic,graph=F)
a <- table(cut(final$five_intergenic,seq))

bins <- quantile(final$three_intergenic,prob=seq(0,1,0.05),na.rm=T)

a <- table(cut(final$five_intergenic,bins))
b <- table(cut(final$three_intergenic,bins))
data_out <- cbind(a,b)
data_out
write.csv(data_out,file="count_bins.csv")

library(ggplot2)
plot(final$five_intergenic,final$three_intergenic)

cont_table <- table(cut(final$five_intergenic,bins),cut(final$three_intergenic,bins))
xt <- as.vector(bins)
xt
library("plotly")
fig <- plot_ly(x=xt,y=xt,z=cont_table,type="contour")
png(file="a.png")
fig
dev.off()
bins
cont_table
