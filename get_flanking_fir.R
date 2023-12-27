#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error.
if (length(args) == 0 || length(args) == 1) {
 stop("At leat two arguments must be supplied", call. = FALSE)
}

library(dplyr)

data <- read.csv(file = args[1],header = FALSE,col.names = c("chr","start","end","strand","gene_name"))

# distance to the previous gene for all genes regardless of strand so that
# a gene have inter value which is distance b/n its start and end of previous gene
# rank was used to remove gene in the edges of the chromosome

flanking_data <- data %>% group_by(chr) %>% arrange(start,end,.by_group = TRUE) %>% mutate(inter=start-lag(end,default=first(end))) %>%
                 mutate(rank = 1:length(chr)) %>% filter(rank < max(rank)) %>%
                 filter(rank > min(rank)) %>% mutate(rank = NULL)

# negative distance is overlap so that is changed to zero length of flanking distance.
flanking_data[which(flanking_data$inter < 0),"inter"] <- 0

# assign true 5' and 3' FIR for all genes
# for strand "-" inter of next gene is 5' FIR and its inter is 3' FIR.
# for strand "+" its inter is 5' FIR and inter of next gene is 3' FIR.
final_output <- flanking_data %>% group_by(chr) %>% arrange(start,end,.by_group = TRUE) %>%
                mutate(five_p = if_else(strand == "+",inter,lead(inter))) %>%
                mutate(three_p = if_else(strand == "+",lead(inter),inter)) %>%
                ungroup() %>% select(gene_name,five_p,three_p)

#rename_columns
colnames(final_output) <- c("gene","five_p_fir","three_p_fir")

write.csv(final_output,file = args[2],quote = FALSE,row.names = FALSE)
