### This R code calculates
# 1) read number-normalized counts
# 2) fpkm
# 3) TMM-normalized counts 
### Written by Xin Guan (github.com/x-guan)
library(tidyverse)
library(edgeR)

rm(list=ls())
setwd("/Users/Guan/Xin/R/Github/seq_de_novo")

# input count info
count_info <- read.csv("data_demo/s0_mca_count.csv", header=T)
# input total read number of each sample
read_num <- read.csv("data_demo/s0_mca_read_num.csv", header = T)
read_num_pair <- read_num[seq(1,47,2),]
co_eff <- with(read_num_pair, read_num_2/(sum(read_num_2)/24))

### calculate read number-normalized counts
count_norm <- cbind.data.frame(count_info[,1:3], 
                                 round(sapply(1:24, function(i){count_info[,4:27][,i]/co_eff[i]})))
colnames(count_norm) <- colnames(count_info)
write.csv(count_norm, "data_demo/s5_count_norm.csv", row.names=F)

### calculate fpkm
fpk <- count_info[,4:27]/count_info$length*1000
fpkm <- cbind.data.frame(count_info[,1:3], 
                         sapply(1:24, function(i){fpk[, i]/read_num_pair[,2][i]*1000000}))
colnames(fpkm) <- colnames(count_info)
write.csv(fpkm, "data_demo/s5_fpkm.csv", row.names=F)

### calculate tmm-normalized counts
tmm_coeff <- count_info %>%
  column_to_rownames("contig_id_orf") %>%
  select(-c(contig, length_orf)) %>%
  calcNormFactors
count_tmm <- cbind.data.frame(count_info[,1:3],
                              round(sapply(1:24, function(i){count_info[,4:27][, i]/tmm_coeff[i]})))
colnames(count_tmm) <- colnames(count_info)
write.csv(count_tmm, file = "data_demo/s5_count_tmm.csv", row.names = F)


