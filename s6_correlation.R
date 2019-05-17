### This R code calculates Pearson correlation of the gene expression data
### Written by Xin Guan (github.com/x-guan)
library(tidyverse)

rm(list=ls())
setwd("/Users/Guan/Xin/R/Github/seq_de_novo")

### calculate the correlation co-efficient and ranking between each gene of interest (e.g. mc1, mc2 and mc3) and all contigs
mca_fpkm <- read.csv("data_demo/s5_fpkm.csv", header=T)
row.names(mca_fpkm) <- mca_fpkm$contig_id
mca_fpkm_entire <- mca_fpkm[,4:27] %>% t
mca_fpkm_partial <- mca_fpkm[mca_fpkm$contig_id_orf %in% c("mc1", "mc2", "mc3"),4:27] %>% t
mca_cor <- cor(mca_fpkm_partial, mca_fpkm_entire) %>% t %>% round(digits=2)
mca_cor <- data.frame(contig_id = row.names(mca_cor), mca_cor)
mca_cor_2 <- mca_cor %>%
  gather(contig_id_cand, r, -contig_id)

# extract r2 and ranking info for each gene of interest
contig_cand <- unique(mca_cor_2$contig_id_cand) %>% as.vector
for (i in 1:length(contig_cand)){
  cor_tmp <- mca_cor_2 %>%
    filter(contig_id_cand == contig_cand[i]) %>%
    arrange(-r) %>%
    select(-contig_id_cand)
  cor_tmp$ranking <- 1:nrow(cor_tmp)
  write.csv(cor_tmp, paste0("data_demo/s6_cor_", contig_cand[i], ".csv"), row.names=F)
}
