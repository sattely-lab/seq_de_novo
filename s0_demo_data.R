### this R code generate demo data (mc1:mc100) from a large data set
### Written by Xin Guan (github.com/x-guan)
library(tidyverse)
library(seqinr)

rm(list=ls())
setwd("/Users/Guan/Xin/R/Github/seq_de_novo")

dat_index_1 <- read.csv("/Users/Guan/Xin/Unix/SCG4_server/mca/180113_mca_getorf/mc_id.csv", header = T) %>%
  filter(contig_id %in% paste0("mc", 1:100)) %>%
  rename(contig_orf = contig) %>%
  rename(contig_id_orf = contig_id) %>%
  rename(length_orf = length)
tmp <- str_split(dat_index_1$contig_orf, "_", simplify=T)
dat_index_1 <- dat_index_1 %>%
  mutate(contig = paste(tmp[,1], tmp[,2], tmp[,3], tmp[,4], tmp[,5], sep="_")) %>%
  mutate(start_orf = tmp[,7]) %>%
  mutate(end_orf = tmp[,8])
dat_index_2 <- read.csv("/Users/Guan/Xin/Unix/SCG4_server/mca/180112_mca/mca_id.csv", header = T)
dat_index_3 <- dat_index_1 %>%
  left_join(dat_index_2, by="contig")
write.csv(dat_index_3, "data_demo/s0_mca_id.csv", row.names=F)

### blastp data
blast_re <- read.delim("/Users/Guan/Xin/Unix/SCG4_server/mca/result_mca/blast_mca/blastp_at2mca/at2mca.blastp", header = F)
contig_orf_demo <- read.csv("data_demo/s0_mca_id.csv", header=T) %>%
  select(contig_orf) %>% t %>% as.vector
blast_demo <- blast_re %>%
  filter(V2 %in% contig_orf_demo)
write.table(blast_demo, "data_demo/s0_at2mca.blastp", sep="\t", row.names=F, col.names=F)

### pfam data
pfam_re <- read.table("/Users/Guan/Xin/Unix/SCG4_server/mca/result_mca/hmmer_mca/mca_getorf_hmmer.domtblout", header = F, fill = T)
contig_orf_demo <- read.csv("data_demo/s0_mca_id.csv", header=T) %>%
  select(contig_orf) %>% t %>% as.vector
pfam_demo <- pfam_re %>%
  filter(V1 %in% contig_orf_demo)
write.table(pfam_demo, "data_demo/s0_mca.domtblout", row.names=F, col.names=F)

### sequence data
seq_re <- read.fasta("/Users/Guan/Xin/Unix/SCG4_server/mca/result_mca/bioawk_mca/mca_cd_600.fasta", seqtype="DNA")
contig_demo <- read.csv("data_demo/s0_mca_id.csv", header=T) %>%
  select(contig) %>% t %>% as.vector
seq_demo <- seq_re[getName(seq_re) %in% contig_demo]
write.fasta(sequences = seq_demo, names = getName(seq_demo), file.out = "data_demo/s0_mca_seq.fasta")

### count data
count_re <- read.csv("/Users/Guan/Xin/Unix/SCG4_server/mca/180113_mca_getorf/mc_count.csv", header=T)
count_re <- count_re[1:100,] %>%
  rename(contig_id_orf = contig_id, length_orf = length)
write.csv(count_re, "data_demo/s0_mca_count.csv", row.names=F)




