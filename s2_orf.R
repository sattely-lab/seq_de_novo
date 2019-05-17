### This R code transfer annotation (BLAST and Pfam) into numbers, including
# start/end of homologous sections (homo_start and homo_end), 
# start/end of open reading frame (orf_start and orf_end),
# start/end of coding sequence (cds_start and cds_end), and
# whether or not each orf contains start/stop codons (y - present, n1 - absent due to orf edge, n2 absent due to a upstream stop codon)
### Written by Xin Guan (github.com/x-guan)
library(tidyverse)
library(seqinr)

rm(list=ls())
setwd("/Users/Guan/Xin/R/Github/seq_de_novo")

##### A function to select orf according to blast/pfam results
CalORF <- function(dat_annot, dat_fasta, file_name){
  # calculate query_start and query_end
  dat_annot$query_start <- with(dat_annot, ifelse(getorf_start < getorf_end, 
                                                  getorf_start + (contig_start - 1)*3,
                                                  getorf_start - (contig_start - 1)*3))
  dat_annot$query_end <- with(dat_annot, ifelse(getorf_start < getorf_end, 
                                                getorf_start - 1 + contig_end*3,
                                                getorf_start + 1 - contig_end*3))
  # for loop to search against each component in the dat_blast1 fasta file
  orf_blast <- NULL
  for (i in 1:nrow(dat_annot)){
    blast_tmp <- dat_annot[i,]
    # collect info about "contig", "length", "contig_orf", "contig_id_orf", "length_orf",  "query_start", "query_end"
    # "desc_1", "desc_2", "desc_3"
    contig <- as.character(blast_tmp$contig)
    length <- as.integer(blast_tmp$length)
    contig_orf <- as.character(blast_tmp$contig_orf)
    contig_id_orf <- as.character(blast_tmp$contig_id_orf)
    length_orf <- as.character(blast_tmp$length_orf)
    homo_start <- blast_tmp$query_start
    homo_end <- blast_tmp$query_end
    desc_1 <- blast_tmp$desc_1
    desc_2 <- blast_tmp$desc_2
    desc_3 <- blast_tmp$desc_3
    
    # if/else to work on forward and reverse direction respectively
    if (blast_tmp$query_start < blast_tmp$query_end) {
      #########################################################
      ### contig was translated in the forward direction
      ### calculate orf
      left_del <- (homo_start - 1) %% 3   # num of base to remove at the 5' end
      right_del <- (length - left_del) %% 3   # num of base to remove at the 3' end
      orf_start <- left_del + 1   # the 1st base of orf
      orf_end <- length - right_del   # the last base of orf
      ### calculate cds
      fasta1 <- dat_fasta[attributes(dat_fasta)$name %in% blast_tmp$contig]
      fasta1 <- fasta1[[1]]
      # search for all stop codons (position of the 1st base) upstream of the homologous sequence
      stop_codon_us <- NULL
      for (j in orf_start:homo_start){
        if (((j - orf_start) %% 3 == 0 & fasta1[j] == "t" & fasta1[j + 1] == "a" & fasta1[j + 2] == "a") |
            ((j - orf_start) %% 3 == 0 & fasta1[j] == "t" & fasta1[j + 1] == "a" & fasta1[j + 2] == "g") |
            ((j - orf_start) %% 3 == 0 & fasta1[j] == "t" & fasta1[j + 1] == "g" & fasta1[j + 2] == "a"))
          stop_codon_us <- c(stop_codon_us, j)
      }
      # search for all stop codons (position of the 3rd base) downstream of the homologous sequence
      stop_codon_ds <- NULL
      for (k in homo_end:orf_end){
        if (((orf_end - k) %% 3 == 0 & fasta1[k - 2] == "t" & fasta1[k - 1] == "a" & fasta1[k] == "a") |
            ((orf_end - k) %% 3 == 0 & fasta1[k - 2] == "t" & fasta1[k - 1] == "a" & fasta1[k] == "g") |
            ((orf_end - k) %% 3 == 0 & fasta1[k - 2] == "t" & fasta1[k - 1] == "g" & fasta1[k] == "a"))
          stop_codon_ds <- c(stop_codon_ds, k)
      }
      # search for all start codons (position of the 1st base) upstream of the homologous sequence
      start_codon_us <- NULL
      for (p in orf_start:homo_start){
        if ((p - orf_start) %% 3 == 0 & fasta1[p] == "a" & fasta1[p + 1] == "t" & fasta1[p + 2] == "g")
          start_codon_us <- c(start_codon_us, p)
      }
      # search for the 1st start codon behind the upstream stop codon
      if (length(stop_codon_us)){
        if (length(start_codon_us[start_codon_us > max(stop_codon_us)])){
          cds_start <- min(start_codon_us[start_codon_us > max(stop_codon_us)])
          det_start <- "y"   # cds starts with a start codon
        } else{
          cds_start <- max(stop_codon_us)
          det_start <- "stop"   # start codon is missing between upstream stop codon and homologous sequence
        }
      } else {
        if (length(start_codon_us)){
          cds_start <- min(start_codon_us)
          det_start <- "y"   # cds starts with a start codon
        } else{
          cds_start <- orf_start
          det_start <- "orf"   # start codon is missing due to the 5' edge of orf
        }
      }
      # search for cds_end (the 1st start codon after the homologous portion)
      if (length(stop_codon_ds)){
        cds_end <- min(stop_codon_ds)
        det_stop <- "y"   # cds ends with a stop codon
      } else {
        cds_end <- orf_end
        det_stop <- "orf"   # stop codon is missing due to the 3' edge of orf
      }
      #########################################################
    } else {
      #########################################################
      ### contig was translated in the reverse direction
      ### calculate orf
      left_del <- (length - homo_start) %% 3   # num of base to remove at the 5' end
      right_del <- (homo_end - 1) %% 3   # num of base to remove at the 3' end
      orf_start <- length - left_del   # the 1st base of orf
      orf_end <- right_del + 1   # the last base of orf
      ### calculate cds
      fasta1 <- dat_fasta[attributes(dat_fasta)$name %in% blast_tmp$contig]
      fasta1 <- fasta1[[1]]
      # search for all stop codons (position of the 1st base) upstream of the homologous sequence
      stop_codon_us <- NULL
      for (j in homo_start:orf_start){
        if (((orf_start - j) %% 3 == 0 & fasta1[j] == "a" & fasta1[j - 1] == "t" & fasta1[j - 2] == "t") |
            ((orf_start - j) %% 3 == 0 & fasta1[j] == "a" & fasta1[j - 1] == "t" & fasta1[j - 2] == "c") |
            ((orf_start - j) %% 3 == 0 & fasta1[j] == "a" & fasta1[j - 1] == "c" & fasta1[j - 2] == "t"))
          stop_codon_us <- c(stop_codon_us, j)
      }
      # search for all stop codons (position of the 3rd base) downstream of the homologous sequence
      stop_codon_ds <- NULL
      for (k in orf_end:homo_end){
        if (((k - orf_end) %% 3 == 0 & fasta1[k + 2] == "a" & fasta1[k + 1] == "t" & fasta1[k] == "t") |
            ((k - orf_end) %% 3 == 0 & fasta1[k + 2] == "a" & fasta1[k + 1] == "t" & fasta1[k] == "c") |
            ((k - orf_end) %% 3 == 0 & fasta1[k + 2] == "a" & fasta1[k + 1] == "c" & fasta1[k] == "t"))
          stop_codon_ds <- c(stop_codon_ds, k)
      }
      # search for the 1st start codon behind the upstream stop codon
      start_codon_us <- NULL
      for (p in homo_start:orf_start){
        if ((orf_start - p) %% 3 == 0 & fasta1[p] == "t" & fasta1[p - 1] == "a" & fasta1[p - 2] == "c")
          start_codon_us <- c(start_codon_us, p)
      }
      # search for cds_start (the 1st start codon after the last stop codon)
      if (length(stop_codon_us)){
        if (length(start_codon_us[start_codon_us < min(stop_codon_us)])){
          cds_start <- max(start_codon_us[start_codon_us < min(stop_codon_us)])
          det_start <- "y"   # cds starts with a start codon
        } else{
          cds_start <- min(stop_codon_us)
          det_start <- "stop"   # start codon is missing between upstream stop codon and homologous sequence
        }
      } else {
        if (length(start_codon_us)){
          cds_start <- max(start_codon_us)
          det_start <- "y"   # cds starts with a start codon
        } else{
          cds_start <- orf_start
          det_start <- "orf"   # start codon is missing due to the 5' edge of orf
        }
      }
      # search for cds_end (the 1st start codon after the homologous portion)
      if (length(stop_codon_ds)){
        cds_end <- max(stop_codon_ds)
        det_stop <- "y"   # cds ends with a stop codon
      } else {
        cds_end <- orf_end
        det_stop <- "orf"   # stop codon is missing due to the 3' edge of orf
      }
      #########################################################
    }
    ### determine if a contig contains full length cds (only if both start codon and stop codon occur)
    if (det_start == "y" & det_stop =="y"){
      det_full_cds <- "y"   # the contig contains full length cds
    } else {
      if (det_start == "stop") {
        det_full_cds <- "n2"   # the contig must be a product of mis-assembly
      } else
        det_full_cds <- "n1"   # the contig could be a gene fragment
    } 
    # merge data
    orf_blast_tmp <- data.frame(contig = contig, length = length, contig_orf = contig_orf, contig_id_orf = contig_id_orf, length_orf = length_orf,
                                desc_1 = desc_1, desc_2 = desc_2,  desc_3 = desc_3,
                                homo_start = homo_start, homo_end = homo_end, 
                                orf_start = orf_start, orf_end = orf_end, 
                                cds_start = cds_start, cds_end = cds_end, 
                                start_codon = det_start, stop_codon = det_stop, full_cds = det_full_cds)
    orf_blast <- rbind.data.frame(orf_blast, orf_blast_tmp, make.row.names = F)
    # check progress
    if(i%%100 == 0) print(i)
  }
  write.csv(orf_blast, file = file_name, row.names = F)
}

##########################################################################################################################
### Calculate ORF using BLAST data
# columns required in dat_annot:
# "contig", "contig_id",  "length", "getorf_start", "getorf_end",
# "contig_start", "contig_end", "desc_1", "desc_2", "desc_3"

# input blastp results
blast_hit1 <- read.csv("data_demo/s1_at2mca_hit1_subject.csv", header = T)
colnames(blast_hit1)[c(1,2,7:10)] <- c("at_homolog", "contig_orf", "at_start", "at_end", "contig_start", "contig_end")
# input other required info
mca_id <- read.csv("data_demo/s0_mca_id.csv", header=T)
dat_annot <- mca_id %>%
  left_join(blast_hit1, by="contig_orf") %>%
  select(contig, length, contig_orf, contig_id_orf, length_orf, getorf_start = start_orf, getorf_end = end_orf, 
         at_homolog, contig_start, contig_end, identity_perc, alignment_len, expect, score)
# add 3 descriptions: 1) Arabidopsis homolog, 2) % of identity, and 3) evalue
dat_annot$desc_1 <- dat_annot$at_homolog
dat_annot$desc_2 <- with(dat_annot, round((alignment_len*identity_perc/100)/(length_orf/3)*100))
dat_annot$desc_3 <- dat_annot$expect
dat_annot <- dat_annot[!is.na(dat_annot$at_homolog),]
# input 2/3: fasta info (original nucl sequences)
dat_fasta <- read.fasta("data_demo/s0_mca_seq.fasta")
# input 3/3: name of the output csv file
file_name <- "data_demo/s2_orf_blast.csv"
# execute
CalORF(dat_annot, dat_fasta, file_name)

##########################################################################################################################
### Calculate ORF using Pfam dat
pfam_hit1 <- read.csv("data_demo/s1_mca_pfam_hit1.csv", header = T)
colnames(pfam_hit1)[c(1)] <- c("contig_orf")
# input other required info
mca_id <- read.csv("data_demo/s0_mca_id.csv", header=T)
dat_annot <- mca_id %>%
  left_join(pfam_hit1, by="contig_orf") %>%
  select(contig, length, contig_orf, contig_id_orf, length_orf,
         getorf_start = start_orf, getorf_end = end_orf, contig_start, contig_end, 
         pfam_profile, pfam_accession, expect = i_evalue_dom, score = score_dom)
# add 3 descriptions: 1) Arabidopsis homolog, 2) % of identity, and 3) evalue
dat_annot$desc_1 <- dat_annot$pfam_profile
dat_annot$desc_2 <- dat_annot$pfam_accession
dat_annot$desc_3 <- dat_annot$expect
dat_annot <- dat_annot[!is.na(dat_annot$pfam_profile),]
# input 2/3: fasta info (original nucl sequences)
dat_fasta <- read.fasta("data_demo/s0_mca_seq.fasta")
# input 3/3: name of the output csv file
file_name <- "data_demo/s2_orf_pfam.csv"
# execute
CalORF(dat_annot, dat_fasta, file_name)

##########################################################################################################################
### merge blast hit and pfam hit
dat_blast <- read.csv("data_demo/s2_orf_blast.csv", header = T)
colnames(dat_blast)[6:17] <- paste(colnames(dat_blast)[6:17], "bt", sep = "_")
dat_pfam <- read.csv("data_demo/s2_orf_pfam.csv", header = T)
colnames(dat_pfam)[6:17] <- paste(colnames(dat_pfam)[6:17], "pm", sep = "_")
dat_bt_pm <- merge(dat_blast, dat_pfam, by = c("contig", "length", "contig_orf", "contig_id_orf", "length_orf"), all.x = T, all.y = T)
write.csv(dat_bt_pm, "data_demo/s2_orf_merge.csv", row.names = F)


