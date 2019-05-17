### This R code gathers the ORF info calculated from CalORF function,
# determines if BLAST or Pfam annotations suggested different CDS, and
# generates prot/nucl sequence info 
### Written by Xin Guan (github.com/x-guan)
library(tidyverse)
library(seqinr)

rm(list=ls())
setwd("/Users/Guan/Xin/R/Github/seq_de_novo")

GetORF <- function(InputORF, InputFASTA, OutputNucl, OutputProt){
  orf_mca <- read.csv(InputORF, header = T)
  orf_blast_pfam <- orf_mca
  cds_info <- NULL
  for (i in 1:nrow(orf_blast_pfam)){
    orf_tmp <- orf_blast_pfam[i,]
    if (is.na(orf_tmp$cds_start_bt)){
      cds_tmp <- cbind.data.frame(contig = orf_tmp$contig, contig_orf = orf_tmp$contig_orf, contig_id_orf = orf_tmp$contig_id_orf,
                                  indicator = "p", start = orf_tmp$cds_start_pm, end = orf_tmp$cds_end_pm,
                                  start_codon = orf_tmp$start_codon_pm, stop_codon = orf_tmp$stop_codon_pm,
                                  full_cds = orf_tmp$full_cds_pm)
    } else {
      if (is.na(orf_tmp$cds_start_pm)) {
        cds_tmp <- cbind.data.frame(contig = orf_tmp$contig, contig_orf = orf_tmp$contig_orf, contig_id_orf = orf_tmp$contig_id_orf,
                                    indicator = "b", start = orf_tmp$cds_start_bt, end = orf_tmp$cds_end_bt, 
                                    start_codon = orf_tmp$start_codon_bt, stop_codon = orf_tmp$stop_codon_bt,
                                    full_cds = orf_tmp$full_cds_bt)
        
      } else {
        if (orf_tmp$cds_start_bt == orf_tmp$cds_start_pm & orf_tmp$cds_end_bt == orf_tmp$cds_end_pm){
          cds_tmp <- cbind.data.frame(contig = orf_tmp$contig, contig_orf = orf_tmp$contig_orf, contig_id_orf = orf_tmp$contig_id_orf,
                                      indicator = "bp",  start = orf_tmp$cds_start_bt, end = orf_tmp$cds_end_bt, 
                                      start_codon = orf_tmp$start_codon_bt, stop_codon = orf_tmp$stop_codon_bt,
                                      full_cds = orf_tmp$full_cds_bt)
        } else {
          cds_tmp_bt <- cbind.data.frame(contig = orf_tmp$contig, contig_orf = orf_tmp$contig_orf, contig_id_orf = orf_tmp$contig_id_orf, 
                                         indicator = "bx", start = orf_tmp$cds_start_bt, end = orf_tmp$cds_end_bt, 
                                         start_codon = orf_tmp$start_codon_bt, stop_codon = orf_tmp$stop_codon_bt,
                                         full_cds = orf_tmp$full_cds_bt)
          cds_tmp_pm <- cbind.data.frame(contig = orf_tmp$contig, contig_orf = orf_tmp$contig_orf, contig_id_orf = orf_tmp$contig_id_orf, 
                                         indicator = "xp", start = orf_tmp$cds_start_pm, end = orf_tmp$cds_end_pm, 
                                         start_codon = orf_tmp$start_codon_pm, stop_codon = orf_tmp$stop_codon_pm,
                                         full_cds = orf_tmp$full_cds_pm)
          cds_tmp <- rbind.data.frame(cds_tmp_bt, cds_tmp_pm)
        }
      }
    } 
    cds_info <- rbind.data.frame(cds_info, cds_tmp)
  }
  cds_info$contig_id_2 <- with(cds_info, paste0(contig_id_orf, ifelse(indicator == "bx", "bx",
                                                                      ifelse(indicator == "xp", "xp", ""))))
  cds_info$aa_num <- with(cds_info, (abs(end - start) + 1)%/%3)
  cds_info$msa_name <- with(cds_info, paste(contig_id_2, indicator, aa_num, full_cds, sep = "|"))
  cds_info$fasta_name <- with(cds_info, paste(contig_id_2, indicator, aa_num, start_codon, stop_codon, full_cds, sep = " | "))
  write.csv(cds_info, OutputCDS, row.names = F)
  
  ###########################################################
  ### extract nucl/prot sequences
  fastadt_mca <- read.fasta(InputFASTA)
  fasta_nucl <- list()
  fasta_prot <- list()
  for (i in 1:nrow(cds_info)){
    fasta_tmp <- fastadt_mca[getName(fastadt_mca) == cds_info$contig[i]]   ### stopped here, need full-length contig
    if (cds_info$start[i] < cds_info$end[i]){
      fasta_tmp_cds <- getFrag(fasta_tmp, cds_info$start[i], cds_info$end[i])
      fasta_tmp_cds <- fasta_tmp_cds[[1]]
      fasta_tmp_prot <- translate(fasta_tmp_cds)
    } else {
      fasta_tmp_cds <- getFrag(fasta_tmp, cds_info$start[i], cds_info$end[i])
      fasta_tmp_cds <- comp(fasta_tmp_cds[[1]])
      fasta_tmp_prot <- translate(fasta_tmp_cds)
    }
    fasta_nucl <- c(fasta_nucl, list(fasta_tmp_cds))
    fasta_prot <- c(fasta_prot, list(fasta_tmp_prot))
  }
  write.fasta(sequences = fasta_nucl, names = cds_info$fasta_name, file.out = OutputNucl)
  write.fasta(sequences = fasta_prot, names = cds_info$fasta_name, file.out = OutputProt)
}

##########################################################################################################################
### get cds info, nucl seq and prot seq
# input: path to contig_getorf id info (do not change)
InputORF <- "data_demo/s2_orf_merge.csv"
InputFASTA <- "data_demo/s0_mca_seq.fasta"   # input nucl seq
OutputCDS <- "data_demo/s3_mca_cds.csv"   # define output file containing cds info
OutputNucl <- "data_demo/s3_mca_nucl.fasta"   # define output nucl seq file 
OutputProt <- "data_demo/s3_mca_prot.fasta"   # define output prot seq file
# execute
GetORF(InputORF, InputFASTA, OutputNucl, OutputProt)

