### This R code calculates Open Reading Frame (ORF) of each assembled contig 
### Written by Xin Guan (github.com/x-guan)
library(tidyverse)
library(seqinr)

rm(list=ls())
setwd("/Users/Guan/Xin/R/Github/seq_de_novo")

GetPrimer <- function(InputSeq, OutputName, InsertType){
  annot_1 <- getAnnot(InputSeq)
  annot_2 <- t(as.data.frame(str_split(annot_1, " | ")))
  row.names(annot_2) <- c()
  annot_2[,1] <- str_replace(annot_2[,1], ">", "")
  input_annot <- annot_2[,c(1,7,9)]
  colnames(input_annot) <- c("contig_id_2", "start_codon", "stop_codon")
  
  # design F primers
  primer_f <- NULL
  for (i in 1:length(InputSeq)){
    fasta_tmp <- InputSeq[[i]]
    for (j in 1:9){
      primer <- getFrag(fasta_tmp, 1, 22 + (j - 1))
      gc_cont <- GC(primer)
      if (gc_cont > 0.45 - (j - 1) / 100 * 2) break
    }
    primer_info <- cbind.data.frame(seq_f = str_to_upper(paste(primer, collapse = "")), 
                                    len_f = 22 + (j - 1), gc_cont_f = gc_cont)
    primer_f <- rbind.data.frame(primer_f, primer_info)
  }
  
  # design R primers 
  primer_r <- NULL
  for (i in 1:length(InputSeq)){
    fasta_tmp <- InputSeq[[i]]
    fasta_tmp_cr <- comp(rev(fasta_tmp))
    for (j in 1:9){
      primer <- getFrag(fasta_tmp_cr, 1, 22 + (j - 1))
      gc_cont <- GC(primer)
      if (gc_cont > 0.45 - (j - 1) / 100 * 2) break
    }
    primer_info <- cbind.data.frame(seq_r = str_to_upper(paste(primer, collapse = "")), 
                                    len_r = 22 + (j - 1), gc_cont_r = gc_cont)
    primer_r <- rbind.data.frame(primer_r, primer_info)
  }
  
  input_annot <- cbind.data.frame(input_annot, primer_f, primer_r)
  
  # generate primer seq to be synthesized
  if (InsertType == "Gibson_pEAQ_AgeI/XhoI"){
    linker_f <- "tattctgcccaaattcgcgaccggt"   # 5' linker used to insert a gene into pEAQ-HT
    linker_r <- "tgaaaccagagttaaaggcctcgag"   # 3' linker used to insert a gene into pEAQ-HT
    input_annot$primer_f <- with(input_annot, ifelse(start_codon == "y", 
                                                     paste0(linker_f, seq_f),
                                                     paste0(linker_f, "atg", seq_f)))
    input_annot$primer_r <- with(input_annot, ifelse(stop_codon == "y", 
                                                     paste0(linker_r, seq_r),
                                                     paste0(linker_r, "tta", seq_r)))
  }
    
  # output
  write.csv(input_annot, OutputName, row.names = F)
}

##########################################################################################################################
### design primers
# input
InputSeq <- read.fasta("data_demo/s3_mca_nucl.fasta", seqtype = "DNA")   # input nucl seq
OutputName <- "data_demo/s4_mca_primer.csv"   # define output primer file
InsertType <- "Gibson_pEAQ_AgeI/XhoI"
# execute
GetPrimer(InputSeq, OutputName, InsertType)



