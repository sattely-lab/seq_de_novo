### This R code calculates the best hit for each contig in the blast or pfam datasets
### Written by Xin Guan (github.com/x-guan)

rm(list=ls())
setwd("/Users/Guan/Xin/R/Github/seq_de_novo")

### GetBlastHit1: a function that searches for the best blastp hit for each subject or query
GetBlastHit1 <- function(blast_re, output_csv_file, function_mode){
  colnames(blast_re) <- c("query", "subject", "identity_perc", "alignment_len",
                          "mismatch_num", "gap_num", "query_start", "query_end",
                          "subject_start", "subject_end", "expect", "score")
  if (function_mode == "subject"){
    index_blast <- NULL
    for (i in 1:length(unique(blast_re$subject))){
      index_tmp1 <- with(blast_re, which(subject == unique(subject)[i]))
      index_tmp2 <- with(blast_re, min(index_tmp1[score[index_tmp1] == max(score[index_tmp1])]))
      index_blast <- c(index_blast, index_tmp2)
      if(i%%100 == 0) print(i)
    }
  } else {
    if (function_mode == "query"){
      index_blast <- NULL
      for (i in 1:length(unique(blast_re$query))){
        index_tmp1 <- with(blast_re, which(query == unique(query)[i]))
        index_tmp2 <- with(blast_re, min(index_tmp1[score[index_tmp1] == max(score[index_tmp1])]))
        index_blast <- c(index_blast, index_tmp2)
        if(i%%100 == 0) print(i)
      }
    } else {
      print ("Function_mode was not recognized!")
    }
  }
  blast_re_hit1 <- blast_re[index_blast,]
  # save
  write.csv(blast_re_hit1, file = output_csv_file, row.names = F)
}

# # GetBlastHit1s: a function that is faster than GetBlastHit1 on large scale data
# GetBlastHit1s <- function(blast_re, output_csv_file, function_mode){
#   # initial input datasets
#   colnames(blast_re) <- c("query", "subject", "identity_perc", "alignment_len",
#                           "mismatch_num", "gap_num", "query_start", "query_end",
#                           "subject_start", "subject_end", "expect", "score")
#   target_finish <- NULL
#   # define output dataset
#   blast_re_hit1 <- blast_re[0,]
#   # sub-fraction dataset based on e-values
#   blast_cutoff <- c(1000:1)
#   
#   # determine if subject or query is evaluated
#   if (function_mode == "subject"){
#     target_remain <- unique(blast_re$subject)
#     for (j in 1:length(blast_cutoff)){
#       blast_re_tmp1 <- blast_re[blast_re$score > blast_cutoff[j],]
#       print (paste0(j, "/", length(blast_cutoff), ", ", blast_cutoff[j], ", ", length(unique(blast_re_tmp1$subject))))
#       if (nrow(blast_re_tmp1)){
#         index_tmp <- NULL
#         for (i in 1:length(unique(blast_re_tmp1$subject))){
#           index_tmp1 <- with(blast_re_tmp1, which(subject == unique(subject)[i]))
#           if (length(index_tmp1) == 1){
#             index_tmp <- c(index_tmp, index_tmp1)
#           } else {
#             index_tmp2 <- with(blast_re_tmp1, min(index_tmp1[score[index_tmp1] == max(score[index_tmp1])]))
#             index_tmp <- c(index_tmp, index_tmp2)
#           }
#         }
#         blast_re_tmp2 <- blast_re_tmp1[index_tmp,]
#         target_tmp <- as.vector(blast_re_tmp2$subject)
#         # save the extracted dataset
#         blast_re_hit1 <- rbind.data.frame(blast_re_hit1, blast_re_tmp2)
#         # re-set initial dataset for the next round of computation
#         target_remain <- setdiff(target_remain, target_tmp)
#         target_finish <- c(target_finish, target_tmp)
#         blast_re <- blast_re[blast_re$subject %in% target_remain,]
#       }
#     }
#   } else {
#     if (function_mode == "query"){
#       target_remain <- unique(blast_re$query)
#       for (j in 1:length(blast_cutoff)){
#         blast_re_tmp1 <- blast_re[blast_re$score > blast_cutoff[j],]
#         print (paste0(j, "/", length(blast_cutoff), ", ", blast_cutoff[j], ", ", length(unique(blast_re_tmp1$subject))))
#         if (nrow(blast_re_tmp1)){
#           index_tmp <- NULL
#           for (i in 1:length(unique(blast_re_tmp1$query))){
#             index_tmp1 <- with(blast_re_tmp1, which(query == unique(query)[i]))
#             if (length(index_tmp1) == 1){
#               index_tmp <- c(index_tmp, index_tmp1)
#             } else {
#               index_tmp2 <- with(blast_re_tmp1, min(index_tmp1[score[index_tmp1] == max(score[index_tmp1])]))
#               index_tmp <- c(index_tmp, index_tmp2)
#             }
#           }
#           blast_re_tmp2 <- blast_re_tmp1[index_tmp,]
#           target_tmp <- as.vector(blast_re_tmp2$query)
#           # save the extracted dataset
#           blast_re_hit1 <- rbind.data.frame(blast_re_hit1, blast_re_tmp2)
#           # re-set initial dataset for the next round of computation
#           target_remain <- setdiff(target_remain, target_tmp)
#           target_finish <- c(target_finish, target_tmp)
#           blast_re <- blast_re[blast_re$query %in% target_remain,]
#         }
#       }
#     } else {
#       print ("Function_mode was not recognized!")
#     }
#   }
#   # save
#   write.csv(blast_re_hit1, file = output_csv_file, row.names=F)
# }

### GetPfamHit1: a function that searches for the best Pfam hit for each contig
GetPfamHit1 <- function(pfam_re, output_csv_file){
  pfam_re <- pfam_re[,c(1,3:22)]
  colnames(pfam_re) <- c("contig", "contig_len_aa", "pfam_profile", "pfam_accession", 
                         "pfam_len_aa", "evalue_fseq", "score_fseq", "bias_fseq", 
                         "id_dom", "num_dom", "c_evalue_dom", "i_evalue_dom", 
                         "score_dom", "bias_dom", "pfam_start", "pfam_end", 
                         "contig_start", "contig_end", "env_start", "env_end", "acc")
  index_pfam <- NULL
  for (i in 1:length(unique(pfam_re$contig))){
    index_tmp1 <- with(pfam_re, which(contig == unique(contig)[i]))
    index_tmp2 <- with(pfam_re, min(index_tmp1[score_dom[index_tmp1] == max(score_dom[index_tmp1])]))
    index_pfam <- c(index_pfam, index_tmp2)
    if(i%%100 == 0) print(i)
  }
  data_pfam <- pfam_re[index_pfam,]
  write.csv(data_pfam, file = output_csv_file, row.names = F)
}

#################################################################
### calculate the best BLAST hit for each contig
# define the input blastp file
# sample data are a blastp result, in which benzylisoquinoline alkaloid (bia) 
# proteins are query, while Arabidopsis thaliana (at) proteins are subject
blast_re <- read.delim("data_demo/s0_at2mca.blastp", header=F)

### calculate the best BLAST hit for each query
# define the output csv file
output_csv_file <- "data_demo/s1_at2mca_hit1_query.csv"
# define hit1 for each query ("query") or for each subject (""subject)
function_mode <- "query"
# execute
GetBlastHit1(blast_re, output_csv_file, function_mode)
# GetBlastHit1s(blast_re, output_csv_file, function_mode)

### calculate the best BLAST hit for each subject
# define the output csv file
output_csv_file <- "data_demo/s1_at2mca_hit1_subject.csv"
# define hit1 for each query ("query") or for each subject (""subject)
function_mode <- "subject"
# execute
GetBlastHit1(blast_re, output_csv_file, function_mode)
# GetBlastHit1s(blast_re, output_csv_file, function_mode)

### calculate the best Pfam hit for each subject
# define the input hmmer domtblout file
pfam_re <- read.table("data_demo/s0_mca.domtblout", header=F, fill=T)
# define the output csv file
output_csv_file <- "data_demo/s1_mca_pfam_hit1.csv"
# execute
GetPfamHit1(pfam_re, output_csv_file)
