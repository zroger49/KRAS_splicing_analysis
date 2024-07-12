#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro
# @Description: Code to combine the RSEM output from all samples into a single file
# @software version: R=4.4

#setwd("D:/PhD/KRAS_splicing"

library(tidyverse)

# Load files one by one
isoform_tpm <- data.frame()

# Save file names
file_names <- c()

# For isoforms
number_of_gene <- 0
for (file in dir("quant/")){
  if (endsWith(file, "isoforms.results")){
    quant <- read.csv(paste0("quant/", file), sep = "\t") %>% select(transcript_id, TPM) %>% column_to_rownames("transcript_id")
    file_name <- str_split(file, "\\.")[[1]][1]
    file_names <- c(file_names, file_name)
    if (number_of_gene == 0){
      isoform_tpm <- quant
      number_of_gene <- nrow(quant)
    }
    else{
      isoform_tpm <- cbind(isoform_tpm, quant)
    }
    
    if (number_of_gene != nrow(quant)){
      print("WARNING : Number of isoforms not the same")
    }
    
    
  }
}

colnames(isoform_tpm) <- file_names

write.table(isoform_tpm, "quant/tpm_isoform.csv", quote = F, sep = "\t")


# For Gene
genes_tpm <- data.frame()
number_of_gene <- 0
file_names <- c()
for (file in dir("quant/")){
  if (endsWith(file, "genes.results")){
    quant <- read.csv(paste0("quant/", file), sep = "\t") %>% select(gene_id, TPM) %>% column_to_rownames("gene_id")
    file_name <- str_split(file, "\\.")[[1]][1]
    file_names <- c(file_names, file_name)
    if (number_of_gene == 0){
      genes_tpm <- quant
      number_of_gene <- nrow(quant)
    }
    else{
      genes_tpm <- cbind(genes_tpm, quant)
    }
    
    if (number_of_gene != nrow(quant)){
      print("WARNING : Number of genes not the same")
    }
    
    
  }
}

colnames(genes_tpm) <- file_names

write.table(genes_tpm, "quant/tpm_genes.csv", quote = F, sep = "\t")
