#!/usr/bin/env Rscript
# @Author: Rogério Eduardo Ramos Ribeiro
# @E-mail: rogerio.e.ramos.ribeiro
# @Description: Code to analyse the PSI data
# @software version: R=4.4

#Biological question: Is there more alternative splicing for KRAS inhibited cell lines

library(data.table)
library(tidyverse)
library(ggplot2)

# Load data----
# Load Transcripts PSI
transcript_psi <- fread("SUPPA/psiPerIsoform_isoform.psi")
# Load biotype
biotype <- read.csv("genome/biotype.csv", header = F)
# Load expression data
expression_data_per_gene <- fread("quant/tpm_genes.csv")
expression_data_per_isofrm <- fread("quant/tpm_isoform.csv")


# Transform into a transcript - gene psi table
transcript_psi <- transcript_psi %>%
  separate(V1, into = c("Gene_ID", "Transcript_ID"), sep = ";")

# Filter to keep only genes in the pcod and lincRNA 

biotype.pcod <- biotype %>% filter(V2 == "protein_coding")

transcript_psi.filtered <- transcript_psi %>% 
  filter(Gene_ID %in% biotype.pcod$V1)

# Formula for SHANNON entropy
# H = sum(p(x) * log2(px))
# If p(x) = 0; 0 * log2(0) = 0
shannon <- function(vector) -sum(vector * log2(vector) / log2(length(vector)), na.rm =  TRUE)


# Analysis Pipeline:

# 1 - Subset analysis samples (4 cellines, each with control vs KRAS inhibited). 2 cell lines are KRAS dependent, 2 are indepedent (HCT116 + SW480 are KRAS dependent, LS174T and SW837 are KRAS independent)
# KRAS independent and KRAS dependent lines can be concatenated in the same analysis
# 2 - Filter for genes with minimal expression in the samples (transcript with minimal expression were already pre-filtered)
# 3 - Remove genes with only 1 transcript
# 4 - Compute average PSI per genes in the group
# 5 - Compute Shannon in both groups. Compare with the formula: shannon / log(D) where D is the number of isoforms
# 6 - Test if Shannon distribution deviates from 0. Since the sample size is to large maybe a simulation is a better idea.

transcript_shannon_comparison <- function(transcript_psi.filtered, expression_data_per_gene, expression_data_per_isofrm, set, min.gen.exp){
  # Filter samples
  samples_ctrl <- comparison_to_run[[set]][["ctrl"]]
  samples_krinb <- comparison_to_run[[set]][["krasInhibited"]]
  
  transcript_psi.filtered_samples  <- transcript_psi.filtered %>% select(Gene_ID, Transcript_ID, all_of(c(samples_ctrl, samples_krinb)))
  expression_data_per_gene.filtered_samples <- expression_data_per_gene %>% select(V1, all_of(c(samples_ctrl, samples_krinb))) %>% column_to_rownames("V1")
  expression_data_per_isofrm.filtered_samples <- expression_data_per_isofrm %>% select(V1, all_of(c(samples_ctrl, samples_krinb)))  %>% column_to_rownames("V1")
  
  # Filter genes with minimal expression
  keep_genes <- row.names(expression_data_per_gene.filtered_samples)[rowSums(expression_data_per_gene.filtered_samples >= min.gen.exp) > ncol(expression_data_per_gene.filtered_samples)/ 2]
  
  transcript_psi.filtered_genes_isoforms <- transcript_psi.filtered_samples %>% 
    filter(Gene_ID %in% keep_genes)
  
  # Remove genes with only 1 transcript
  n_transcript_per_gene <- transcript_psi.filtered_genes_isoforms %>% 
    group_by(Gene_ID) %>% 
    summarise(n_transcript = n())
  
  genes_with_multiple_isoforms <- n_transcript_per_gene %>% 
    filter(n_transcript != 1)
  
  
  transcript_psi.filtered_single_isoform_genes <- transcript_psi.filtered_genes_isoforms %>% 
    filter(Gene_ID %in% genes_with_multiple_isoforms$Gene_ID)
  
  # Get mean PSI across control samples and KrasSamples
  transcript_psi.ctrl_means <- rowMeans(transcript_psi.filtered_single_isoform_genes[,samples_ctrl], na.rm = TRUE)
  transcript_psi.ctrl_kras <- rowMeans(transcript_psi.filtered_single_isoform_genes[,samples_krinb], na.rm = TRUE)
  
  transcript_psi_mean <- data.frame(transcript_psi.filtered_single_isoform_genes[,c("Gene_ID", "Transcript_ID")], 
                                    "ctr_mean" = transcript_psi.ctrl_means, 
                                    "kras_mean" = transcript_psi.ctrl_kras)
  
  # Compute SHANNON per condition
  
  entropy_per_condition <- transcript_psi_mean %>%
    pivot_longer(cols = c("ctr_mean", "kras_mean"), names_to = "condition", values_to = "psi") %>% 
    group_by(Gene_ID, condition) %>%
    summarize(entropy = shannon(psi), .groups = 'drop') %>% 
    pivot_wider(names_from = condition, values_from = entropy)
  
  relative_change_in_entropy <- entropy_per_condition$kras_mean - entropy_per_condition$ctr_mean
  test <- wilcox.test(relative_change_in_entropy)
  mu <- mean(relative_change_in_entropy)
  
  return(list("list" = entropy_per_condition, "mu" = mu, 'test' = test))
}


# Comparisons to run
comparison_to_run <- list(
  "HCT116" = list(
    "ctrl" = c("HCT116_siCTRL1_S1_R1", "HCT116_siCTRL2_S2_R1"),
    "krasInhibited" = c("HCT116_siKRAS1_S3_R1", "HCT116_siKRAS2_S4_R1")
  ),
  "SW480" = list(
    "ctrl" = c("SW480_siCTRL1_S36_R1", "SW480_siCTRL2_S37_R1", "SW480_siCTRL3_S38_R1"),
    "krasInhibited" = c("SW480_siKRAS1_S39_R1", "SW480_siKRAS2_S40_R1", "SW480_siKRAS3_S41_R1")
  ),
  "LS174T" = list(
    "ctrl" = c("LS174T_siCTRL1_S5_R1", "LS174T_siCTRL2_S6_R1"),
    "krasInhibited" = c("LS174T_siKRAS1_S7_R1", "LS174T_siKRAS2_S8_R1")
  ),
  "SW837" = list(
    "ctrl" = c("SW837_siCTRL1_S42_R1", "SW837_siCTRL2_S43_R1", "SW837_siCTRL3_S44_R1"),
    "krasInhibited" = c("SW837_siKRAS1_S45_R1", "SW837_siKRAS2_S46_R1", "SW837_siKRAS3_S47_R1")
  ),
  "HCT116_SW480" = list(
    "ctrl" = c("HCT116_siCTRL1_S1_R1", "HCT116_siCTRL2_S2_R1", "SW480_siCTRL1_S36_R1", "SW480_siCTRL2_S37_R1", "SW480_siCTRL3_S38_R1"),
    "krasInhibited" = c("HCT116_siKRAS1_S3_R1", "HCT116_siKRAS2_S4_R1", "SW480_siKRAS1_S39_R1", "SW480_siKRAS2_S40_R1", "SW480_siKRAS3_S41_R1")
  ),
  "LS174T_SW837" = list(
    "ctrl" = c("LS174T_siCTRL1_S5_R1", "LS174T_siCTRL2_S6_R1", "SW837_siCTRL1_S42_R1", "SW837_siCTRL2_S43_R1", "SW837_siCTRL3_S44_R1"),
    "krasInhibited" = c("LS174T_siKRAS1_S7_R1", "LS174T_siKRAS2_S8_R1", "SW837_siKRAS1_S45_R1", "SW837_siKRAS2_S46_R1", "SW837_siKRAS3_S47_R1")
  ),
  "all" = list(
    "ctrl" = c("HCT116_siCTRL1_S1_R1", "HCT116_siCTRL2_S2_R1", "SW480_siCTRL1_S36_R1", "SW480_siCTRL2_S37_R1", "SW480_siCTRL3_S38_R1", "LS174T_siCTRL1_S5_R1", "LS174T_siCTRL2_S6_R1", "SW837_siCTRL1_S42_R1", "SW837_siCTRL2_S43_R1", "SW837_siCTRL3_S44_R1"),
    "krasInhibited" = c("HCT116_siKRAS1_S3_R1", "HCT116_siKRAS2_S4_R1", "SW480_siKRAS1_S39_R1", "SW480_siKRAS2_S40_R1", "SW480_siKRAS3_S41_R1", "LS174T_siKRAS1_S7_R1", "LS174T_siKRAS2_S8_R1", "SW837_siKRAS1_S45_R1", "SW837_siKRAS2_S46_R1", "SW837_siKRAS3_S47_R1")
  )
)


gene_table <- list()
results <- data.frame()

for (test in names(comparison_to_run)){
  cat("Running: ", test)
  cat("...Sample size is ", length(comparison_to_run[[test]][["ctrl"]]) + length(comparison_to_run[[test]][["krasInhibited"]]))
  
  res <- transcript_shannon_comparison(transcript_psi.filtered, expression_data_per_gene, expression_data_per_isofrm, test, 1)
  gene_table[[test]] <- res$list
  
  results <- rbind(results, c(test, res$mu, res$test$p.value))
  cat("\n")
}

colnames(results) <- c("comparison","mean_shannon", "pval")
results$FDR <- p.adjust(as.numeric(results$pval))

if (!dir.exists("plots/")){
  dir.create("plots")
}


# Making some plots
for (test in names(gene_table)){
  data <- gene_table[[test]]
  mu <- results %>% filter(comparison == test) %>% pull(mean_shannon)
  pvalue <- results %>% filter(comparison == test) %>% pull(pval)
  

  # Reshape data to longer format
  data_long <- data %>%
    pivot_longer(cols = c("ctr_mean", "kras_mean"), names_to = "condition", values_to = "entropy")
  
  # Plot
  p1 <- ggplot(data_long, aes(x = condition, y = entropy)) +
    geom_violin(trim = FALSE, fill = "skyblue", color = "black") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "red") +
    labs(title = paste0("Violin Plot of Entropy: ", test),
         x = "Condition",
         y = "Entropy") +
    theme_minimal() +
    annotate("text", x = 1.5, y = max(data_long$entropy), 
             label = paste("Difference: ", round(as.numeric(mu), 3), "\nP-value: ", format.pval(as.numeric(pvalue), digits = 3)))
  
  pdf(paste0("plots/shannon_", test, ".pdf"))
  print(p1)
  dev.off()
  
  # Number of genes with higher entropy
  cat("There are ",   sum(data$kras_mean > data$ctr_mean), " (", 100 *  round(sum(data$kras_mean > data$ctr_mean) / nrow(data), 2),"%)"," genes with higher entropy in KrasInhibited in ", test, sep  = "")
  cat("\n")
}


## Due to the high number of genes being tested I do not trust the significance of the results
# Randomly sample entropy across all genes in all conditions. 
average_number_of_genes <- round(mean(do.call("c", lapply(gene_table, function(x) c(nrow(x))))), 0)

entropy_distribution <- do.call("c", lapply(gene_table, function(x) c(x$ctr_mean, x$kras_mean)))
entropy_distribution <- unname(entropy_distribution)

set.seed(42) # For reproducteble results

number_of_iterations <- 1000
pvalue_treshold <- 0.05
average_diference <- c()

for (i in 1:number_of_iterations){
  control_sample <- sample(entropy_distribution, average_number_of_genes)
  kras_inhibited_sample <- sample(entropy_distribution, average_number_of_genes)
 
  mean_average <- mean(kras_inhibited_sample - control_sample)
  if (mean_average < 0 ){
    mean_average = -mean_average
  }
  
  average_diference <- c(average_diference, mean_average)
}


max_diference <- max(average_diference)
max_statistical_significance_treshold <- average_diference[order(average_diference, decreasing = T)][pvalue_treshold*number_of_iterations]



results$comparison <- factor(results$comparison, levels = results$comparison[7:1])

## Plot the absolute mean difference in each analysis
pdf(paste0("plots/shannon_mean_with_simulation_treshdold.pdf"))
ggplot(results, aes(y = comparison, x = abs(as.numeric(mean_shannon)))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_vline(xintercept = max_diference, color = "red", linetype = "dashed", size = 1.5) +
  geom_vline(xintercept = max_statistical_significance_treshold, color = "black", linetype = "dotted", size = 1.5) +
  labs(title = "Bar Plot of Mean Shannon by Comparison",
       x = "Mean Shannon",
       y = "Comparison") +
  theme_minimal() +
  theme(
    text = element_text(size = 14, color = "black"),  # Increase text size and set color to black
    axis.title = element_text(size = 16, color = "black"),# Increase axis title size and set color to blac
    axis.text = element_text(size = 16, colour = "black"),
    plot.title = element_text(size = 18, hjust = 0.5, color = "black")  # Increase plot title size, center it, and set color to black
  )
dev.off()
