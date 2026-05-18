# ==============================================================================
# Script: Differential Gene Expression Analysis & Visualizations (DESeq2)
# Laboratory: High-Throughput Transcriptomic Profiling
# Pipeline Stage: Phase 2 (Downstream Visualizations via R)
# ==============================================================================

#------------------------------------------------------------------------------#
# Install required packages (run once)
#------------------------------------------------------------------------------#
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "Rsubread"), ask = FALSE)
install.packages(c("tidyverse", "pheatmap"), dependencies = TRUE)
install.packages("ggrepel")

#------------------------------------------------------------------------------#
# Load libraries
#------------------------------------------------------------------------------#
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)

#------------------------------------------------------------------------------#
# Load data
#------------------------------------------------------------------------------#
counts_df <- read.delim("countdata.txt", header = TRUE, sep = "\t")
rownames(counts_df) <- counts_df$Geneid
counts_df <- counts_df[, -c(1)]  # Remove EntrezGeneID and Length columns

sample_data <- read.delim("factordata.txt", header = TRUE, sep = "\t")

#------------------------------------------------------------------------------#
# Filter lowly expressed genes
#------------------------------------------------------------------------------#
counts_df <- counts_df[rowSums(counts_df >= 10) >= 3, ]  # Fixed: at least 10 reads in >= 3 samples

#------------------------------------------------------------------------------#
# Create DESeq2 object and run analysis
#------------------------------------------------------------------------------#
dds <- DESeqDataSetFromMatrix(
  countData = counts_df,
  colData = data.frame(condition = sample_data$Condition),
  design = ~ condition
)

dds <- DESeq(dds)
res <- as.data.frame(results(dds))

#------------------------------------------------------------------------------#
# Filter significant DEGs (padj < 0.05, |log2FC| > 1)
#------------------------------------------------------------------------------#
sig_genes <- res %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  rownames()

write.table(sig_genes, "DE_genes.txt", quote = FALSE, row.names = FALSE)

#------------------------------------------------------------------------------#
# Gene Annotation (Mapping Ensembl IDs to Gene Symbols)
#------------------------------------------------------------------------------#
library(org.Hs.eg.db)

# Map Ensembl IDs to official Hugo Gene Symbols
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

#------------------------------------------------------------------------------#
# Volcano Plot Data Preparation
#------------------------------------------------------------------------------#
# Classify genes based on log2FC threshold (0.585 = 1.5 fold change)
res <- res %>%
  mutate(regulation = case_when(
    log2FoldChange > 0.585 & padj < 0.05  ~ "Upregulated",
    log2FoldChange < -0.585 & padj < 0.05 ~ "Downregulated",
    TRUE                                  ~ "Not Significant"
  ))

# Plot global transcriptomic landscape
ggplot(res %>% filter(!is.na(padj)), aes(x = log2FoldChange, y = -log10(padj),
                                         color = regulation)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue",
                                "Not Significant" = "grey")) +
  
  # Label highly significant DEGs with gene symbols
  geom_text_repel(
    data = res %>% filter(!is.na(padj) & (log2FoldChange > 2 | log2FoldChange < -2) & padj < 0.01),
    aes(label = symbol), 
    size = 3,
    max.overlaps = 20, 
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black", 
    show.legend = FALSE
  ) +
  
  labs(title = "Volcano Plot with Gene Names", 
       x = "Log2 Fold Change", 
       y = "Adjusted P-value") +
  theme_bw()

#------------------------------------------------------------------------------#
# Heatmaps (Hierarchical Clustering)
#------------------------------------------------------------------------------#
normalized_counts <- counts(dds, normalized = TRUE)

# --- All significant DEGs ---
# Z-score normalization across rows for consistent scaling
degs_scaled <- t(scale(t(normalized_counts[sig_genes, ])))

pheatmap(degs_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Heatmap of All DEGs")

# --- Top 10 DEGs by adjusted p-value ---
top10_genes <- res %>%
  filter(!is.na(padj), padj < 0.05) %>%
  arrange(padj) %>%
  head(10) %>%
  rownames()

top10_scaled <- t(scale(t(normalized_counts[top10_genes, ])))

# Handle missing symbols by replacing them with original Ensembl IDs
top10_symbols <- res[top10_genes, "symbol"]
top10_symbols[is.na(top10_symbols)] <- top10_genes[is.na(top10_symbols)]

my_colors <- colorRampPalette(c("blue", "white", "red"))(50)

# Plot expression profiles for the top 10 most significant genes
pheatmap(top10_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         labels_row = top10_symbols,
         fontsize_row = 10,
         color = my_colors,
         cellwidth = 30,   
         cellheight = 45,
         main = "top 10 DEGs")
