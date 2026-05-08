# Step 5: DESeq2 Differential Expression Analysis
# Load data and setup
# Load required libraries
library(dplyr)
library(ggplot2)
setwd("~/Desktop/Alzheimer_RNAseq_Project/data")

counts <- read.delim("GSE53697_raw_counts_GRCh38.p13_NCBI.tsv", 
                     row.names=1)

setwd("~/Desktop/Alzheimer_RNAseq_Project")

metadata <- data.frame(
  sample_id = colnames(counts),
  condition = c(
    "control", "control", "control", "control",
    "control", "control", "AD", "AD",
    "AD", "AD", "AD", "AD"
  ),
  stringsAsFactors = FALSE
)

# Filter low expression genes
counts_filtered <- counts[rowMeans(counts) >= 5, ]

# STEP 5.1: Run DESeq2 

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = metadata,
  design = ~ condition
)

dds <- DESeq(dds)
results <- results(dds)

print("DESeq2 analysis complete")
print(head(results))

# STEP 5.2: Filter significant DEGs 

sig_genes <- results[
  !is.na(results$padj) & 
    results$padj < 0.05 & 
    abs(results$log2FoldChange) > 0.3,
]

print(paste("Number of significant DEGs:", nrow(sig_genes)))
print(paste("Upregulated genes:", sum(sig_genes$log2FoldChange > 0)))
print(paste("Downregulated genes:", sum(sig_genes$log2FoldChange < 0)))

# STEP 5.3: Save DEG results 

results_df <- as.data.frame(results)
results_df$gene_id <- rownames(results_df)
results_df <- results_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

write.csv(results_df, "results/01_all_results.csv", row.names = FALSE)
print("All results saved to results/01_all_results.csv")

sig_genes_df <- as.data.frame(sig_genes)
sig_genes_df$gene_id <- rownames(sig_genes_df)
sig_genes_df <- sig_genes_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

write.csv(sig_genes_df, "results/02_significant_degs.csv", row.names = FALSE)
print("Significant DEGs saved to results/02_significant_degs.csv")

# STEP 5.4: Create Volcano Plot 

library(ggplot2)

volcano_data <- as.data.frame(results) %>%
  mutate(
    significant = padj < 0.05 & abs(log2FoldChange) > 0.3,
    direction = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated")
  )

p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), 
                                      color = significant, shape = direction)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
  geom_vline(xintercept = c(-0.3, 0.3), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  labs(
    title = "Volcano Plot: Alzheimer's Disease vs Control",
    x = "log2(Fold Change)",
    y = "-log10(Adjusted p-value)",
    color = "Significant",
    shape = "Direction"
  )

print(p_volcano)

ggsave("plots/04_volcano_plot.png", plot = p_volcano, width = 10, height = 8)
print("Volcano plot saved to plots/04_volcano_plot.png")

# STEP 5.5: Create MA Plot 

p_ma <- ggplot(as.data.frame(results), 
               aes(x = log10(baseMean), y = log2FoldChange)) +
  geom_point(alpha = 0.6, size = 2, 
             color = ifelse(results$padj < 0.05 & abs(results$log2FoldChange) > 0.3, 
                            "red", "gray")) +
  geom_hline(yintercept = c(-0.3, 0.3), linetype = "dashed", alpha = 0.5) +
  theme_minimal() +
  theme(text = element_text(size = 12)) +
  labs(
    title = "MA Plot: Alzheimer's Disease vs Control",
    x = "log10(Mean Expression)",
    y = "log2(Fold Change)"
  )

print(p_ma)

ggsave("plots/05_ma_plot.png", plot = p_ma, width = 10, height = 8)
print("MA plot saved to plots/05_ma_plot.png")

print("Step 5 complete!")