# Step 7: Heatmap Visualization


# Load required libraries
library(pheatmap)
library(RColorBrewer)
library(tidyverse)

# Set working directory
setwd("~/Desktop/Alzheimer_RNAseq_Project")

# Create output directories
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)


# 1. Get normalized expression data

# Apply VST (Variance Stabilizing Transformation) for better visualization
dds_normalized <- vst(dds)
normalized_counts <- assay(dds_normalized)

cat("Normalized counts matrix dimensions:", dim(normalized_counts), "\n")


# 2. Select top 30 genes by adjusted p-value

# Get results as data frame
results_df <- results %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

# Select top 30 genes
top_30_genes <- results_df %>%
  head(30) %>%
  pull(gene_id)

cat("Top 30 genes selected for heatmap\n")
cat("First 10 genes:\n")
print(head(top_30_genes, 10))


# 3. Extract expression matrix for top genes


expr_matrix_top30 <- normalized_counts[as.character(top_30_genes), ]

cat("\nExpression matrix for top 30 genes:\n")
cat("Dimensions:", dim(expr_matrix_top30), "\n")
print(head(expr_matrix_top30))


# 4. Prepare sample annotations


sample_anno <- data.frame(
  row.names = colnames(expr_matrix_top30),
  Condition = metadata$condition
)

cat("\nSample annotations:\n")
print(sample_anno)


# 5. Create heatmap for top 30 genes


cat("\nGenerating heatmap for top 30 genes...\n")

pdf("results/figures/heatmap_top30_genes.pdf", width = 10, height = 14)

pheatmap(expr_matrix_top30,
         annotation_col = sample_anno,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-3, 3, length.out = 101),
         show_colnames = TRUE,
         show_rownames = TRUE,
         fontsize_row = 9,
         fontsize_col = 11,
         fontsize = 12,
         main = "Top 30 Differentially Expressed Genes in Alzheimer's Disease",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete")

dev.off()

cat("Heatmap saved to results/figures/heatmap_top30_genes.pdf\n")


# 6. Create heatmap for top 15 genes (zoomed version)


top_15_genes <- results_df %>%
  head(15) %>%
  pull(gene_id)

expr_matrix_top15 <- normalized_counts[as.character(top_15_genes), ]

cat("\nGenerating heatmap for top 15 genes (zoomed)...\n")

pdf("results/figures/heatmap_top15_genes.pdf", width = 10, height = 10)

pheatmap(expr_matrix_top15,
         annotation_col = sample_anno,
         scale = "row",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-3, 3, length.out = 101),
         show_colnames = TRUE,
         show_rownames = TRUE,
         fontsize_row = 11,
         fontsize_col = 12,
         fontsize = 13,
         main = "Top 15 Differentially Expressed Genes in Alzheimer's Disease",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete")

dev.off()

cat("Heatmap saved to results/figures/heatmap_top15_genes.pdf\n")


# 7. Create heatmap for key Alzheimer's genes


# Define Alzheimer's related genes if they exist in our data
ad_related_genes <- c("60401", "100287102", "653635", "102466751")

# Filter for genes that exist in our data
ad_genes_in_data <- ad_related_genes[ad_related_genes %in% rownames(normalized_counts)]

if(length(ad_genes_in_data) > 0) {
  expr_matrix_ad <- normalized_counts[ad_genes_in_data, ]
  
  cat("\nGenerating heatmap for Alzheimer's related genes...\n")
  cat("Genes found:", ad_genes_in_data, "\n")
  
  pdf("results/figures/heatmap_ad_related_genes.pdf", width = 8, height = 6)
  
  pheatmap(expr_matrix_ad,
           annotation_col = sample_anno,
           scale = "row",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-3, 3, length.out = 101),
           show_colnames = TRUE,
           show_rownames = TRUE,
           fontsize_row = 12,
           fontsize_col = 12,
           fontsize = 13,
           main = "Expression of Key Alzheimer's Related Genes",
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           clustering_distance_cols = "euclidean",
           clustering_method = "complete")
  
  dev.off()
  
  cat("Heatmap saved to results/figures/heatmap_ad_related_genes.pdf\n")
} else {
  cat("No specific Alzheimer's genes found in the data\n")
}

# 8. Create sample dendrogram heatmap


cat("\nGenerating sample clustering heatmap...\n")

# Use all significant genes (pvalue < 0.01) for sample clustering
sig_genes_pvalue <- results_df %>%
  filter(pvalue < 0.01) %>%
  pull(gene_id)

if(length(sig_genes_pvalue) > 0) {
  expr_matrix_sig <- normalized_counts[as.character(sig_genes_pvalue[1:min(50, length(sig_genes_pvalue))]), ]
  
  pdf("results/figures/heatmap_sample_clustering.pdf", width = 9, height = 12)
  
  pheatmap(expr_matrix_sig,
           annotation_col = sample_anno,
           scale = "row",
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-3, 3, length.out = 101),
           show_colnames = TRUE,
           show_rownames = FALSE,
           fontsize_col = 12,
           fontsize = 12,
           main = "Sample Clustering based on Top 50 Significant Genes",
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete")
  
  dev.off()
  
  cat("Sample clustering heatmap saved to results/figures/heatmap_sample_clustering.pdf\n")
}


# 9. Create summary statistics

# Calculate mean expression for Control and AD
control_samples <- colnames(normalized_counts)[metadata$condition == "control"]
ad_samples <- colnames(normalized_counts)[metadata$condition == "AD"]

expr_control_mean <- rowMeans(normalized_counts[as.character(top_30_genes), control_samples])
expr_ad_mean <- rowMeans(normalized_counts[as.character(top_30_genes), ad_samples])

summary_df <- data.frame(
  gene_id = top_30_genes,
  mean_control = expr_control_mean,
  mean_ad = expr_ad_mean,
  fold_change = expr_ad_mean - expr_control_mean
) %>%
  arrange(desc(abs(fold_change)))

write.csv(summary_df, 
          "results/tables/heatmap_top30_genes_summary.csv", 
          row.names = FALSE)

cat("\nSummary statistics saved to results/tables/heatmap_top30_genes_summary.csv\n")


# 10. Print summary


cat("\n")
cat("========================================\n")
cat("Heatmap Visualization Complete\n")
cat("========================================\n")
cat("Generated heatmaps:\n")
cat("1. Top 30 genes heatmap\n")
cat("2. Top 15 genes heatmap (zoomed)\n")
cat("3. Alzheimer's related genes heatmap\n")
cat("4. Sample clustering heatmap\n")
cat("5. Summary statistics table\n")
cat("\nAll files saved to results/ directory\n")
cat("========================================\n")

# 11. Print interpretation tips

cat("\nInterpretation Guide:\n")
cat("- Blue color: Lower expression levels\n")
cat("- White color: Middle expression levels\n")
cat("- Red color: Higher expression levels\n")
cat("- Left dendrogram: Gene clustering (similar genes grouped together)\n")
cat("- Top dendrogram: Sample clustering (similar samples grouped together)\n")
cat("========================================\n")