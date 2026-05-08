# Step 6: GO/KEGG Enrichment Analysis
# Load required libraries
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)
library(ggplot2)

# Set working directory
setwd("~/Desktop/Alzheimer_RNAseq_Project")

# Create output directories
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)


# 1. Prepare gene list from DESeq2 results


results_df <- results %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  filter(!is.na(padj)) %>%
  dplyr::select(gene_id, log2FoldChange, padj, pvalue)

results_df$gene_id <- as.numeric(results_df$gene_id)

head(results_df)
cat("Total genes analyzed:", nrow(results_df), "\n")


# 2. Extract significant genes


sig_genes <- results_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 0.3) %>%
  pull(gene_id)

cat("Number of significant genes (padj < 0.05, |log2FC| > 0.3):", length(sig_genes), "\n")


# 3. GO Enrichment - Biological Process


ego_bp <- enrichGO(gene = sig_genes,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pvalueCutoff = 0.05,
                   readable = TRUE)

cat("Number of significant GO BP terms:", nrow(ego_bp@result), "\n")

if(nrow(ego_bp@result) > 0) {
  cat("\nTop GO BP terms:\n")
  print(head(ego_bp@result[, c("Description", "pvalue", "p.adjust", "Count")], 10))
  
  pdf("results/figures/GO_BP_enrichment.pdf", width = 12, height = 8)
  print(barplot(ego_bp, showCategory = 20))
  dev.off()
  
  cat("GO BP enrichment plot saved\n")
} else {
  cat("No significant GO BP terms found\n")
}


# 4. GSEA (Gene Set Enrichment Analysis)


geneList <- setNames(results_df$log2FoldChange, results_df$gene_id)
geneList <- sort(geneList, decreasing = TRUE)

cat("\nTotal genes in ranked list:", length(geneList), "\n")

gse_go <- gseGO(geneList = geneList,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05)

cat("Number of significant GSEA pathways:", nrow(gse_go@result), "\n")

if(nrow(gse_go@result) > 0) {
  cat("\nTop GSEA pathways:\n")
  print(head(gse_go@result[, c("Description", "pvalue", "p.adjust", "NES")], 10))
  
  num_plots <- min(6, nrow(gse_go@result))
  pdf("results/figures/GSEA_GO_result.pdf", width = 12, height = 10)
  print(gseaplot2(gse_go, geneSetID = 1:num_plots, base_size = 10))
  dev.off()
  
  cat("GSEA plot saved\n")
} else {
  cat("No significant GSEA pathways found\n")
}


# 5. KEGG Pathway Enrichment


kk <- enrichKEGG(gene = sig_genes,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)

cat("Number of significant KEGG pathways:", nrow(kk@result), "\n")

if(nrow(kk@result) > 0) {
  cat("\nTop KEGG pathways:\n")
  print(head(kk@result[, c("Description", "pvalue", "p.adjust", "Count")], 10))
  
  pdf("results/figures/KEGG_enrichment.pdf", width = 12, height = 8)
  print(barplot(kk, showCategory = 20))
  dev.off()
  
  cat("KEGG enrichment plot saved\n")
} else {
  cat("No significant KEGG pathways found\n")
}


# 6. Save results to CSV files


if(nrow(ego_bp@result) > 0) {
  write.csv(ego_bp@result, 
            "results/tables/GO_BP_enrichment_results.csv", 
            row.names = FALSE)
}

if(nrow(gse_go@result) > 0) {
  write.csv(gse_go@result, 
            "results/tables/GSEA_GO_results.csv", 
            row.names = FALSE)
}

if(nrow(kk@result) > 0) {
  write.csv(kk@result, 
            "results/tables/KEGG_enrichment_results.csv", 
            row.names = FALSE)
}


# 7. Print summary


cat("\n")
cat("========================================\n")
cat("GO/KEGG Enrichment Analysis Complete\n")
cat("========================================\n")
cat("Significant DEGs analyzed:", length(sig_genes), "\n")
cat("GO BP terms found:", nrow(ego_bp@result), "\n")
cat("GSEA pathways found:", nrow(gse_go@result), "\n")
cat("KEGG pathways found:", nrow(kk@result), "\n")
cat("========================================\n")
cat("Results saved to results/tables/ and results/figures/\n")
cat("========================================\n")