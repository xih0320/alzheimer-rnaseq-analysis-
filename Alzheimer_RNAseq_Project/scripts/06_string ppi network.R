# Step 6: Protein-Protein Interaction (PPI) Network Analysis
# Using STRING Database for top DEGs
# Install required packages if needed
# BiocManager::install("STRINGdb")
library(STRINGdb)
library(igraph)
library(tidyverse)
library(ggplot2)

# Set working directory
setwd("~/Desktop/Alzheimer_RNAseq_Project")

# Section 1: Prepare gene list from DESeq2 results

# Load the DESeq2 results (should be from 03_deseq2_analysis.R)
# If not already loaded, uncomment and run:
# load("results/deseq2_results.RData")

# Get top DEGs (e.g., top 50 significant genes by padj)
top_n_genes <- 50

results_df <- results %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(top_n_genes)

# Get gene symbols (convert from Entrez IDs if needed)
library(org.Hs.eg.db)

gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = results_df$gene_id,
                       keytype = "ENTREZID",
                       column = "SYMBOL",
                       multiVals = "first")

gene_df <- data.frame(
  gene_id = results_df$gene_id,
  gene_symbol = gene_symbols,
  padj = results_df$padj,
  log2FoldChange = results_df$log2FoldChange,
  stringsAsFactors = FALSE
)

# Remove rows with NA gene symbols
gene_df <- gene_df %>%
  filter(!is.na(gene_symbol)) %>%
  arrange(padj)

print(paste("Total genes for PPI analysis:", nrow(gene_df)))
print("Top 10 genes:")
print(head(gene_df, 10))

# Save gene list
write.csv(gene_df, "results/top_genes_for_ppi.csv", row.names = FALSE)

# Section 2: Query STRING Database

print("Initializing STRING database connection...")

# Initialize STRING database for human (species = 9606)
string_db <- STRINGdb$new(version = "12.0",
                          species = 9606,
                          score_threshold = 400,  # Medium confidence (0-1000 scale)
                          input_directory = "")

# Map gene symbols to STRING protein IDs
mapped_genes <- string_db$map(gene_df,
                              removeUnmappedRows = TRUE,
                              my_data_frame_id_col_names = "gene_symbol")

print(paste("Mapped genes:", nrow(mapped_genes), "out of", nrow(gene_df)))
print("Mapped data:")
print(head(mapped_genes))

# Save mapped genes
write.csv(mapped_genes, "results/mapped_genes_to_string.csv", row.names = FALSE)

# Section 3: Extract PPI Network


print("Extracting protein-protein interactions...")

# Get interaction data from STRING
if(nrow(mapped_genes) > 1) {
  interactions <- string_db$get_interactions(mapped_genes$STRING_id)
  
  print(paste("Number of protein-protein interactions:", nrow(interactions)))
  
  if(nrow(interactions) > 0) {
    print("Sample interactions:")
    print(head(interactions, 10))
    
    # Save interaction data
    write.csv(interactions, "results/ppi_interactions.csv", row.names = FALSE)
  } else {
    print("Warning: No interactions found. May need lower score threshold.")
  }
} else {
  print("Not enough genes to create network. Need at least 2 genes.")
  interactions <- NULL
}

# Section 4: Build Network Graph

if(!is.null(interactions) && nrow(interactions) > 0) {
  
  print("Building interaction network...")
  #Convert interactions to dataframe(interactions)
  interactions_df <- as.data.frame(interactions)
  
  # Create edge list from interactions
  edge_list <- data.frame(
    from = interactions_df$from,
    to = interactions_df$to,
    combined_score = interactions_df$combined_score,
    stringsAsFactors = FALSE
  )
  
  # Create igraph object
  g <- graph_from_data_frame(d = edge_list,
                             directed = FALSE,
                             vertices = NULL)
  
  # Add gene annotations to vertices
  V(g)$gene_symbol <- mapped_genes$gene_symbol[match(names(V(g)), mapped_genes$STRING_id)]
  V(g)$padj <- mapped_genes$padj[match(names(V(g)), mapped_genes$STRING_id)]
  V(g)$log2FC <- mapped_genes$log2FoldChange[match(names(V(g)), mapped_genes$STRING_id)]
  
  print(paste("Network summary:"))
  print(paste("  Nodes:", vcount(g)))
  print(paste("  Edges:", ecount(g)))
  print(paste("  Density:", round(edge_density(g), 4)))
  

  # Section 5: Calculate Network Centrality Measures
  
  print("Calculating network centrality measures...")
  
  # Degree centrality
  V(g)$degree <- degree(g)
  
  # Betweenness centrality
  V(g)$betweenness <- betweenness(g)
  
  # Closeness centrality
  V(g)$closeness <- closeness(g)
  
  # Eigenvector centrality (if graph is connected)
  if(is.connected(g)) {
    V(g)$eigenvector <- evcent(g)$vector
  } else {
    V(g)$eigenvector <- 0
  }
  
  # Create node attributes dataframe
  node_centrality <- data.frame(
    protein_id = names(V(g)),
    gene_symbol = V(g)$gene_symbol,
    padj = V(g)$padj,
    log2FC = V(g)$log2FC,
    degree = V(g)$degree,
    betweenness = V(g)$betweenness,
    closeness = V(g)$closeness,
    eigenvector = V(g)$eigenvector,
    stringsAsFactors = FALSE
  )
  
  # Rank hub genes by degree
  hub_genes <- node_centrality %>%
    arrange(desc(degree)) %>%
    head(15)
  
  print("Top 15 hub genes by degree:")
  print(hub_genes)
  
  # Save centrality measures
  write.csv(node_centrality, "results/network_centrality_measures.csv", row.names = FALSE)
  write.csv(hub_genes, "results/hub_genes.csv", row.names = FALSE)
  
  # Section 6: Identify Network Modules (Communities)

  print("Detecting functional modules in the network...")
  
  # Use Louvain community detection algorithm
  communities <- cluster_louvain(g)
  
  V(g)$module <- membership(communities)
  
  print(paste("Number of modules detected:", length(communities)))
  print("Module sizes:")
  print(sizes(communities))
  
  # Create module membership dataframe
  module_membership <- data.frame(
    protein_id = names(V(g)),
    gene_symbol = V(g)$gene_symbol,
    module = V(g)$module,
    degree = V(g)$degree,
    stringsAsFactors = FALSE
  ) %>%
    arrange(module, desc(degree))
  
  print("Module membership (showing top genes per module):")
  for(m in 1:length(communities)) {
    cat(paste("\nModule", m, ":\n"))
    module_genes <- module_membership %>%
      filter(module == m) %>%
      head(5)
    print(module_genes)
  }
  
  # Save module membership
  write.csv(module_membership, "results/network_modules.csv", row.names = FALSE)
  
  # Section 7: Visualization - Network Graph

  print("Creating network visualizations...")
  
  # Set layout
  layout <- layout_with_fr(g)  # Fruchterman-Reingold layout
  
  # Prepare node sizes based on degree
  node_size <- V(g)$degree * 2 + 3
  
  # Prepare node colors based on log2FC (red = upregulated, blue = downregulated)
  node_colors <- ifelse(V(g)$log2FC > 0, "red", "blue")
  
  # PDF: Full network
  pdf("results/figures/ppi_network_full.pdf", width = 14, height = 12)
  
  plot(g,
       layout = layout,
       vertex.size = node_size,
       vertex.color = node_colors,
       vertex.frame.color = "black",
       vertex.label = V(g)$gene_symbol,
       vertex.label.cex = 0.6,
       vertex.label.dist = 0,
       edge.width = 0.5,
       edge.color = "#CCCCCC",
       main = "Alzheimer's Disease: Protein-Protein Interaction Network\n(Top 50 DEGs)",
       sub = "Node color: Red=upregulated, Blue=downregulated\nNode size=degree centrality")
  
  # Add legend
  legend("topleft",
         legend = c("Upregulated", "Downregulated"),
         col = c("red", "blue"),
         pch = 19,
         bty = "n",
         cex = 1.2)
  
  dev.off()
  print("Saved: ppi_network_full.pdf")
  

  # Section 8: Visualization - Hub Genes Network

  
  # Create subgraph of top hub genes (degree >= 3)
  hub_threshold <- 3
  hub_nodes <- V(g)[degree >= hub_threshold]
  
  if(length(hub_nodes) > 1) {
    g_hub <- induced_subgraph(g, hub_nodes)
    
    layout_hub <- layout_with_fr(g_hub)
    node_size_hub <- V(g_hub)$degree * 3 + 5
    node_colors_hub <- ifelse(V(g_hub)$log2FC > 0, "#FF6B6B", "#4ECDC4")
    
    pdf("results/figures/ppi_network_hub_genes.pdf", width = 12, height = 10)
    
    plot(g_hub,
         layout = layout_hub,
         vertex.size = node_size_hub,
         vertex.color = node_colors_hub,
         vertex.frame.color = "black",
         vertex.label = V(g_hub)$gene_symbol,
         vertex.label.cex = 0.8,
         vertex.label.dist = 0.3,
         edge.width = 1,
         edge.color = "#999999",
         main = "Hub Genes Network (Degree >= 3)\nAlzheimer's Disease",
         sub = "Highly connected proteins and their interactions")
    
    dev.off()
    print("Saved: ppi_network_hub_genes.pdf")
  }
  

  # Section 9: Visualization - Degree Distribution

  
  pdf("results/figures/degree_distribution.pdf", width = 10, height = 8)
  
  degree_dist <- data.frame(degree = V(g)$degree)
  
  ggplot(degree_dist, aes(x = degree)) +
    geom_histogram(binwidth = 1, fill = "#3498db", color = "black", alpha = 0.7) +
    geom_vline(aes(xintercept = mean(degree)), color = "red", linetype = "dashed", size = 1) +
    labs(title = "Degree Distribution in PPI Network",
         subtitle = "Alzheimer's Disease DEGs",
         x = "Degree (number of interactions)",
         y = "Number of proteins") +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12))
  
  dev.off()
  print("Saved: degree_distribution.pdf")
  

  # Section 10: Visualization - Hub Genes Bar Plot
  
  pdf("results/figures/hub_genes_barplot.pdf", width = 12, height = 8)
  
  hub_plot <- hub_genes %>%
    mutate(gene_symbol = fct_reorder(gene_symbol, degree)) %>%
    ggplot(aes(x = gene_symbol, y = degree, fill = log2FC)) +
    geom_col() +
    scale_fill_gradient2(low = "#4ECDC4", mid = "white", high = "#FF6B6B",
                         limits = c(-max(abs(hub_genes$log2FC)), max(abs(hub_genes$log2FC)))) +
    coord_flip() +
    labs(title = "Top 15 Hub Genes in PPI Network",
         subtitle = "Ranked by degree centrality",
         x = "Gene",
         y = "Degree (number of interactions)",
         fill = "log2FC") +
    theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.text = element_text(size = 11),
          axis.title = element_text(size = 12))
  
  print(hub_plot)
  
  dev.off()
  print("Saved: hub_genes_barplot.pdf")
  

  # Section 11: Summary Report
  
  print("\n================== NETWORK ANALYSIS SUMMARY ==================")
  print(paste("Total proteins in network:", vcount(g)))
  print(paste("Total interactions:", ecount(g)))
  print(paste("Network density:", round(edge_density(g), 4)))
  print(paste("Average degree:", round(mean(V(g)$degree), 2)))
  print(paste("Max degree:", max(V(g)$degree)))
  print(paste("Number of modules:", length(communities)))
  print(paste("Network modularity:", round(modularity(communities), 3)))
  
  print("\nTop 5 hub genes by degree:")
  print(hub_genes[1:5, c("gene_symbol", "degree", "padj", "log2FC")])
  
  print("\nFiles saved:")
  print("  - results/top_genes_for_ppi.csv")
  print("  - results/ppi_interactions.csv")
  print("  - results/network_centrality_measures.csv")
  print("  - results/hub_genes.csv")
  print("  - results/network_modules.csv")
  print("  - results/figures/ppi_network_full.pdf")
  print("  - results/figures/ppi_network_hub_genes.pdf")
  print("  - results/figures/degree_distribution.pdf")
  print("  - results/figures/hub_genes_barplot.pdf")
  
  print("=========================================================\n")
  
} else {
  print("Unable to create network: insufficient interaction data.")
  print("Consider lowering the score threshold or using more genes.")
}
