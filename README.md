# Alzheimer's Disease RNA-seq Analysis: Neuroinflammatory Hub Proteins via PPI Network

## Project Overview
Bulk RNA-seq analysis (GEO: GSE53697) identifying disease-associated hub proteins in protein-protein interaction networks from Alzheimer's disease patients.

## Dataset
- **Accession:** GSE53697 - RNAseq in Alzheimer's Disease patients
- **Tissue:** Dorsolateral prefrontal cortex (dlPFC)
- **Source:** Mount Sinai Brain Bank
- **Samples:** 8 control + 9 AD
- **Platform:** Illumina HiSeq 2500

## Analysis Pipeline

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_data_loading_qc.R` | Load counts, exploratory analysis, sample QC |
| 2 | `02_filtering_and_pca.R` | Gene filtering (CoV), normalization, PCA |
| 3 | `03_deseq2_analysis.R` | Differential expression with batch correction |
| 4 | `04_go_kegg_enrichment.R` | GO/KEGG pathway analysis |
| 5 | `05_heatmap_visualization.R` | Top DEGs heatmap, sample clustering |
| 6 | `06_string_ppi_network.R` | STRING PPI mapping, hub gene identification |

## Key Findings

### Hub Proteins Identified
| Gene | Degree | log2FC | Function |
|------|--------|--------|----------|
| **CX3CR1** | 2 | -0.91 | Chemokine receptor; neuroinflammation |
| **CLEC9A** | 2 | -1.42 | Dendritic cell receptor; antigen presentation |
| SAMD11 | 2 | -0.06 | Regulatory protein |
| KLHL17 | 2 | +0.04 | Regulatory protein |

### Biological Interpretation
**CX3CR1 downregulation** impairs microglial activation and Aβ clearance (Liu et al., Aging Cell 2025).  
**CLEC9A downregulation** compromises dendritic cell antigen presentation to T cells.

**Conclusion:** AD involves immune dysregulation characterized by impaired immune clearance, not excessive inflammation.

## Results
- `results/figures/` - All visualizations (heatmaps, volcano plots, PPI networks, hub gene barplot)
- `results/hub_genes.csv` - Hub gene metrics (degree, centrality, log2FC)
- `results/network_centrality_measures.csv` - Network statistics

## Author
Kitty (Xinyue Hu)  
Master's in Computer Science, University of Tulsa  
GitHub: [@xih0320](https://github.com/xih0320)  
Email: xinyuehu100@gmail.com

## References
- Scheckel C, et al. Elife. 2016;5:e21426. PMID: 26894958
- Liu PP, et al. Aging Cell. 2025 Feb. doi: 10.1111/acel.14393
  
