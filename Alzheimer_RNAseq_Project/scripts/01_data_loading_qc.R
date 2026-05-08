#Step 1
setwd("~/Desktop/Alzheimer_RNAseq_Project")
counts <- read.delim("GSE53697_raw_counts_GRCh38.p13_NCBI.tsv", 
                     row.names=1)
head(counts[,1:5])
dim(counts)

#Step 2 
metadata <- data.frame(
  sample_id = colnames(counts),
  condition = c(
    "control", "control", "control", "control",
    "control", "control", "AD", "AD",
    "AD", "AD", "AD", "AD"
  ),
  stringsAsFactors = FALSE
)

metadata
table(metadata$condition)