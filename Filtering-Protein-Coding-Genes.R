if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

library(biomaRt)

# Connect to Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 104)  # ✅ This is where 'version' goes)

# Get list of protein-coding genes
protein_coding_genes <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype"),
  filters = "biotype",
  values = "protein_coding",
  mart = ensembl)

#Since the Ensembl IDs have version numbers (e.g. "ENSG00000141510.14"), you may need to remove them: 
mydata <- read.csv('Data/TCGA-BRCA.star_tpm.tsv', sep = '\t')
  
rownames(mydata) <- sub("\\..*", "", rownames(mydata))
