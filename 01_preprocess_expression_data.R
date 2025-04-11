# ---- Install and load required packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
if (!requireNamespace("RSQLite", quietly = TRUE)) install.packages("RSQLite")
if (!requireNamespace("DBI", quietly = TRUE)) install.packages("DBI")

library(biomaRt)
library(DBI)
library(RSQLite)

# ---- Connect to Ensembl using a specific version ----
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 104)

# ---- Retrieve protein-coding genes from Ensembl ----
protein_coding_genes <- getBM(
  attributes = c("ensembl_gene_id", "gene_biotype"),
  filters = "biotype",
  values = "protein_coding",
  mart = ensembl
)

# ---- Load your gene expression data ----
# First, try loading without setting row.names in case the gene IDs are in a named column
mydata <- read.csv("Data/TCGA-BRCA.star_tpm.tsv", sep = '\t')

# Detect whether gene IDs are already in the first column
gene_col_name <- colnames(mydata)[1]
if (grepl("^ENSG", mydata[[gene_col_name]][1])) {
  
  # Clean Ensembl IDs (remove version numbers)
  mydata[[gene_col_name]] <- sub("\\..*", "", mydata[[gene_col_name]])
  
  # Remove duplicate gene IDs — keep the first occurrence
  mydata <- mydata[!duplicated(mydata[[gene_col_name]]), ]
  
  # Set as rownames
  rownames(mydata) <- mydata[[gene_col_name]]
  mydata[[gene_col_name]] <- NULL
  
  message("✅ Gene IDs cleaned and duplicates removed.")
  
} else {
  stop("❌ Gene IDs do not appear to be in the first column. Please check the file format.")
}


# ---- (Optional) Keep only protein-coding genes ----
# Uncomment if you want to filter down to protein-coding genes
mydata <- mydata[rownames(mydata) %in% protein_coding_genes$ensembl_gene_id, ]

# ---- Move rownames (Gene IDs) into a column for database compatibility ----
mydata$GeneID <- rownames(mydata)
mydata <- mydata[, c(ncol(mydata), 1:(ncol(mydata)-1))]  # GeneID to first column

write.table(mydata, file = './Data/mydata.txt', sep = '\t')

# ---- Save to SQLite ----
con <- dbConnect(RSQLite::SQLite(), "./Data/BRCA_GeneExpression.db")
dbWriteTable(con, "Gene_Expression", mydata, overwrite = TRUE)
dbDisconnect(con)

cat("✅ Gene expression data saved to SQLite.\n")
