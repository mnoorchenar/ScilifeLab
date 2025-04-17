# 🧬 Breast Cancer Gene Network Analysis Pipeline

This pipeline analyzes breast cancer RNA-Seq data to identify key regulatory genes through a combination of gene co-expression network construction, network-based feature engineering, and sparsity-driven unsupervised feature selection (SPNFSR). The pipeline is modular, tunable, and fully traceable through SQLite databases.

---

## 📂 Project Structure

- `functions.R`: Modular functions for preprocessing, network construction, feature extraction, and ranking  
- `main.R`: Master script to run the entire pipeline step-by-step  
- `Data/`: Input files (TSV) and generated output databases (SQLite)  
- `SPNFSR_Tuning_Results_*.pdf`: Visual summaries of parameter tuning results  

---

## 📦 Requirements

Install the following packages before running the pipeline:

```r
install.packages(c(
  "igraph", "RSQLite", "DBI", "Matrix", "HiClimR", "ggplot2",
  "pheatmap", "RColorBrewer", "cluster", "factoextra", "patchwork", "scales", "dplyr"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biomaRt")
```

---

## 🔄 Pipeline Overview

### 🧼 Step 1: Data Preparation
- Filters for protein-coding genes using Ensembl v104
- Maps Ensembl IDs to HGNC gene symbols
- Removes duplicates and version suffixes
- Stores cleaned expression matrix to SQLite

### 🔗 Step 2: Gene Co-expression Network Construction
- Calculates gene–gene correlation using `fastCor()`
- Converts correlations to mutual information (MI)
- Applies a user-defined MI threshold (e.g., 0.7 or 0.8)
- Saves the resulting edge list to SQLite

### 📊 Step 3: Node Feature Extraction
- Computes 9 features per gene:
  - **Network**: PageRank, Betweenness, Closeness, Eigenvector, Degree, Strength, Entropy
  - **Expression**: Mean and SD
- Optionally restricts analysis to the largest connected component
- Saves results to a features table in SQLite

### 🚀 Step 4: SPNFSR (Sparse Projection for Feature Selection)
- Constructs an RBF-based similarity graph
- Applies SPNFSR with grid search over α, β, k, and σ²
- Evaluates unsupervised clustering quality using silhouette score
- Ranks genes based on their projection weights
- Saves both rankings and all parameter tuning results to SQLite

### 📈 Step 5: Visualization
- Generates a summary PDF with:
  - Top 15 configurations by silhouette score
  - Influence of each parameter (α, β, k, σ²) on performance

### 🧠 Step 6: Top Gene & Neighborhood Summary
- For the top-ranked genes:
  - Lists neighbors in the gene network
  - Reports mean expression levels and number of connections
- Useful for biological interpretation and validation

### 🔍 Step 7: Network Connectivity Check
- Analyzes how connected the gene network is under a given MI threshold
- Reports:
  - Number of connected components
  - Size and percentage of the largest component
  - Total nodes and edges
- Helps assess how MI threshold affects global network structure

---

## ▶️ Running the Pipeline

Run the pipeline by modifying and executing `main.R`. Set the desired mutual information threshold and call each function step-by-step using your selected input and database paths. Output tables are stored in SQLite for ease of inspection and reuse.

---

## 📌 Summary

This pipeline integrates network theory and unsupervised feature learning to identify key driver genes from RNA-Seq data in breast cancer. It is designed to be modular, parameter-tunable, and database-driven for reproducibility and downstream analysis.

The method highlights how network fragmentation varies with MI threshold, and how topological and expression-based features can be used to prioritize candidate genes in an interpretable, data-driven way.
