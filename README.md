<h2>🧬 Breast Cancer Gene Network Analysis Pipeline</h2>
<p>This pipeline analyzes RNA-Seq data for breast cancer to identify influential genes using network science and a sparsity-regularized projection technique. It includes preprocessing, network construction, feature extraction, feature selection, and visualization.</p>

<h3>📂 Project Structure</h3>
<ul>
  <li><code>functions.R</code> – Contains all modular functions</li>
  <li><code>main.R</code> – Executes the full pipeline step-by-step</li>
  <li><code>Data/</code> – Input and output files (TSV, SQLite)</li>
  <li><code>SPNFSR_Tuning_Results_*.pdf</code> – Auto-generated visual reports</li>
</ul>

<h3>📦 Requirements</h3>
<p>Install these R packages before running the pipeline:</p>
<pre><code>
install.packages(c(
  "igraph", "RSQLite", "DBI", "Matrix", "HiClimR", "ggplot2",
  "pheatmap", "RColorBrewer", "cluster", "factoextra", "patchwork", "scales"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biomaRt")
</code></pre>

<h3>📁 Input Data</h3>
<ul>
  <li><strong>Gene Expression File:</strong> <code>Data/TCGA-BRCA.star_tpm.tsv</code></li>
  <li>Format: Rows = Ensembl Gene IDs, Columns = Sample TPMs</li>
</ul>

<h3>🔄 Pipeline Overview</h3>

<h4>Step 1: 🧼 Data Preparation</h4>
<ul>
  <li>Filters for protein-coding genes using Ensembl v104</li>
  <li>Removes Ensembl version suffixes</li>
  <li>Saves to SQLite: <code>Gene_Expression</code></li>
</ul>

<h4>Step 2: 🔗 Gene Network Construction</h4>
<ul>
  <li>Computes Pearson correlation (HiClimR)</li>
  <li>Converts to Mutual Information (MI)</li>
  <li>Retains edges with MI ≥ threshold (e.g., 0.7)</li>
  <li>Outputs: <code>Gene_Network_Edges_0.7</code> (SQLite)</li>
</ul>

<h4>Step 3: 📊 Node Feature Extraction</h4>
<ul>
  <li>Extracts 9 features per gene node:</li>
  <ul>
    <li>PageRank, Betweenness, Closeness, Eigenvector</li>
    <li>Degree, Strength, Entropy</li>
    <li>Expression Mean, Expression SD</li>
  </ul>
  <li>Outputs: <code>Gene_AllFeatures_0.7</code> (SQLite)</li>
</ul>

<h4>Step 4: 🚀 SPNFSR with Unsupervised Tuning</h4>
<ul>
  <li>Constructs similarity graph (RBF kernel)</li>
  <li>Applies sparsity-constrained projection</li>
  <li>Ranks genes based on projection magnitude</li>
  <li>Grid search over:</li>
  <ul>
    <li><code>α ∈ {0.01, 0.1, 0.5, 1, 5, 10, 25}</code></li>
    <li><code>β ∈ {0.001, 0.01, 0.05, 0.1, 0.5, 1, 2}</code></li>
    <li><code>k ∈ {3, 5, 7, 9}</code></li>
    <li><code>σ² ∈ {25, 50, 100, 150, 200, 300, 500}</code></li>
  </ul>
  <li>Evaluation metric: <strong>Silhouette score</strong></li>
  <li>Outputs:</li>
  <ul>
    <li><code>BRCA_ranked_genes_SPNFSR_CV_0.7</code> (SQLite)</li>
    <li><code>SPNFSR_CV_Results_0.7</code> (SQLite)</li>
  </ul>
</ul>

<h4>Step 5: 📊 Visualization</h4>
<ul>
  <li>Generates combined PDF showing:</li>
  <ul>
    <li>Top 15 configurations (scaled)</li>
    <li>Effect of each parameter on silhouette score</li>
  </ul>
  <li>Output: <code>SPNFSR_Tuning_Results_0.7.pdf</code></li>
</ul>

<h3>▶️ Running the Pipeline</h3>
<p>After setting your input file, run the pipeline step-by-step from <code>main.R</code>:</p>
<pre><code class="r">
source("functions.R")

# Run full analysis for mi_threshold = 0.7
mi_threshold <- 0.7
db_path <- "./Data/BRCA_GeneExpression.db"

prepare_expression_data(input_path = "Data/TCGA-BRCA.star_tpm.tsv", db_path = db_path)
build_gene_network(db_path = db_path, mi_threshold = mi_threshold)
extract_node_features(db_path = db_path, mi_threshold = mi_threshold)
run_spnfsr_cv(db_path = db_path, mi_threshold = mi_threshold)
plot_spnfsr_results(db_path = db_path, mi_threshold = mi_threshold)
</code></pre>

<h3>📌 Summary</h3>
<p>
  This pipeline identifies <strong>key regulatory genes</strong> from RNA-Seq data using a graph-based, sparsity-driven projection method. It integrates network topology and expression patterns with robust unsupervised feature ranking.
</p>
<p>
  Designed for reproducibility, scalability, and interpretability, the method supports threshold tuning and modular re-use across datasets.
</p>
