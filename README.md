<h2>🧬 Breast Cancer Gene Network Analysis Pipeline</h2>
<p>This pipeline analyzes RNA-Seq data for breast cancer to identify influential genes using network science and a sparsity-regularized projection technique. It involves data cleaning, network construction, feature extraction, and feature selection with unsupervised parameter tuning.</p>

<h3>📦 Requirements</h3>
<p>Install these R packages before running the pipeline:</p>
<pre><code>
install.packages(c("igraph", "RSQLite", "DBI", "Matrix", "HiClimR", "ggplot2", "pheatmap", "RColorBrewer", "cluster", "factoextra"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("biomaRt")
</code></pre>

<h3>📁 Input Data</h3>
<ul>
  <li><strong>Gene Expression File:</strong> <code>Data/TCGA-BRCA.star_tpm.tsv</code></li>
  <li>Format: Rows = Ensembl gene IDs, Columns = Samples</li>
</ul>

<h3>🔄 Pipeline Steps</h3>

<h4>Step 1: Data Cleaning and Filtering</h4>
<ul>
  <li>Connects to Ensembl v104 using <code>biomaRt</code></li>
  <li>Filters for <strong>protein-coding genes</strong> only</li>
  <li>Removes version numbers from Ensembl gene IDs</li>
  <li>Saves the expression matrix to a local SQLite database</li>
</ul>
<p><strong>✅ Output:</strong> Saved in <code>BRCA_GeneExpression.db</code> under <code>Gene_Expression</code></p>

<h4>Step 2: Gene Network Construction (Mutual Information)</h4>
<ul>
  <li>Computes fast Pearson correlation using <code>HiClimR::fastCor</code></li>
  <li>Converts correlations to <strong>Mutual Information (MI)</strong> using the formula:</li>
  <code>MI(i, j) = -0.5 × log(1 - r<sub>ij</sub><sup>2</sup>)</code>
  <li>Constructs an undirected weighted graph using edges with MI ≥ 0.7</li>
</ul>
<p><strong>✅ Outputs:</strong></p>
<ul>
  <li><code>BRCA_network_edges.txt</code></li>
  <li>Stored in SQLite under <code>Gene_Network_Edges</code></li>
</ul>

<h4>Step 3: Node-Level Feature Extraction</h4>
<ul>
  <li>Calculates 9 features for each gene:</li>
  <ul>
    <li>PageRank</li>
    <li>Betweenness Centrality</li>
    <li>Closeness Centrality</li>
    <li>Eigenvector Centrality</li>
    <li>Degree</li>
    <li>Strength (sum of edge weights)</li>
    <li>Entropy (of connected edge weights)</li>
    <li>Expression Mean</li>
    <li>Expression Standard Deviation</li>
  </ul>
</ul>
<p><strong>✅ Output:</strong> <code>BRCA_node_features.txt</code>, stored in SQLite as <code>Gene_Network_Features</code></p>

<h4>Step 4: Feature Selection via SPNFSR with Unsupervised Tuning</h4>
<ul>
  <li>Constructs a similarity matrix using RBF kernel over Euclidean distances</li>
  <li>Computes graph Laplacian \( L = I - S - S^T + SS^T \)</li>
  <li>Decomposes \( M = X^T L X \) into positive and negative parts</li>
  <li>Iteratively learns sparse feature weights \( W \in \mathbb{R}^{n \times 1} \)</li>
  <li>Projects genes using \( s_i = \| x_i^T W \|^2 \) to rank their importance</li>
  <li>Performs grid search over parameters:</li>
  <ul>
    <li>\( \alpha \in \{0.1, 1, 10\} \)</li>
    <li>\( \beta \in \{0.01, 0.1, 1\} \)</li>
    <li>\( k \in \{3, 5, 7\} \)</li>
    <li>\( \sigma^2 \in \{50, 100, 200\} \)</li>
  </ul>
  <li>Each configuration is evaluated using <strong>silhouette score</strong> from K-means clustering on the projected space</li>
  <li>The best-performing parameter combination is selected for final ranking</li>
</ul>
<p><strong>✅ Outputs:</strong></p>
<ul>
  <li><code>BRCA_ranked_genes_SPNFSR_CV.txt</code>: Ranked gene list</li>
  <li>Stored in SQLite as <code>BRCA_ranked_genes_SPNFSR_CV</code></li>
  <li>Cross-validation results stored as <code>SPNFSR_CV_Results</code></li>
</ul>

<h3>📊 Visualizations (Optional)</h3>
<ol>
  <li><strong>Top 20 Genes:</strong> Horizontal barplot</li>
  <li><strong>Correlation Heatmap:</strong> For all node-level features</li>
</ol>

<h3>📌 Summary</h3>
<p>
  This pipeline identifies <strong>key regulatory genes</strong> from RNA-Seq data using a graph-based, sparsity-driven projection method. By combining mutual information, graph topology, and unsupervised model tuning, the method produces a robust and interpretable ranking of genes that are most influential in the breast cancer co-expression network.
</p>
