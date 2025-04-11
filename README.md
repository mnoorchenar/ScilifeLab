<h2>🧬 Breast Cancer Gene Network Analysis Pipeline</h2>
<p>This pipeline analyzes RNA-Seq data for breast cancer to identify influential genes using network science and unsupervised learning techniques. It involves data cleaning, network construction, feature extraction, and gene ranking.</p>

<h3>📦 Requirements</h3>
<p>Install these R packages before running the pipeline:</p>
<pre><code>
install.packages(c("igraph", "RSQLite", "DBI", "Matrix", "HiClimR", "ggplot2", "pheatmap", "RColorBrewer"))
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
  <li>Filters for <strong>protein-coding genes</strong> (optional)</li>
  <li>Cleans gene IDs (removes version numbers, removes duplicates)</li>
  <li>Saves expression matrix to SQLite database</li>
</ul>
<p><strong>✅ Output:</strong> Gene expression saved in <code>BRCA_GeneExpression.db</code> under <code>Gene_Expression</code></p>

<h4>Step 2: Gene Network Construction (Mutual Information)</h4>
<ul>
  <li>Computes fast Pearson correlation with <code>HiClimR::fastCor</code></li>
  <li>Converts correlation matrix to <strong>Mutual Information (MI)</strong></li>
  <li>Builds gene network with edges MI ≥ 0.7</li>
  <li>Constructs undirected, weighted graph</li>
</ul>
<p><strong>✅ Outputs:</strong></p>
<ul>
  <li><code>BRCA_network_edges.txt</code></li>
  <li>Saved to SQLite as <code>Gene_Network_Edges</code></li>
</ul>

<h4>Step 3: Node-Level Feature Extraction</h4>
<ul>
  <li>Calculates:</li>
  <ul>
    <li>PageRank</li>
    <li>Betweenness Centrality</li>
    <li>Closeness Centrality</li>
    <li>Entropy (of edge weights)</li>
  </ul>
</ul>
<p><strong>✅ Output:</strong> <code>BRCA_node_features.txt</code>, stored in SQLite as <code>Gene_Network_Features</code></p>

<h4>Step 4: Gene Ranking via SPNFSR-Inspired LS Score</h4>
<ul>
  <li>Standardizes all node features (Z-scores)</li>
  <li>Calculates feature importance weights based on variance:</li>
  <ul>
    <li>
      Let \( Z_{ij} \) be the standardized value of gene \( g_i \)'s \( j \)-th feature.
    </li>
    <li>
      Let \( \sigma_j^2 \) be the variance of feature \( j \).
    </li>
    <li>
      Weight for feature \( j \): \( w_j = \frac{\sigma_j^2}{\sum_k \sigma_k^2} \)
    </li>
    <li>
      LS Score for each gene:
      <br><code>LS(g<sub>i</sub>) = Σ<sub>j</sub> (w<sub>j</sub> × Z<sub>ij</sub>)</code>
    </li>
  </ul>
  <li>Ranks genes by descending LS Score</li>
</ul>

<p><strong>✅ Output:</strong> <code>BRCA_ranked_genes_LS.txt</code></p>

<h3>📌 Summary</h3>
<p>
  This pipeline identifies <strong>key regulatory genes</strong> in a gene interaction network derived from RNA-Seq breast cancer data. It relies on network topology, mutual information, and unsupervised learning — making it ideal for <em>biomarker discovery</em>, <em>drug target prioritization</em>, and <em>systems biology research</em>.
</p>
