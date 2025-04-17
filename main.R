# 📥 Load all your defined functions
source("functions.R")

# Set shared parameters
mi_threshold <- 0.7
db_path <- paste0("./Data/BRCA_GeneExpression_",mi_threshold,".db")

# 1️⃣ Prepare expression data and save to SQLite----
prepare_expression_data(
  input_path = "Data/TCGA-BRCA.star_tpm.tsv",
  db_path = db_path,
  table_name = "Gene_Expression",
  keep_protein_coding = TRUE
)

# 2️⃣ Build gene co-expression network from expression data----
build_gene_network(
  db_path = db_path,
  table_expr = "Gene_Expression",
  mi_threshold = mi_threshold,
  edge_table_prefix = "Gene_Network_Edges"
)

# 3️⃣ Extract graph-based and biological features for each gene node -----
extract_node_features(
  db_path = db_path,
  edge_table_prefix = "Gene_Network_Edges",
  expression_table = "Gene_Expression",
  output_table = "Gene_AllFeatures"
)

# 4️⃣ Run SPNFSR with unsupervised feature selection and ranking ----
run_spnfsr_cv(
  db_path = db_path,
  features_table_prefix = "Gene_AllFeatures",
  output_table_prefix = "BRCA_ranked_genes_SPNFSR_CV",
  result_table_prefix = "SPNFSR_CV_Results",
  N = 100  # or increase to 200–300 for better convergence
)

# 5️⃣ Plot the tuning results as a PDF ----
plot_spnfsr_results(
  db_path = db_path,
  mi_threshold = mi_threshold,
  output_pdf = paste0("SPNFSR_Custom_Report",mi_threshold, ".pdf"),
  scale_scores = TRUE,
  top_n = 20,
  show_plot = TRUE,
  pdf_width = 14,
  pdf_height = 10,
  fill_colors = list(
    top15 = "#e69f00",
    alpha = "#56B4E9",
    beta = "#009E73",
    k_thresh = "#CC79A7",
    sigma2 = "#F0E442"
  )
)

# 6️⃣ summarize_top_gene_neighbors ----
summarize_top_gene_neighbors(
  db_path = db_path,
  ranked_table = "BRCA_ranked_genes_SPNFSR_CV",
  edges_table = "Gene_Network_Edges",
  expression_table = "Gene_Expression",
  output_table = "TopGene_Neighbor_Summary",
  top_n = 10
)

# 7️⃣ check_network_connectivity ----
check_network_connectivity(
  mi_threshold = 0.7,
  db_path = db_path,
  edge_table_prefix = "Gene_Network_Edges"
)

check_network_connectivity(
  mi_threshold = 0.8,
  db_path = db_path,
  edge_table_prefix = "Gene_Network_Edges"
)
