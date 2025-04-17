# ---- Function 1: prepare_expression_data ----
prepare_expression_data <- function(
    input_path = "Data/TCGA-BRCA.star_tpm.tsv",
    db_path = "./Data/BRCA_GeneExpression.db",
    table_name = "Gene_Expression",
    keep_protein_coding = TRUE
) {
  # 📦 Required Packages 
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
  if (!requireNamespace("RSQLite", quietly = TRUE)) install.packages("RSQLite")
  if (!requireNamespace("DBI", quietly = TRUE)) install.packages("DBI")
  
  library(biomaRt)
  library(DBI)
  library(RSQLite)
  
  # 🌐 Connect to Ensembl
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 104)
  
  # 📥 Get Gene Info (optionally filter to protein-coding)
  if (keep_protein_coding) {
    gene_info <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters = "biotype",
      values = "protein_coding",
      mart = ensembl
    )
    message("✅ Retrieved protein-coding gene list from Ensembl.")
  } else {
    gene_info <- getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      mart = ensembl
    )
    message("ℹ️ Retrieved all gene types from Ensembl (no filtering).")
  }
  
  # 🧹 Clean: remove genes without HGNC symbols
  gene_info <- gene_info[gene_info$hgnc_symbol != "", ]
  gene_info <- unique(gene_info)
  
  # 📄 Load expression data
  mydata <- read.csv(input_path, sep = '\t')
  gene_col_name <- colnames(mydata)[1]
  
  if (grepl("^ENSG", mydata[[gene_col_name]][1])) {
    # 🧽 Clean Ensembl IDs (remove version numbers)
    mydata[[gene_col_name]] <- sub("\\..*", "", mydata[[gene_col_name]])
    mydata <- mydata[!duplicated(mydata[[gene_col_name]]), ]
    
    # 🔧 Set Ensembl IDs as rownames temporarily
    rownames(mydata) <- mydata[[gene_col_name]]
    mydata[[gene_col_name]] <- NULL
    
    message("✅ Ensembl gene IDs cleaned and duplicates removed.")
  } else {
    stop("❌ Gene IDs do not appear to be in the first column or not in ENSG format.")
  }
  
  # 🔗 Map Ensembl to HGNC Symbols
  mydata$ensembl_gene_id <- rownames(mydata)
  mydata <- merge(mydata, gene_info[, c("ensembl_gene_id", "hgnc_symbol")], by = "ensembl_gene_id")
  
  if (!"hgnc_symbol" %in% colnames(mydata)) {
    stop("❌ Merge failed — HGNC symbol column missing.")
  }
  
  # 🧼 Filter out missing or duplicate symbols
  mydata <- mydata[mydata$hgnc_symbol != "", ]
  mydata <- mydata[!duplicated(mydata$hgnc_symbol), ]
  
  # 📊 Set HGNC symbols as rownames
  rownames(mydata) <- mydata$hgnc_symbol
  mydata$ensembl_gene_id <- NULL
  mydata$hgnc_symbol <- NULL
  
  message(paste0("✅ Ensembl IDs replaced with HGNC gene symbols. Final gene count: ", nrow(mydata)))
  
  # 💾 Prepare for DB Save
  mydata$GeneSymbol <- rownames(mydata)
  mydata <- mydata[, c(ncol(mydata), 1:(ncol(mydata)-1))]  # Move GeneSymbol to first column
  
  # 💾 Save to SQLite
  con <- dbConnect(SQLite(), db_path)
  dbWriteTable(con, table_name, mydata, overwrite = TRUE)
  dbDisconnect(con)
  
  cat(paste0("✅ Expression data saved to table '", table_name, "' in ", db_path, "\n"))
}



# ---- Function 2: build_gene_network ----
build_gene_network <- function(
    db_path = "./Data/BRCA_GeneExpression.db",
    table_expr = "Gene_Expression",
    mi_threshold = 0.8,
    edge_table_prefix = "Gene_Network_Edges"
) {
  # 📦 Required Packages
  packages <- c("igraph", "RSQLite", "DBI", "Matrix", "HiClimR")
  to_install <- setdiff(packages, rownames(installed.packages()))
  if (length(to_install) > 0) install.packages(to_install)
  lapply(packages, require, character.only = TRUE)
  
  # 🔄 Load expression data
  con <- dbConnect(SQLite(), db_path)
  mydata <- dbReadTable(con, table_expr)
  dbDisconnect(con)
  
  # ✅ Safety check for expected column
  if (!"GeneSymbol" %in% colnames(mydata)) {
    stop("❌ Expected column 'GeneSymbol' not found in the expression table.")
  }
  
  # 🧬 Set rownames and prepare matrix (samples x genes)
  rownames(mydata) <- mydata$GeneSymbol
  data_mat <- t(as.matrix(mydata[, setdiff(names(mydata), "GeneSymbol")]))  # Drop GeneSymbol column
  
  # ⚡ Fast correlation
  cor_mat <- fastCor(data_mat, nSplit = 10, upperTri = FALSE, verbose = TRUE)
  cor_mat[is.na(cor_mat)] <- 0
  cor_mat <- pmin(pmax(cor_mat, -0.9999), 0.9999)
  
  # 🔁 Convert correlation to mutual information
  gene_names <- colnames(data_mat)
  rownames(cor_mat) <- gene_names
  colnames(cor_mat) <- gene_names
  
  mi_mat <- -0.5 * log(1 - cor_mat^2)
  diag(mi_mat) <- 0
  rownames(mi_mat) <- gene_names
  colnames(mi_mat) <- gene_names
  
  # 🔎 Extract gene pairs above MI threshold
  edge_indices <- which(mi_mat >= mi_threshold, arr.ind = TRUE)
  edge_indices <- edge_indices[edge_indices[, 1] < edge_indices[, 2], , drop = FALSE]
  
  if (nrow(edge_indices) == 0) {
    stop(paste0("❌ No gene-gene edges passed the MI threshold (", mi_threshold, "). Try lowering it."))
  }
  
  # 📊 Build edge list
  edges <- data.frame(
    from = rownames(mi_mat)[edge_indices[, 1]],
    to = colnames(mi_mat)[edge_indices[, 2]],
    weight = mi_mat[cbind(edge_indices[, 1], edge_indices[, 2])]
  )
  
  # 🧠 Create graph
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # 💾 Save edge list to SQLite DB
  con <- dbConnect(SQLite(), db_path)
  dbWriteTable(con, edge_table_prefix, edges, overwrite = TRUE)
  dbDisconnect(con)
  
  # ✅ Output summary
  cat("✅ Network built with", vcount(g), "nodes and", ecount(g), "edges.\n")
  cat("📁 Edge list saved to table: '", edge_table_prefix, "'\n")
  
  return(invisible(g))
}



# ---- Function 3: extract_node_features -----
extract_node_features <- function(
    db_path = "./Data/BRCA_GeneExpression.db",
    edge_table_prefix = "Gene_Network_Edges",
    expression_table = "Gene_Expression",
    output_table = "Gene_AllFeatures"
) {
  # 📦 Load Required Libraries
  packages <- c("igraph", "DBI", "RSQLite")
  installed <- rownames(installed.packages())
  to_install <- setdiff(packages, installed)
  if (length(to_install) > 0) install.packages(to_install)
  lapply(packages, require, character.only = TRUE)
  
  # 📥 Step 1: Load Network Edge List from SQLite
  con <- dbConnect(SQLite(), db_path)
  edges <- dbReadTable(con, paste0(edge_table_prefix))
  dbDisconnect(con)
  
  # 🔗 Step 2: Build Graph from Edge List
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # 🧱 Optional: Use Only the Largest Connected Component
  comp <- components(g)
  giant <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
  g <- giant
  
  # 📊 Step 3: Compute Topological Features
  pr_scores  <- page_rank(g)$vector
  btw_scores <- betweenness(g, normalized = TRUE)
  cls_scores <- closeness(g, normalized = TRUE)
  eig_scores <- eigen_centrality(g, weights = E(g)$weight)$vector
  deg_scores <- degree(g)
  str_scores <- strength(g, weights = E(g)$weight)
  
  entropy_scores <- sapply(V(g), function(v) {
    incident_edges <- incident(g, v, mode = "all")
    weights <- E(g)[incident_edges]$weight
    p <- weights / sum(weights)
    -sum(p * log2(p + 1e-12))  # Avoid log(0)
  })
  
  # 🧬 Step 4: Load Gene Expression Data for Biological Features
  con <- dbConnect(SQLite(), db_path)
  expr_data <- dbReadTable(con, expression_table)
  dbDisconnect(con)
  
  rownames(expr_data) <- expr_data$GeneSymbol
  expr_mat <- as.matrix(expr_data[, -1])
  genes_in_graph <- names(pr_scores)
  expr_mat <- expr_mat[rownames(expr_mat) %in% genes_in_graph, ]
  
  expr_mean <- rowMeans(expr_mat)
  expr_sd   <- apply(expr_mat, 1, sd)
  
  # 🧾 Step 5: Combine All Features into a Data Frame
  node_features <- data.frame(
    Gene        = genes_in_graph,
    PageRank    = pr_scores[genes_in_graph],
    Betweenness = btw_scores[genes_in_graph],
    Closeness   = cls_scores[genes_in_graph],
    Eigenvector = eig_scores[genes_in_graph],
    Degree      = deg_scores[genes_in_graph],
    Strength    = str_scores[genes_in_graph],
    Entropy     = entropy_scores[genes_in_graph],
    ExprMean    = expr_mean[genes_in_graph],
    ExprSD      = expr_sd[genes_in_graph]
  )
  
  # 💾 Step 6: Save Features to SQLite
  con <- dbConnect(SQLite(), db_path)
  dbWriteTable(con, paste0(output_table), node_features, overwrite = TRUE)
  dbDisconnect(con)
  
  cat("✅ All features extracted and saved to", paste0(output_table), "\n")
}




# ---- Function 4: run_spnfsr_cv
run_spnfsr_cv <- function(
    db_path = "./Data/BRCA_GeneExpression.db",
    features_table_prefix = "Gene_AllFeatures",
    output_table_prefix = "BRCA_ranked_genes_SPNFSR_CV",
    result_table_prefix = "SPNFSR_CV_Results",
    N = 100
) {
  # 📦 Load Libraries
  packages <- c("DBI", "RSQLite", "Matrix", "cluster", "stats", "factoextra")
  to_install <- setdiff(packages, rownames(installed.packages()))
  if (length(to_install) > 0) install.packages(to_install)
  lapply(packages, require, character.only = TRUE)
  
  # 🧬 Load Feature Matrix
  con <- dbConnect(SQLite(), db_path)
  node_features <- dbReadTable(con, paste0(features_table_prefix))
  dbDisconnect(con)
  
  features <- node_features[, c(
    "PageRank", "Betweenness", "Closeness", "Eigenvector",
    "Degree", "Strength", "Entropy", "ExprMean", "ExprSD"
  )]
  rownames(features) <- node_features$Gene
  X_full <- scale(as.matrix(features))
  
  # 🎯 Define Parameter Grid
  alpha_list <- c(0.01, 0.1, 0.5, 1, 5, 10, 25)
  beta_list <- c(0.001, 0.01, 0.05, 0.1, 0.5, 1, 2)
  k_list <- c(3, 5, 7, 9)
  sigma2_list <- c(25, 50, 100, 150, 200, 300, 500)
  
  # 🧰 Helper Functions
  construct_Q <- function(W) diag(as.vector(1 / (2 * sqrt(W^2 + 1e-10))))
  construct_R <- function(W, X) {
    temp <- X %*% W
    diag(1 / (2 * sqrt(as.vector(temp)^2 + 1e-10)))
  }
  compute_similarity <- function(X, k_thresh, sigma2) {
    m <- nrow(X)
    O <- as.matrix(dist(X))
    S <- matrix(0, m, m)
    for (i in 1:m) {
      for (j in 1:m) {
        if (O[i, j] <= k_thresh) {
          S[i, j] <- exp(-O[i, j]^2 / sigma2)
        }
      }
    }
    S
  }
  evaluate_silhouette <- function(X_proj, k = 3) {
    if (any(is.na(X_proj)) || any(is.infinite(X_proj))) return(NA)
    if (nrow(X_proj) < k || length(unique(rowSums(X_proj))) < k) return(NA)
    km <- tryCatch(kmeans(X_proj, centers = k, nstart = 10), error = function(e) NULL)
    if (!is.null(km)) {
      sil <- silhouette(km$cluster, dist(X_proj))
      return(mean(sil[, 3]))
    } else {
      return(NA)
    }
  }
  
  # 🔁 Grid Search Loop
  results <- list()
  param_id <- 1
  
  for (alpha in alpha_list) {
    for (beta in beta_list) {
      for (k_thresh in k_list) {
        for (sigma2 in sigma2_list) {
          cat("🔧 alpha =", alpha, " beta =", beta, " k =", k_thresh, " sigma2 =", sigma2, "\n")
          
          X <- X_full
          m <- nrow(X)
          n <- ncol(X)
          
          S <- compute_similarity(X, k_thresh, sigma2)
          I_m <- diag(m)
          L <- I_m - S - t(S) + S %*% t(S)
          
          M <- t(X) %*% L %*% X
          M1 <- (abs(M) + M) / 2
          M2 <- (abs(M) - M) / 2
          
          W <- matrix(runif(n), n, 1)
          Q <- diag(n)
          R <- diag(m)
          
          for (iter in 1:N) {
            XTRX <- t(X) %*% (R %*% X)
            a <- alpha * M2 %*% W + XTRX %*% W
            b <- (XTRX + beta * Q + alpha * M1) %*% W
            W <- W * (a / (b + 1e-10))
            Q <- construct_Q(W)
            R <- construct_R(W, X)
          }
          
          # Evaluate projection
          X_proj <- X %*% W
          sil_score <- evaluate_silhouette(X_proj)
          
          results[[param_id]] <- list(
            alpha = alpha,
            beta = beta,
            k_thresh = k_thresh,
            sigma2 = sigma2,
            silhouette = sil_score,
            W = W
          )
          param_id <- param_id + 1
        }
      }
    }
  }
  
  # 🎯 Filter valid results
  silhouettes <- sapply(results, function(r) r$silhouette)
  valid_idxs <- which(!is.na(silhouettes))
  
  if (length(valid_idxs) == 0) {
    stop("❌ All parameter combinations failed to produce valid silhouette scores.")
  }
  
  best_idx <- valid_idxs[which.max(silhouettes[valid_idxs])]
  best_result <- results[[best_idx]]
  
  # ✅ Final Gene Ranking
  W_opt <- best_result$W
  W_norm <- W_opt / max(abs(W_opt))
  W_final <- t(W_norm)^2
  
  X_proj_final <- X_full %*% t(W_final)
  gene_scores <- rowSums(X_proj_final^2)
  
  ranked_genes <- data.frame(
    Gene = rownames(X_full),
    SPNFSR_Score = gene_scores
  )
  ranked_genes <- ranked_genes[order(-ranked_genes$SPNFSR_Score), ]
  
  # 💾 Save Ranked Genes
  con <- dbConnect(SQLite(), db_path)
  dbWriteTable(con, paste0(output_table_prefix), ranked_genes, overwrite = TRUE)
  dbDisconnect(con)
  
  cat("✅ Cross-validated SPNFSR completed.\n")
  cat(paste0("Best parameters: alpha = ", best_result$alpha,
             ", beta = ", best_result$beta,
             ", k = ", best_result$k_thresh,
             ", sigma2 = ", best_result$sigma2, "\n"))
  
  # 💾 Save All Grid Search Results
  cv_results_df <- do.call(rbind, lapply(results, function(res) {
    data.frame(
      alpha = res$alpha,
      beta = res$beta,
      k_thresh = res$k_thresh,
      sigma2 = res$sigma2,
      silhouette = res$silhouette
    )
  }))
  
  con <- dbConnect(SQLite(), db_path)
  dbWriteTable(con, paste0(result_table_prefix), cv_results_df, overwrite = TRUE)
  dbDisconnect(con)
  
  cat("✅ All SPNFSR cross-validation results saved to DB.\n")
}


# ---- Function 5: plot_spnfsr_results -----
plot_spnfsr_results <- function(
    db_path = "./Data/BRCA_GeneExpression.db",
    result_table_prefix = "SPNFSR_CV_Results",
    mi_threshold = 0.8,
    output_pdf = "SPNFSR_Tuning_Results.pdf",
    scale_scores = TRUE,
    top_n = 15,
    show_plot = TRUE,
    pdf_width = 12,
    pdf_height = 10,
    fill_colors = list(
      top15 = "#e69f00",
      alpha = "#56B4E9",
      beta = "#009E73",
      k_thresh = "#CC79A7",
      sigma2 = "#F0E442"
    )
) {
  # 📦 Load Required Libraries
  packages <- c("DBI", "RSQLite", "dplyr", "ggplot2", "patchwork", "scales")
  installed <- rownames(installed.packages())
  to_install <- setdiff(packages, installed)
  if (length(to_install) > 0) install.packages(to_install)
  lapply(packages, require, character.only = TRUE)
  
  # 📥 Load CV Results Table
  con <- dbConnect(SQLite(), db_path)
  node_features <- dbReadTable(con, paste0(result_table_prefix))
  dbDisconnect(con)
  
  # 🔢 Convert silhouette to numeric and clean invalid
  node_features$silhouette <- as.numeric(node_features$silhouette)
  invalid_rows <- sum(!is.finite(node_features$silhouette))
  if (invalid_rows > 0) {
    message("⚠️ Removed ", invalid_rows, " rows with invalid silhouette scores.")
  }
  node_features <- node_features %>%
    filter(is.finite(silhouette))
  
  # 📈 Scale if required
  if (scale_scores) {
    node_features <- node_features %>%
      mutate(silhouette_scaled = rescale(silhouette))
    yval <- "silhouette_scaled"
    ylabel <- "Scaled Silhouette Score (0–1)"
  } else {
    node_features <- node_features %>%
      mutate(silhouette_scaled = silhouette)
    yval <- "silhouette"
    ylabel <- "Silhouette Score"
  }
  
  # 🧼 Define Clean Plot Theme
  theme_clean <- theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.caption = element_text(hjust = 1, face = "italic", size = 9)
    )
  
  # 📊 Top N Configurations by Silhouette
  top_configs <- node_features %>%
    arrange(desc(!!sym(yval))) %>%
    head(top_n)
  
  p1 <- ggplot(top_configs, aes(
    x = reorder(paste0("α=", alpha, ", β=", beta, ", k=", k_thresh, ", σ²=", sigma2), !!sym(yval)),
    y = !!sym(yval)
  )) +
    geom_bar(stat = "identity", fill = fill_colors$top15) +
    coord_flip() +
    labs(title = paste0("Top ", top_n, " Configurations"),
         x = NULL, y = ylabel) +
    theme_clean
  
  # 📈 Parameter Effect Plots
  p2 <- ggplot(node_features, aes(x = as.factor(alpha), y = !!sym(yval))) +
    stat_summary(fun = mean, geom = "bar", fill = fill_colors$alpha) +
    labs(title = "Alpha vs Score", x = "α", y = NULL) +
    theme_clean
  
  p3 <- ggplot(node_features, aes(x = as.factor(beta), y = !!sym(yval))) +
    stat_summary(fun = mean, geom = "bar", fill = fill_colors$beta) +
    labs(title = "Beta vs Score", x = "β", y = NULL) +
    theme_clean
  
  p4 <- ggplot(node_features, aes(x = as.factor(k_thresh), y = !!sym(yval))) +
    stat_summary(fun = mean, geom = "bar", fill = fill_colors$k_thresh) +
    labs(title = "k_thresh vs Score", x = "k", y = NULL) +
    theme_clean
  
  if ("sigma2" %in% colnames(node_features)) {
    p5 <- ggplot(node_features, aes(x = as.factor(sigma2), y = !!sym(yval))) +
      stat_summary(fun = mean, geom = "bar", fill = fill_colors$sigma2) +
      labs(title = "Sigma² vs Score", x = "σ²", y = NULL) +
      theme_clean
  } else {
    p5 <- ggplot() + theme_void() + labs(title = "Sigma² Not Available")
  }
  
  # 🖼️ Combine All Plots
  final_plot <- (p1 / (p2 | p3 | p4 | p5)) +
    plot_annotation(title = paste0("SPNFSR Parameter Tuning (MI ≥ ", mi_threshold, ")"))
  
  # 💾 Save to PDF
  ggsave(output_pdf, plot = final_plot, width = pdf_width, height = pdf_height, units = "in")
  
  # 📤 Show in Viewer (optional)
  if (show_plot) {
    print(final_plot)
  }
  
  cat("✅ Plot saved to", output_pdf, "\n")
}

# ---- Function 6: summarize_top_gene_neighbors ----
summarize_top_gene_neighbors <- function(
    db_path = "./Data/BRCA_GeneExpression.db",
    ranked_table = "BRCA_ranked_genes_SPNFSR_CV",
    edges_table = "Gene_Network_Edges",
    expression_table = "Gene_Expression",
    output_table = "TopGene_Neighbor_Summary",
    top_n = 10
) {
  # 📦 Load Required Packages
  packages <- c("DBI", "RSQLite", "igraph", "dplyr")
  installed <- rownames(installed.packages())
  to_install <- setdiff(packages, installed)
  if (length(to_install) > 0) install.packages(to_install)
  lapply(packages, require, character.only = TRUE)
  
  # 🔗 Connect to DB and Load Data
  con <- dbConnect(SQLite(), db_path)
  top_genes <- dbReadTable(con, ranked_table) %>% head(top_n)
  edges <- dbReadTable(con, edges_table)
  expr_data <- dbReadTable(con, expression_table)
  dbDisconnect(con)
  
  # 🧠 Build Network Graph
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # 🧬 Calculate Mean Expression (TPM)
  expr_data <- expr_data %>%
    mutate(ExprMean = rowMeans(across(where(is.numeric)))) %>%
    select(GeneSymbol, ExprMean)
  
  # 🧾 Summarize Neighbors for Top Genes
  summary_df <- top_genes %>%
    rename(GeneSymbol = Gene) %>%
    left_join(expr_data, by = "GeneSymbol") %>%
    rowwise() %>%
    mutate(
      Neighbors = list(neighbors(g, GeneSymbol) %>% names()),
      NumNeighbors = length(Neighbors),
      NeighborNames = paste(Neighbors, collapse = ", ")
    ) %>%
    ungroup() %>%
    select(GeneSymbol, SPNFSR_Score, ExprMean, NumNeighbors, NeighborNames)
  
  # 💾 Save Result to SQLite DB
  con <- dbConnect(SQLite(), db_path)
  dbWriteTable(con, output_table, summary_df, overwrite = TRUE)
  dbDisconnect(con)
  
  cat("✅ Neighbor summary saved to table:", output_table, "\n")
  
  return(invisible(summary_df))
}

# ---- Function 7: check_network_connectivity ----
check_network_connectivity <- function(
    mi_threshold,
    db_path = NULL,
    edge_table_prefix = "Gene_Network_Edges"
) {
  # 📦 Load required packages
  packages <- c("igraph", "DBI", "RSQLite")
  to_install <- setdiff(packages, rownames(installed.packages()))
  if (length(to_install) > 0) install.packages(to_install)
  lapply(packages, require, character.only = TRUE)
  
  # 🗂️ Set DB path if not provided
  if (is.null(db_path)) {
    db_path <- paste0("./Data/BRCA_GeneExpression_", mi_threshold, ".db")
  }
  
  # 🧬 Load edge list
  con <- dbConnect(SQLite(), db_path)
  table_name <- paste0(edge_table_prefix)  # e.g., "Gene_Network_Edges"
  edges <- dbReadTable(con, table_name)
  dbDisconnect(con)
  
  if (nrow(edges) == 0) {
    stop("❌ Edge list is empty — cannot assess connectivity.")
  }
  
  # 🔗 Build graph
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # 🧩 Check connected components
  comps <- components(g)
  num_components <- comps$no
  largest_size <- max(comps$csize)
  percent_giant <- round(100 * largest_size / vcount(g), 2)
  
  # 🖨️ Output summary
  cat("📊 Network Connectivity Summary (MI Threshold =", mi_threshold, ")\n")
  cat(" - Total nodes:              ", vcount(g), "\n")
  cat(" - Total edges:              ", ecount(g), "\n")
  cat(" - Number of components:     ", num_components, "\n")
  cat(" - Size of largest component:", largest_size, "\n")
  cat(" - % of nodes in largest:    ", percent_giant, "%\n")
  
  # Return for further use if needed
  return(invisible(list(
    graph = g,
    components = comps,
    threshold = mi_threshold,
    total_nodes = vcount(g),
    total_edges = ecount(g),
    num_components = num_components,
    largest_component_size = largest_size,
    percent_in_giant = percent_giant
  )))
}

