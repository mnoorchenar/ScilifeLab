# ---- Function 1: prepare_expression_data ----
prepare_expression_data <- function(
    input_path = "Data/TCGA-BRCA.star_tpm.tsv",
    db_path = "./Data/BRCA_GeneExpression.db",
    table_name = "Gene_Expression",
    keep_protein_coding = TRUE
) {
  # Required Packages 
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
  if (!requireNamespace("RSQLite", quietly = TRUE)) install.packages("RSQLite")
  if (!requireNamespace("DBI", quietly = TRUE)) install.packages("DBI")
  
  library(biomaRt)
  library(DBI)
  library(RSQLite)
  
  # Connect to Ensembl 
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 104)
  
  # Get Protein-Coding Gene List 
  # c("ensembl_gene_id", "hgnc_symbol", "gene_biotype")
  protein_coding_genes <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "biotype",
    values = "protein_coding",
    mart = ensembl
  )
  
  # Load Expression File 
  mydata <- read.csv(input_path, sep = '\t')
  
  gene_col_name <- colnames(mydata)[1]
  if (grepl("^ENSG", mydata[[gene_col_name]][1])) {
    
    # Remove version numbers from Ensembl IDs
    mydata[[gene_col_name]] <- sub("\\..*", "", mydata[[gene_col_name]])
    
    # Remove duplicates
    mydata <- mydata[!duplicated(mydata[[gene_col_name]]), ]
    
    # Set rownames
    rownames(mydata) <- mydata[[gene_col_name]]
    mydata[[gene_col_name]] <- NULL
    
    message("✅ Gene IDs cleaned and duplicates removed.")
  } else {
    stop("❌ Gene IDs do not appear to be in the first column.")
  }
  
  # Filter Protein-Coding Only
  if (keep_protein_coding) {
    mydata <- mydata[rownames(mydata) %in% protein_coding_genes$ensembl_gene_id, ]
    message("✅ Filtered to protein-coding genes.")
  }
  
  # Prepare for DB Save 
  mydata$GeneID <- rownames(mydata)
  mydata <- mydata[, c(ncol(mydata), 1:(ncol(mydata)-1))]
  
  # Save to SQLite 
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
  # 📦 Install and load required packages
  packages <- c("igraph", "RSQLite", "DBI", "Matrix", "HiClimR")
  installed <- rownames(installed.packages())
  to_install <- setdiff(packages, installed)
  if (length(to_install) > 0) install.packages(to_install)
  lapply(packages, require, character.only = TRUE)
  
  # 🔄 Load expression data
  con <- dbConnect(SQLite(), db_path)
  mydata <- dbReadTable(con, table_expr)
  dbDisconnect(con)
  
  # 🧬 Prepare expression matrix
  rownames(mydata) <- mydata$GeneID
  data_mat <- t(as.matrix(mydata[, -1]))  # Genes = columns, Samples = rows
  
  # ⚡ Fast correlation
  cor_mat <- fastCor(data_mat, nSplit = 10, upperTri = FALSE, verbose = TRUE)
  cor_mat[is.na(cor_mat)] <- 0
  cor_mat <- pmin(pmax(cor_mat, -0.9999), 0.9999)
  
  # 🔁 Convert to mutual information
  mi_mat <- -0.5 * log(1 - cor_mat^2)
  diag(mi_mat) <- 0
  
  # 🔎 Filter edges by MI threshold
  edge_indices <- which(mi_mat >= mi_threshold, arr.ind = TRUE)
  edge_indices <- edge_indices[edge_indices[, 1] < edge_indices[, 2], ]
  
  edges <- data.frame(
    from = rownames(mi_mat)[edge_indices[, 1]],
    to = colnames(mi_mat)[edge_indices[, 2]],
    weight = mi_mat[edge_indices]
  )
  
  # 🧠 Build graph (optional but useful if you want to return it)
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # 💾 Save to database
  con <- dbConnect(SQLite(), db_path)
  dbWriteTable(con, paste0(edge_table_prefix), edges, overwrite = TRUE)
  dbDisconnect(con)
  
  # ✅ Console output
  cat("✅ Network built with", vcount(g), "nodes and", ecount(g), "edges.\n")
  cat("📁 Edge list saved to table:", paste0(edge_table_prefix), "\n")
  
  # Optionally return the graph object
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
  
  rownames(expr_data) <- expr_data$GeneID
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



# ---- Function 4: run_spnfsr_cv -----
run_spnfsr_cv <- function(
    db_path = "./Data/BRCA_GeneExpression.db",
    features_table_prefix = "Gene_AllFeatures",
    output_table_prefix = "BRCA_ranked_genes_SPNFSR_CV",
    result_table_prefix = "SPNFSR_CV_Results",
    N = 100
) {
  # 📦 Load Libraries
  packages <- c("DBI", "RSQLite", "Matrix", "cluster", "stats", "factoextra")
  installed <- rownames(installed.packages())
  to_install <- setdiff(packages, installed)
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
  construct_Q <- function(W) {
    diag_vec <- 1 / (2 * sqrt(W^2 + 1e-10))
    diag(as.vector(diag_vec))
  }
  construct_R <- function(W, X) {
    temp <- X %*% W
    diag_vec <- 1 / (2 * sqrt(as.vector(temp)^2 + 1e-10))
    diag(diag_vec)
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
    km <- tryCatch(kmeans(X_proj, centers = k, nstart = 10), error = function(e) NULL)
    if (!is.null(km)) {
      sil <- silhouette(km$cluster, dist(X_proj))
      mean(sil[, 3])
    } else {
      NA
    }
  }
  
  # 🔁 Grid Search Loop
  results <- list()
  param_id <- 1
  
  for (alpha in alpha_list) {
    for (beta in beta_list) {
      for (k_thresh in k_list) {
        for (sigma2 in sigma2_list) {
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
  
  # 🏆 Find Best Configuration
  silhouettes <- sapply(results, function(r) r$silhouette)
  best_idx <- which.max(silhouettes)
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
  
  # 💾 Save Ranked Genes to Database
  con <- dbConnect(SQLite(), db_path)
  dbWriteTable(con, paste0(output_table_prefix), ranked_genes, overwrite = TRUE)
  dbDisconnect(con)
  
  cat("✅ Cross-validated SPNFSR completed.\n")
  cat(paste0("Best parameters: alpha = ", best_result$alpha,
             ", beta = ", best_result$beta,
             ", k = ", best_result$k_thresh,
             ", sigma2 = ", best_result$sigma2, "\n"))
  
  # 📊 Save All Grid Search Results
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


# ---- Function 4: plot_spnfsr_results -----
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
