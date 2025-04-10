# 📦 Load required packages (install if missing)
packages <- c("igraph", "RSQLite", "DBI", "Matrix", "HiClimR")
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if (length(to_install) > 0) install.packages(to_install)
library(igraph)
library(DBI)
library(RSQLite)
library(Matrix)
library(HiClimR)  # New package to speed up correlation

# ---- Load gene expression data from SQLite ----
con <- dbConnect(SQLite(), "./Data/BRCA_GeneExpression.db")
mydata <- dbReadTable(con, "Gene_Expression")
dbDisconnect(con)

# ---- Reconstruct the matrix ----
rownames(mydata) <- mydata$GeneID
data_mat <- t(as.matrix(mydata[, -1]))  # Transpose: samples = rows

# ---- Compute Pearson correlation FAST ----
# Use fastCor for performance (uses C++ backend)
cor_mat <- fastCor(data_mat, nSplit = 10, upperTri = FALSE, verbose = TRUE)
cor_mat[is.na(cor_mat)] <- 0

# ---- Convert to mutual information approximation ----
cor_mat <- pmin(pmax(cor_mat, -0.9999), 0.9999)  # Prevent log(0)
mi_mat <- -0.5 * log(1 - cor_mat^2)
diag(mi_mat) <- 0

# ---- Filter by MI threshold ----
mi_threshold <- 0.7
edge_indices <- which(mi_mat >= mi_threshold, arr.ind = TRUE)
edge_indices <- edge_indices[edge_indices[,1] < edge_indices[,2], ]

# ---- Build edges data frame ----
edges <- data.frame(
  from = rownames(mi_mat)[edge_indices[,1]],
  to = colnames(mi_mat)[edge_indices[,2]],
  weight = mi_mat[edge_indices]
)

# ---- Build graph ----
g <- graph_from_data_frame(edges, directed = FALSE)

# ---- Save edge list to text ----
write.table(edges, file = "./Data/BRCA_network_edges.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# ---- Save edge list to SQLite ----
con <- dbConnect(SQLite(), "./Data/BRCA_GeneExpression.db")
dbWriteTable(con, "Gene_Network_Edges", edges, overwrite = TRUE)
dbDisconnect(con)

cat("✅ Network built with", vcount(g), "nodes and", ecount(g), "edges.\n")
cat("📁 Edge list saved to text and SQLite database.\n")
