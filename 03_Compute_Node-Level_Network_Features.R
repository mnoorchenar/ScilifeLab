# ---- Step 3: Load Network from SQLite and Compute Node-Level Features ----

library(igraph)
library(DBI)
library(RSQLite)

# Load edge list from database
con <- dbConnect(SQLite(), "./Data/BRCA_GeneExpression.db")
edges <- dbReadTable(con, "Gene_Network_Edges")
dbDisconnect(con)

# Build graph from edges
g <- graph_from_data_frame(edges, directed = FALSE)

# ---- Compute Features ----

# PageRank
pr_scores <- page_rank(g)$vector

# Betweenness Centrality
btw_scores <- betweenness(g, normalized = TRUE)

# Closeness Centrality
cls_scores <- closeness(g, normalized = TRUE)

# ---- Compute entropy per node based on connected edge weights ----
entropy_scores <- sapply(V(g), function(v) {
  incident_edges <- incident(g, v, mode = "all")
  weights <- E(g)[incident_edges]$weight
  p <- weights / sum(weights)
  -sum(p * log2(p + 1e-12))  # Add small value to avoid log(0)
})


# ---- Combine features into a data frame ----
node_features <- data.frame(
  Gene = names(pr_scores),
  PageRank = pr_scores,
  Betweenness = btw_scores,
  Closeness = cls_scores,
  Entropy = entropy_scores
)

# ---- Save results ----
write.table(node_features, file = "./Data/BRCA_node_features.txt", sep = "\t", row.names = FALSE, quote = FALSE)

con <- dbConnect(SQLite(), "./Data/BRCA_GeneExpression.db")
dbWriteTable(con, "Gene_Network_Features", node_features, overwrite = TRUE)
dbDisconnect(con)

cat("✅ Node features computed from loaded network and saved.\n")
