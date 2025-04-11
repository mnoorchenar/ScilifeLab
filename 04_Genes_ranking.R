# ---- LS Score Gene Ranking for BRCA Network ----
# Description: Rank genes using SPNFSR-inspired LS Score and generate visualizations

# 📦 Load Libraries
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# 📂 Load Node Features
node_features <- read.table("./Data/BRCA_node_features.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# 🧮 Prepare Feature Matrix (Topological + Informative)
features <- node_features[, c("PageRank", "Betweenness", "Closeness", "Entropy")]
rownames(features) <- node_features$Gene

# ⚖️ Scale Features
features_scaled <- scale(features)

# ---- 🧠 LS Score Calculation ----
scale_values <- attr(features_scaled, "scaled:scale")
feature_variance <- scale_values^2
feature_weights <- feature_variance / sum(feature_variance)

print(feature_weights)


ls_scores <- as.matrix(features_scaled) %*% feature_weights

ls_ranked_genes <- data.frame(
  Gene = rownames(features_scaled),
  LS_Score = ls_scores
)
ls_ranked_genes <- ls_ranked_genes[order(-ls_ranked_genes$LS_Score), ]


# 💾 Save Ranked List
write.table(ls_ranked_genes, file = "./Data/BRCA_ranked_genes_LS.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cat("✅ LS Score-based gene ranking completed and saved.\n")

# ---- 📈 Visualization ----

# 1. Top 20 Genes by LS Score (Barplot)
top20_ls <- head(ls_ranked_genes, 20)
plot_ls <- ggplot(top20_ls, aes(x = reorder(Gene, LS_Score), y = LS_Score)) +
  geom_bar(stat = "identity", fill = "mediumseagreen") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 Genes by LS Score", x = "Gene", y = "LS Score")

# 2. Heatmap of Feature Correlations
cor_matrix <- cor(features_scaled)
plot_heatmap <- pheatmap(
  cor_matrix,
  main = "Correlation Heatmap of Node-Level Features",
  color = colorRampPalette(brewer.pal(9, "Blues"))(100),
  display_numbers = TRUE,
  fontsize_number = 10
)

# ---- 🖼 Save Visualizations ----
pdf(file = "./Results/BRCA_LS_score_visualizations.pdf", width = 8, height = 6)

print(plot_ls)

# Note: pheatmap is already plotted on generation
pheatmap(
  cor_matrix,
  main = "Correlation Heatmap of Node-Level Features",
  color = colorRampPalette(brewer.pal(9, "Blues"))(100),
  display_numbers = TRUE,
  fontsize_number = 10
)

dev.off()

cat("📄 Visualizations saved to './Results/BRCA_LS_score_visualizations.pdf'\n")
