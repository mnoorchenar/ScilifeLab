# ---- Step 4: Rank Genes Using Unsupervised Learning + Visualizations ----
# Author: [Your Name]
# Description: Rank genes from a BRCA network using PCA and Composite Scoring, with visualizations.

# Load libraries
library(stats)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# ---- Load Node Features ----
node_features <- read.table("./Data/BRCA_node_features.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# ---- Prepare Feature Matrix ----
features <- node_features[, c("PageRank", "Betweenness", "Closeness", "Entropy")]
rownames(features) <- node_features$Gene

# ---- Scale Features ----
features_scaled <- scale(features)

# ---- PCA-Based Ranking ----
pca <- prcomp(features_scaled)
pca_scores <- abs(pca$x[, 1])

pca_ranked_genes <- data.frame(
  Gene = rownames(features_scaled),
  PCA_Score = pca_scores
)
pca_ranked_genes <- pca_ranked_genes[order(-pca_ranked_genes$PCA_Score), ]

# ---- Composite Score Ranking ----
composite_score <- rowMeans(features_scaled)

composite_ranked_genes <- data.frame(
  Gene = rownames(features_scaled),
  Composite_Score = composite_score
)
composite_ranked_genes <- composite_ranked_genes[order(-composite_ranked_genes$Composite_Score), ]

# ---- Save Results ----
write.table(pca_ranked_genes, file = "./Data/BRCA_ranked_genes_PCA.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(composite_ranked_genes, file = "./Data/BRCA_ranked_genes_Composite.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cat("✅ Gene ranking completed and saved.\n")

# ---- 📊 Visualization Section ----

# 1. Scree Plot (Variance Explained by Each Principal Component)
scree_data <- data.frame(PC = paste0("PC", 1:length(pca$sdev)),
                         Variance = (pca$sdev)^2 / sum(pca$sdev^2))

ggplot(scree_data, aes(x = PC, y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Scree Plot - Variance Explained by PCs", y = "Proportion of Variance", x = "Principal Component")

# 2. PCA Biplot (Top 2 PCs)
pca_df <- as.data.frame(pca$x)
pca_df$Gene <- rownames(pca_df)

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "PCA Biplot of Genes", x = "PC1", y = "PC2")

# 3. Top 20 Genes by PCA Score (Barplot)
top20_pca <- head(pca_ranked_genes, 20)
ggplot(top20_pca, aes(x = reorder(Gene, PCA_Score), y = PCA_Score)) +
  geom_bar(stat = "identity", fill = "tomato") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 Genes by PCA Score", x = "Gene", y = "PCA Score")

# 4. Top 20 Genes by Composite Score (Barplot)
top20_composite <- head(composite_ranked_genes, 20)
ggplot(top20_composite, aes(x = reorder(Gene, Composite_Score), y = Composite_Score)) +
  geom_bar(stat = "identity", fill = "darkcyan") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 Genes by Composite Score", x = "Gene", y = "Composite Score")

# 5. Heatmap of Feature Correlation
cor_matrix <- cor(features_scaled)
pheatmap(cor_matrix,
         main = "Correlation Heatmap of Node-Level Features",
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         display_numbers = TRUE,
         fontsize_number = 10)

cat("📈 Visualizations generated successfully.\n")
