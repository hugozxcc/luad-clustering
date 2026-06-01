library(pheatmap)

set.seed(42)

logcpm <- readRDS("02_data_processed/TCGA_LUAD_norm_logCPM.rds")
cl <- read.csv("02_data_processed/TCGA_LUAD_cluster_assignment.csv")
cl <- cl[match(colnames(logcpm), cl$sample_id), ]

de <- read.csv("04_results/de/DE_limma_cluster2_vs_cluster1_full.csv")

# top 50 UP genes that exist in logcpm
candidate <- de[order(-de$logFC), "gene"]
candidate_ok <- candidate[candidate %in% rownames(logcpm)]
top50 <- unique(candidate_ok)[1:50]

mat <- logcpm[top50, ]

# z-score per gene
mat_scaled <- t(scale(t(mat)))
mat_scaled[mat_scaled > 2] <- 2
mat_scaled[mat_scaled < -2] <- -2

# ---- Select samples: all cluster2 + random subset of cluster1 ----
idx_c2 <- which(cl$cluster == 2)
idx_c1 <- which(cl$cluster == 1)

n_c1_show <- 60  # bisa 50-100
idx_c1_sub <- sample(idx_c1, n_c1_show)

idx_show <- c(idx_c1_sub, idx_c2)  # cluster1 subset first, then all cluster2

mat_show <- mat_scaled[, idx_show]
cl_show <- cl[idx_show, ]

annotation_col <- data.frame(Cluster = factor(cl_show$cluster))
rownames(annotation_col) <- cl_show$sample_id

png("04_results/heatmap_top50_cluster2_readable.png", width = 1400, height = 900)
pheatmap(
  mat_show,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  cluster_cols = TRUE,                 # boleh TRUE karena kolom sudah sedikit
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Top 50 Upregulated Genes (All Cluster 2 + 60 samples from Cluster 1)"
)
dev.off()
