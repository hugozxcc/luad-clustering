############################################
# STEP 6: Expression Heatmap of DE Genes
# NOTE ON INTERPRETATION (important):
# The genes shown here are selected BECAUSE they are differentially expressed
# between cluster 1 and cluster 2. Displaying them grouped by cluster will
# therefore always show visual separation -- this figure is ILLUSTRATIVE of the
# DE signature, NOT independent validation that the clusters are distinct.
# The caption/title states this explicitly.
############################################

library(pheatmap)
set.seed(42)

logcpm <- readRDS("02_data_processed/TCGA_LUAD_norm_logCPM.rds")
cl <- read.csv("02_data_processed/TCGA_LUAD_cluster_assignment.csv")
cl <- cl[match(colnames(logcpm), cl$sample_id), ]

# ---- Select genes from the SIGNIFICANCE-FILTERED lists (adj.P < 0.05, |logFC|>1)
# Earlier versions sorted the full DE table by raw logFC with no significance
# filter, which could admit high-FC but non-significant genes.
sig_up   <- read.csv("04_results/de/DE_sig_up_cluster2.csv")
sig_down <- read.csv("04_results/de/DE_sig_down_cluster2.csv")

top_up   <- head(sig_up[order(-sig_up$logFC), "gene"], 25)
top_down <- head(sig_down[order(sig_down$logFC), "gene"], 25)
genes <- unique(c(top_up, top_down))
genes <- genes[genes %in% rownames(logcpm)]
cat("Heatmap genes:", length(genes), "(", length(top_up), "up +", length(top_down), "down )\n")

mat <- logcpm[genes, ]

# z-score per gene, clip to +/-2
mat_scaled <- t(scale(t(mat)))
mat_scaled[mat_scaled >  2] <-  2
mat_scaled[mat_scaled < -2] <- -2

# ---- Balanced sample selection: equal random draw from BOTH clusters ----
# Clusters are now 349 vs 191, so we sample the same number from each rather
# than "all of one + a subset of the other" (which only made sense when one
# cluster had 13 samples).
idx_c1 <- which(cl$cluster == 1)
idx_c2 <- which(cl$cluster == 2)
n_show <- min(80, length(idx_c1), length(idx_c2))   # equal per cluster
idx_show <- c(sample(idx_c1, n_show), sample(idx_c2, n_show))

mat_show <- mat_scaled[, idx_show]
cl_show  <- cl[idx_show, ]

annotation_col <- data.frame(Cluster = factor(cl_show$cluster))
rownames(annotation_col) <- cl_show$sample_id

png("04_results/heatmap_top50_cluster2_readable.png", width = 1400, height = 900)
pheatmap(
  mat_show,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  cluster_cols = TRUE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = paste0("DE-selected genes (top 25 up + 25 down), ", n_show,
                " samples/cluster\n",
                "Illustrative of the cluster DE signature, not independent validation")
)
dev.off()
cat("Heatmap saved.\n")