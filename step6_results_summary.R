############################################
# STEP 6: Results Summary Tables
# Exports top DE genes (significance-filtered) for the results section.
############################################

de_full  <- read.csv("04_results/de/DE_limma_cluster2_vs_cluster1_full.csv")

# Use significance-filtered genes (adj.P.Val < 0.05 & |logFC| > 1) so the
# "top" tables cannot include high-fold-change but non-significant genes.
sig_up   <- subset(de_full, adj.P.Val < 0.05 & logFC >  1)
sig_down <- subset(de_full, adj.P.Val < 0.05 & logFC < -1)

# Top 20 UP in cluster 2 (largest positive logFC among significant genes)
top_up <- head(sig_up[order(-sig_up$logFC), c("gene", "logFC", "adj.P.Val")], 20)
write.csv(top_up, "04_results/de/top20_up_cluster2.csv", row.names = FALSE)

# Top 20 DOWN in cluster 2 (most negative logFC among significant genes)
top_down <- head(sig_down[order(sig_down$logFC), c("gene", "logFC", "adj.P.Val")], 20)
write.csv(top_down, "04_results/de/top20_down_cluster2.csv", row.names = FALSE)

cat("Saved top20 up/down DE tables.\n")
cat("Significant UP:", nrow(sig_up), " DOWN:", nrow(sig_down), "\n")