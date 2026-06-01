############################################
# STEP 5: Differential Expression + GO/KEGG Enrichment + Volcano
# Project: TCGA-LUAD clustering
# Author: Fadil Nugroho
#
# REVISIONS (vs previous version):
#  - limma-trend: eBayes(trend=TRUE, robust=TRUE) for logCPM input
#  - neutral group labels (clusters are now 349 vs 191, not 13-sample minority)
#  - guarded enrichGO/enrichKEGG/dotplot against empty gene sets
#  - merged former step5b (volcano + downregulated enrichment) into one script;
#    removed hardcoded Windows setwd() for portability
############################################

# ========== 0. Setup ==========
rm(list = ls())
setwd(".")

# ========== 1. Load libraries ==========
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "limma", "edgeR", "clusterProfiler",
  "org.Hs.eg.db", "AnnotationDbi", "enrichplot"
), ask = FALSE, update = FALSE)

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")

library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(ggrepel)

# ========== 2. Create output directories ==========
dir.create("04_results", showWarnings = FALSE)
dir.create("04_results/de", showWarnings = FALSE)
dir.create("04_results/go", showWarnings = FALSE)

# ========== 3. Load data ==========
logcpm <- readRDS("02_data_processed/TCGA_LUAD_norm_logCPM.rds")
cat("Expression matrix:", dim(logcpm), "\n")

cl <- read.csv("02_data_processed/TCGA_LUAD_cluster_assignment.csv")
cat("Cluster table:\n")
print(table(cl$cluster))

# Ensure correct order
stopifnot(all(cl$sample_id %in% colnames(logcpm)))
cl <- cl[match(colnames(logcpm), cl$sample_id), ]
stopifnot(all(cl$sample_id == colnames(logcpm)))

# ========== 4. Differential Expression (limma-trend) ==========
# Neutral labels: with the Ward k=2 split the groups are 349 vs 191,
# so "minor/major" no longer applies.
group <- factor(cl$cluster,
                levels = c(1, 2),
                labels = c("Cluster1", "Cluster2"))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(logcpm, design)

contrast.matrix <- makeContrasts(
  C2_vs_C1 = Cluster2 - Cluster1,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
# limma-trend: accounts for the mean-variance relationship in logCPM data.
# (Full limma-voom would require raw counts; trend is the correct choice for a
#  precomputed logCPM matrix.)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

de <- topTable(fit2, coef = "C2_vs_C1", number = Inf, sort.by = "P")
de$gene <- rownames(de)

write.csv(de, "04_results/de/DE_limma_cluster2_vs_cluster1_full.csv",
          row.names = FALSE)

# ========== 5. Select significant genes ==========
sig_up   <- subset(de, adj.P.Val < 0.05 & logFC >  1)
sig_down <- subset(de, adj.P.Val < 0.05 & logFC < -1)

cat("Significant UP genes:",   nrow(sig_up),   "\n")
cat("Significant DOWN genes:", nrow(sig_down), "\n")

write.csv(sig_up,   "04_results/de/DE_sig_up_cluster2.csv",   row.names = FALSE)
write.csv(sig_down, "04_results/de/DE_sig_down_cluster2.csv", row.names = FALSE)

# ---- helper: map ENSEMBL -> ENTREZID, de-duplicated ----
map_to_entrez <- function(ens_ids) {
  if (length(ens_ids) == 0) return(data.frame(ENSEMBL=character(0), ENTREZID=character(0), SYMBOL=character(0)))
  m <- AnnotationDbi::select(
    org.Hs.eg.db, keys = ens_ids, keytype = "ENSEMBL",
    columns = c("ENTREZID", "SYMBOL")
  )
  m <- m[!is.na(m$ENTREZID), ]
  m <- m[!duplicated(m$ENTREZID), ]
  m
}

# ---- helper: run + save GO BP and KEGG with empty-set guards ----
run_enrichment <- function(entrez, tag, title_suffix) {
  if (length(entrez) == 0) {
    cat("No genes to enrich for", tag, "- skipping.\n")
    return(invisible(NULL))
  }
  # GO BP
  ego <- enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                  ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  ego_df <- as.data.frame(ego)
  write.csv(ego_df, sprintf("04_results/go/GO_BP_enrichment_cluster2_%s.csv", tag),
            row.names = FALSE)
  if (!is.null(ego) && nrow(ego_df) > 0) {
    png(sprintf("04_results/go/GO_BP_dotplot_cluster2_%s.png", tag),
        width = 1200, height = 800)
    print(dotplot(ego, showCategory = 15) +
            ggtitle(paste("GO Biological Process Enrichment", title_suffix)))
    dev.off()
  } else cat("GO BP empty for", tag, "- no dotplot.\n")

  # KEGG
  ekegg <- enrichKEGG(gene = entrez, organism = "hsa", pvalueCutoff = 0.05)
  ekegg_df <- as.data.frame(ekegg)
  write.csv(ekegg_df, sprintf("04_results/go/KEGG_enrichment_cluster2_%s.csv", tag),
            row.names = FALSE)
  if (!is.null(ekegg) && nrow(ekegg_df) > 0) {
    png(sprintf("04_results/go/KEGG_dotplot_cluster2_%s.png", tag),
        width = 1200, height = 800)
    print(dotplot(ekegg, showCategory = 15) +
            ggtitle(paste("KEGG Pathway Enrichment", title_suffix)))
    dev.off()
  } else cat("KEGG empty for", tag, "- no dotplot.\n")
}

# ========== 6. Enrichment: UPregulated ==========
mapped_up <- map_to_entrez(sig_up$gene)
write.csv(mapped_up, "04_results/go/mapping_cluster2_up_ensembl_to_entrez.csv",
          row.names = FALSE)
run_enrichment(mapped_up$ENTREZID, "up", "(Cluster 2 - Upregulated)")

# ========== 7. Enrichment: DOWNregulated ==========
mapped_down <- map_to_entrez(sig_down$gene)
write.csv(mapped_down, "04_results/go/mapping_cluster2_down_ensembl_to_entrez.csv",
          row.names = FALSE)
run_enrichment(mapped_down$ENTREZID, "down", "(Cluster 2 - Downregulated)")

# ========== 8. Volcano Plot ==========
de$sig_label <- "Not Significant"
de$sig_label[de$adj.P.Val < 0.05 & de$logFC >  1] <- "Upregulated"
de$sig_label[de$adj.P.Val < 0.05 & de$logFC < -1] <- "Downregulated"
de$sig_label <- factor(de$sig_label,
                       levels = c("Upregulated", "Downregulated", "Not Significant"))

top_up    <- head(sig_up[order(-sig_up$logFC), ], 10)
top_down  <- head(sig_down[order(sig_down$logFC), ], 10)
labels_df <- rbind(top_up, top_down)

p_volcano <- ggplot(de, aes(x = logFC, y = -log10(adj.P.Val), color = sig_label)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c(
    "Upregulated"     = "#E74C3C",
    "Downregulated"   = "#2980B9",
    "Not Significant" = "#BDC3C7"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed", color = "grey40") +
  geom_text_repel(data = labels_df, aes(label = gene),
                  size = 2.5, max.overlaps = 15, color = "black") +
  labs(title = "Volcano Plot: Cluster 2 vs Cluster 1 (limma-trend)",
       x = expression(log[2]~"Fold Change"),
       y = expression(-log[10]~"(adj. p-value)"),
       color = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

ggsave("04_results/de/volcano_cluster2_vs_cluster1.png",
       plot = p_volcano, width = 9, height = 6, dpi = 150)

cat("\nSTEP 5 FINISHED.\n")