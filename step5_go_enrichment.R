############################################
# STEP 5: Differential Expression + GO Enrichment
# Project: TCGA-LUAD clustering
# Author: Fadil Nugroho
############################################

# ========== 0. Setup ==========
rm(list = ls())

# Working directory (pastikan ini folder thesis)
setwd(".")

# ========== 1. Load libraries ==========
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "limma",
  "edgeR",
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot"
), ask = FALSE, update = FALSE)

library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)


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

# ========== 4. Differential Expression (limma) ==========
group <- factor(cl$cluster,
                levels = c(1, 2),
                labels = c("Cluster1_major", "Cluster2_minor"))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design


fit <- lmFit(logcpm, design)

contrast.matrix <- makeContrasts(
  Minor_vs_Major = Cluster2_minor - Cluster1_major,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

de <- topTable(
  fit2,
  coef = "Minor_vs_Major",
  number = Inf,
  sort.by = "P"
)

de$gene <- rownames(de)


write.csv(de,
          "04_results/de/DE_limma_cluster2_vs_cluster1_full.csv",
          row.names = FALSE)

# ========== 5. Select significant genes ==========
sig_up <- subset(de, adj.P.Val < 0.05 & logFC > 1)
sig_down <- subset(de, adj.P.Val < 0.05 & logFC < -1)

cat("Significant UP genes:", nrow(sig_up), "\n")
cat("Significant DOWN genes:", nrow(sig_down), "\n")

write.csv(sig_up, "04_results/de/DE_sig_up_cluster2.csv", row.names = FALSE)
write.csv(sig_down, "04_results/de/DE_sig_down_cluster2.csv", row.names = FALSE)

# ========== 6. Gene ID mapping ==========
ens_up <- sig_up$gene

mapped_up <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ens_up,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
)

mapped_up <- mapped_up[!is.na(mapped_up$ENTREZID), ]
mapped_up <- mapped_up[!duplicated(mapped_up$ENTREZID), ]

entrez_up <- mapped_up$ENTREZID

write.csv(mapped_up,
          "04_results/go/mapping_cluster2_up_ensembl_to_entrez.csv",
          row.names = FALSE)

# ========== 7. GO Enrichment (Biological Process) ==========
ego_bp <- enrichGO(
  gene          = entrez_up,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_bp_df <- as.data.frame(ego_bp)
write.csv(go_bp_df,
          "04_results/go/GO_BP_enrichment_cluster2_up.csv",
          row.names = FALSE)

png("04_results/go/GO_BP_dotplot_cluster2_up.png", width = 1200, height = 800)
dotplot(ego_bp, showCategory = 15) +
  ggtitle("GO Biological Process Enrichment (Cluster 2 - Upregulated)")
dev.off()

# ========== 8. KEGG Enrichment ==========
ekegg <- enrichKEGG(
  gene = entrez_up,
  organism = "hsa",
  pvalueCutoff = 0.05
)

kegg_df <- as.data.frame(ekegg)
write.csv(kegg_df,
          "04_results/go/KEGG_enrichment_cluster2_up.csv",
          row.names = FALSE)

png("04_results/go/KEGG_dotplot_cluster2_up.png", width = 1200, height = 800)
dotplot(ekegg, showCategory = 15) +
  ggtitle("KEGG Pathway Enrichment (Cluster 2 - Upregulated)")
dev.off()

head(read.csv("04_results/go/GO_BP_enrichment_cluster2_up.csv")[, c("Description","p.adjust","Count")])