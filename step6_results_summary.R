de_full <- read.csv("04_results/de/DE_limma_cluster2_vs_cluster1_full.csv")

# Top 20 genes UP (cluster 2)
top_up <- de_full[order(-de_full$logFC), ][1:20, c("gene","logFC","adj.P.Val")]
write.csv(top_up, "04_results/de/top20_up_cluster2.csv", row.names = FALSE)

# Top 20 genes DOWN (cluster 2)
top_down <- de_full[order(de_full$logFC), ][1:20, c("gene","logFC","adj.P.Val")]
write.csv(top_down, "04_results/de/top20_down_cluster2.csv", row.names = FALSE)

head(read.csv("04_results/go/KEGG_enrichment_cluster2_up.csv")[, c("Description","p.adjust","Count")], 10)
clinical <- read.csv("01_data_raw/TCGA_LUAD_clinical.csv")
names(clinical)[1:30]