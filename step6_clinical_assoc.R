cluster <- read.csv("02_data_processed/TCGA_LUAD_cluster_assignment.csv")
clinical <- read.csv("01_data_raw/TCGA_LUAD_clinical.csv")

merged <- merge(cluster, clinical, by.x = "sample_id", by.y = "submitter_id")

# Stage vs cluster
tab_stage <- table(merged$cluster, merged$ajcc_pathologic_stage)
print(tab_stage)
write.csv(tab_stage, "04_results/clinical_cluster_vs_stage.csv")

# Chi-square test (remove empty columns if any)
tab_stage2 <- tab_stage[, colSums(tab_stage) > 0, drop = FALSE]

test_stage <- chisq.test(tab_stage2)
print(test_stage)

sink("04_results/clinical_stage_chisq.txt")
print(tab_stage2)
print(test_stage)
sink()

cat("Clinical association finished.\n")