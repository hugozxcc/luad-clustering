############################################
# STEP 6.3: Clinical Association Analysis
# TCGA-LUAD: Cluster vs AJCC Pathologic Stage
# Note: clinical is patient-level; expression/cluster is sample-level.
############################################

rm(list = ls())
setwd(".")

# ---------- Load data ----------
cluster  <- read.csv("02_data_processed/TCGA_LUAD_cluster_assignment.csv", stringsAsFactors = FALSE)
clinical <- read.csv("01_data_raw/TCGA_LUAD_clinical.csv", stringsAsFactors = FALSE)

cat("cluster rows:",  nrow(cluster),  "\n")
cat("clinical rows:", nrow(clinical), "\n")

# ---------- Convert sample_id -> patient_id ----------
# TCGA patient barcode = first 12 chars (e.g., TCGA-XX-YYYY)
cluster$patient_id <- substr(cluster$sample_id, 1, 12)

# De-duplicate: some patients have >1 sample; keep one row per patient so the
# contingency table counts patients, not aliquots.
n_before <- nrow(cluster)
cluster  <- cluster[!duplicated(cluster$patient_id), ]
cat("Dropped", n_before - nrow(cluster), "duplicate-patient rows; kept",
    nrow(cluster), "\n")

overlap_patient <- sum(cluster$patient_id %in% clinical$submitter_id)
cat("overlap patient_id:", overlap_patient, "\n")
if (overlap_patient == 0) {
  stop("No overlap between cluster patient_id and clinical submitter_id. Check ID formats.")
}

# ---------- Merge ----------
merged <- merge(cluster, clinical, by.x = "patient_id", by.y = "submitter_id")
cat("merged rows:", nrow(merged), "\n")
if (nrow(merged) == 0) stop("Merged table is empty. Merge failed.")

# ---------- Stage vs Cluster (with NA) ----------
tab_stage_all <- table(merged$cluster, merged$ajcc_pathologic_stage, useNA = "ifany")
cat("\nStage table (including NA):\n"); print(tab_stage_all)

dir.create("04_results", showWarnings = FALSE)
write.csv(tab_stage_all, "04_results/clinical_cluster_vs_stage_withNA.csv")

# ---------- Clean: remove NA/blank stage ----------
merged_stage <- merged[!is.na(merged$ajcc_pathologic_stage) &
                         merged$ajcc_pathologic_stage != "", ]

# ---------- Collapse substages (IA/IB -> I, IIA/IIB -> II, ...) ----------
# A 2x4 main-stage table has far better expected cell counts than the
# fragmented ~8-column substage table, making the association test valid.
merged_stage$stage_main <- gsub("Stage ", "", merged_stage$ajcc_pathologic_stage)
merged_stage$stage_main <- gsub("[ABC]$", "", merged_stage$stage_main)  # strip A/B/C suffix

tab_stage <- table(merged_stage$cluster, merged_stage$stage_main)
tab_stage <- tab_stage[, colSums(tab_stage) > 0, drop = FALSE]
cat("\nStage table (collapsed main stages):\n"); print(tab_stage)
write.csv(tab_stage, "04_results/clinical_cluster_vs_stage.csv")

# ---------- Statistical test ----------
sink("04_results/clinical_stage_tests.txt")
cat("Clinical association: cluster vs AJCC pathologic stage (main stages)\n\n")
cat("Merged rows (patients):", nrow(merged), "\n")
cat("Rows after stage cleaning:", nrow(merged_stage), "\n\n")
cat("Contingency table (collapsed):\n"); print(tab_stage); cat("\n")

# Fisher's exact (Monte Carlo) - PRIMARY reported test.
# Robust to small expected counts in a 2x4 table.
cat("Fisher's Exact Test (Monte Carlo, B=10000) -- PRIMARY:\n")
set.seed(42)
fisher_res <- fisher.test(tab_stage, simulate.p.value = TRUE, B = 10000)
print(fisher_res)
cat("\n")

# Chi-square - reference only (expected counts may still be modest).
cat("Pearson Chi-squared test (reference only):\n")
chi_res <- suppressWarnings(chisq.test(tab_stage))
print(chi_res)
cat("\nExpected counts summary:\n")
print(summary(as.vector(chi_res$expected)))
sink()

cat("\nSaved:\n")
cat("- 04_results/clinical_cluster_vs_stage_withNA.csv\n")
cat("- 04_results/clinical_cluster_vs_stage.csv\n")
cat("- 04_results/clinical_stage_tests.txt\n")
cat("\nSTEP 6.3 FINISHED SUCCESSFULLY\n")