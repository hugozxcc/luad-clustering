# ===== 0. Setup packages =====
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor packages
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment"), ask = FALSE, update = FALSE)

# CRAN packages
install.packages(c("tidyverse", "data.table"), dependencies = TRUE)

# Load
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(data.table)

project_dir <- getwd()

raw_dir <- "01_data_raw"
dir.create(raw_dir, showWarnings = FALSE)

# 2) Query expression data (raw counts)
query_exp <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

# cek jumlah file
query_exp

# 3) Download
GDCdownload(
  query = query_exp,
  method = "api",
  files.per.chunk = 5
)

# 4) Prepare
data_se <- GDCprepare(query_exp)
data_se

saveRDS(data_se, file = file.path(raw_dir, "TCGA_LUAD_SE.rds"))

count_matrix <- assay(data_se)  # genes x samples
dim(count_matrix)

# optional tapi good
write.csv(count_matrix, file = file.path(raw_dir, "TCGA_LUAD_counts.csv"))

# 5) Clinical data
clinical <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
dim(clinical)
head(clinical, 3)

write.csv(clinical, file = file.path(raw_dir, "TCGA_LUAD_clinical.csv"), row.names = FALSE)

# 6) Sanity check
# A) cek dimensi
cat("Counts dimension:", dim(count_matrix), "\n")

# B) cek ada NA?
cat("NA in counts:", any(is.na(count_matrix)), "\n")

# C) cek ringkas distribusi counts
summary(as.vector(count_matrix))[1:6]

# D) cek contoh barcode sampel
head(colnames(count_matrix), 3)

# E) cek metadata sampel ada
head(colData(data_se)[, 1:5])