library(data.table)

project_dir <- getwd()
raw_dir <- "01_data_raw"
proc_dir <- "02_data_processed"
dir.create(proc_dir, showWarnings = FALSE)

counts_path <- file.path(raw_dir, "TCGA_LUAD_counts.csv")
counts_df <- fread(counts_path)

gene_id <- counts_df[[1]]
counts_mat <- as.matrix(counts_df[, -1, with = FALSE])
rownames(counts_mat) <- gene_id

storage.mode(counts_mat) <- "numeric"

cat("Raw counts dim:", dim(counts_mat), "\n")
cat("Any NA:", any(is.na(counts_mat)), "\n")

rownames(counts_mat) <- sub("\\..*$", "", rownames(counts_mat))

dup <- duplicated(rownames(counts_mat))
cat("Duplicated genes after trimming:", sum(dup), "\n")

if (sum(dup) > 0) {
  # aggregate by gene id (sum)
  dt <- as.data.table(counts_mat, keep.rownames = "gene")
  dt_agg <- dt[, lapply(.SD, sum), by = gene]
  counts_mat <- as.matrix(dt_agg[, -1])
  rownames(counts_mat) <- dt_agg$gene
}
cat("After de-dup dim:", dim(counts_mat), "\n")

# 3) Library size
lib_size <- colSums(counts_mat)
cat("Library size summary:\n")
print(summary(lib_size))

# 4) CPM
cpm <- t(t(counts_mat) / lib_size * 1e6)

# 5) Filter rule: CPM >= 1 in at least 10% samples
min_prop <- 0.10
min_samples <- ceiling(ncol(counts_mat) * min_prop)

keep <- rowSums(cpm >= 1) >= min_samples
cat("Genes kept:", sum(keep), "out of", nrow(counts_mat), "\n")

counts_filt <- counts_mat[keep, ]

# Saving filtering info
filter_info <- data.frame(
  total_genes = nrow(counts_mat),
  kept_genes  = nrow(counts_filt),
  removed_genes = nrow(counts_mat) - nrow(counts_filt),
  min_samples = min_samples,
  rule = paste0("CPM>=1 in >= ", min_samples, " samples (", min_prop*100, "%)")
)

write.csv(filter_info, file.path(proc_dir, "TCGA_LUAD_gene_filter_info.csv"), row.names = FALSE)

# 6) logCPM
cpm_filt <- t(t(counts_filt) / colSums(counts_filt) * 1e6)
logcpm <- log2(cpm_filt + 1)

cat("logCPM dim:", dim(logcpm), "\n")
cat("logCPM range:", range(logcpm), "\n")

# Save normalized matrix
write.csv(logcpm, file.path(proc_dir, "TCGA_LUAD_norm_logCPM.csv"))
saveRDS(logcpm, file.path(proc_dir, "TCGA_LUAD_norm_logCPM.rds"))

# 7) Sanity checks
cat("Any NA in logCPM:", any(is.na(logcpm)), "\n")
cat("Any Inf in logCPM:", any(is.infinite(logcpm)), "\n")

# cek beberapa statistik
print(summary(as.vector(logcpm)))