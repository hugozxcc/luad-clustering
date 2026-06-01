# LUAD Clustering Analysis

This repository contains scripts and notebooks for a lung adenocarcinoma (LUAD) clustering workflow, including preprocessing, feature selection, clustering, enrichment analysis, and clinical association analysis.

## Repository Structure

- `step1_data_collection.R` — data collection script
- `step2_preprocessing.R` — preprocessing pipeline
- `03_analysis/step3_feature_selection.ipynb` — feature selection notebook
- `03_analysis/step4_hierarchical_clustering.ipynb` — hierarchical clustering notebook
- `step5_go_enrichment.R` — GO and KEGG enrichment analysis
- `step6_clinical_assoc.R` — clinical association analysis
- `step6_heatmap.R` — heatmap generation
- `step6_heatmap_readable.R` — readable heatmap variant
- `step6_results_summary.R` — summary output generation

## Notes

- Raw data, processed data, generated results, manuscript files, and reference PDFs are excluded from version control via `.gitignore`.
- The repository currently tracks code and environment metadata only.

## Requirements

- R
- Python environment defined in `requirements.txt`
- Jupyter Notebook for `.ipynb` files

## Typical Workflow

1. Run `step1_data_collection.R`
2. Run `step2_preprocessing.R`
3. Execute notebooks in `03_analysis/`
4. Run `step5_go_enrichment.R`
5. Run `step6_*.R` scripts for downstream analysis and reporting
