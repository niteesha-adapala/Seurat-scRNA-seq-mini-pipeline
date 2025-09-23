# Single-Cell RNA-seq (Seurat) Mini-Pipeline

A lightweight Seurat workflow for QC → normalization → HVGs → scaling → PCA → clustering → UMAP on 10x Genomics data.

## 📦 Repo contents
- `scripts/scrna_pipeline.R` — end-to-end Seurat script
- `outputs/` — figures etc. (git-ignored)

## 🔧 Requirements
- R ≥ 4.2
- Packages: `Seurat`, `tidyverse`, `ggplot2`  
  ```r
  install.packages(c("Seurat","tidyverse","ggplot2"))


