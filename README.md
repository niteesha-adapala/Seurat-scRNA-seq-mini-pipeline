# Single-Cell RNA-seq (Seurat) Mini-Pipeline

This project applies a **Seurat-based single-cell RNA-seq pipeline** to 10x Genomics data (mouse brain multiome, Alzheimer’s dataset). The workflow covers **QC, filtering, normalization, feature selection, dimensionality reduction, clustering, and visualization**.

## 📦 Repo contents
- `SingleCell_Mouse_AL_RNAseq.R` - end-to-end Seurat script
- `outputs/` - figures etc. (git-ignored)
- `20k_NSCLC_DTC_3p_nextgem_intron_donor_1_count_sample_feature_bc_matrix.h5/` - Input file

## 🔧 Requirements
- R ≥ 4.2
- Packages: `Seurat`, `tidyverse`, `ggplot2`  
---

## 🔬 Workflow & Results

### 1. Quality Control (QC)
We calculated three QC metrics per cell:
- **nFeature_RNA**: number of genes detected per cell  
- **nCount_RNA**: number of UMIs (sequencing depth)  
- **percent.mt**: percentage of mitochondrial genes  

**Results:**  
- The violin plots show wide distributions for nFeature_RNA and nCount_RNA.  
- `percent.mt` was effectively zero in this dataset (likely pre-filtered).  
- Scatter plots confirm a strong correlation between sequencing depth and gene counts (R ≈ 0.88).

**Interpretation:**  
Cells with very low genes or extreme counts are likely low-quality or doublets, and are removed in later filtering.

---

### 2. Filtering
We applied thresholds:
- Keep cells with **200–2500 genes**  
- Exclude cells with **>5% mitochondrial reads**  

This step ensures only high-quality, biologically meaningful cells remain for downstream analysis.

---

### 3. Normalization & Highly Variable Genes (HVGs)
We normalized data using **log normalization** and identified the **top 2000 variable genes**.

**Results:**  
- HVG plot highlights genes with highest variance (e.g., *Bnc2, Vip, Reln*).  
- These features drive clustering and dimensionality reduction.

**Interpretation:**  
Highly variable genes are most informative for distinguishing cell populations.

---

### 4. Scaling & PCA
All genes were scaled, and PCA was run on HVGs.

**Results:**  
- PC heatmaps show gene loadings across cells for the first components.  
- Elbow plot indicates ~10–15 PCs capture most variation.

**Interpretation:**  
We retained the first 15 PCs for clustering and visualization.

---

### 5. Clustering
Using PCA embeddings:
- Constructed a shared nearest-neighbor (SNN) graph.  
- Applied **Louvain clustering** at multiple resolutions.  

**Results:**  
- At resolution = 0.1, ~13 clusters were detected.  
- UMAP visualization shows distinct subpopulations of cells.

**Interpretation:**  
These clusters likely correspond to major brain cell types (neurons, astrocytes, oligodendrocytes, microglia, etc.), but precise annotation requires marker genes.

---

## 📈 Example Plots

- **QC Violin & Scatter:** Cell quality distribution and nFeature–nCount correlation  
- **HVGs:** Highly variable gene selection  
- **PC Heatmap & Elbow Plot:** Dimensionality reduction diagnostics  
- **PCA & UMAP Cluster Maps:** Separation of clusters at different resolutions  

*(See `outputs/` folder for generated figures.)*

---

## 🚧 Next Steps

This project currently stops at **clustering and UMAP visualization**. The next phases will include:

1. **Doublet Detection:**  
   - Use tools like **DoubletFinder** or **scrublet** to remove artificial cell multiplets.

2. **Cluster Annotation via Biomarkers:**  
   - Identify marker genes with `FindAllMarkers()`.  
   - Match cluster-specific markers to known canonical markers (e.g., *Gfap* → astrocytes, *Snap25* → neurons, *Mog* → oligodendrocytes, *C1qa* → microglia).  
   - Rename clusters with `RenameIdents()` for biological interpretability.

3. **Optional Automated Annotation:**  
   - Explore tools like **SingleR** or **Azimuth** for reference-based cell type labeling.

---

## 📚 References
- Hao et al., **Seurat v4**: Integrated analysis of multimodal single-cell data  
- 10x Genomics datasets and documentation  

---



