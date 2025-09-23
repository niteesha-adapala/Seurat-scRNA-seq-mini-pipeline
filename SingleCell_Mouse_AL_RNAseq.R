#install.packages("ggplot2")                   
#install.packages(c("hdf5r","Matrix"))
#install.packages('SeuratDisk')
library(Seurat)
library(tidyverse)
library(ggplot2)

nsclc.sparse.m <- Read10X_h5(filename = "C:\Users/Niteesha/Desktop/mini-pro/20k_NSCLC_DTC_3p_nextgem_intron_donor_1_count_sample_feature_bc_matrix.h5")

## Load the dataset

mats <- Seurat::Read10X_h5(
  "C:/Users/Niteesha/Downloads/Multiome_RNA_ATAC_Mouse_Brain_Alzheimers_AppNote_filtered_feature_bc_matrix.h5"
)
nsclc.sparse.m -> mats
obj <- CreateSeuratObject(counts = if (is.list(mats)) mats[[1]] else mats)
mats -> mouseAL

str(mouseAL)
cts <-  mouseAL$`Gene Expression`
head (cts, n=10)

## Initialize the Seurat object with the raw (non-normalized data).
mouseAL.obj <- CreateSeuratObject(counts = cts, project = "SingleCell", min.cells = 3, min.features = 200)
str(mouseAL.obj)
mouseAL.obj


# 1. QC
View(mouseAL.obj@meta.data)

# % MT reads
mouseAL.obj[["percent.mt"]] <- PercentageFeatureSet(mouseAL.obj, pattern = "^MT-")
View(mouseAL.obj@meta.data)

VlnPlot(mouseAL.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(mouseAL.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')


# 2. Filtering 
mouseAL.obj <- subset(mouseAL.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)
mouseAL.obj

# 3. Normalize data 

mouseAL.obj <- NormalizeData(mouseAL.obj)
str(mouseAL.obj)


# 4. Identify highly variable features 
mouseAL.obj <- FindVariableFeatures(mouseAL.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(mouseAL.obj), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(mouseAL.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


# 5. Scaling 
all.genes <- rownames(mouseAL.obj)
mouseAL.obj <- ScaleData(mouseAL.obj, features = all.genes)

str(mouseAL.obj)

# 6. Perform Linear dimensionality reduction 
mouseAL.obj <- RunPCA(mouseAL.obj, features = VariableFeatures(object = mouseAL.obj))

  # visualize PCA results
print(mouseAL.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(mouseAL.obj, dims = 1, cells = 500, balanced = TRUE)


# determine dimensionality of the data
ElbowPlot(mouseAL.obj)


# 7. Clustering 
mouseAL.obj <- FindNeighbors(mouseAL.obj, dims = 1:15)

# understanding resolution
mouseAL.obj <- FindClusters(mouseAL.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(mouseAL.obj@meta.data)

DimPlot(mouseAL.obj, group.by = "RNA_snn_res.0.1", label = TRUE)

# setting identity of clusters
Idents(mouseAL.obj)
Idents(mouseAL.obj) <- "RNA_snn_res.0.1"
Idents(mouseAL.obj)

# non-linear dimensionality reduction
mouseAL.obj <- RunUMAP(mouseAL.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label

# individual clusters
DimPlot(mouseAL.obj, reduction = "umap")

# Biomarkers:
#cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)


# find markers for all clusters
all.markers <- FindAllMarkers(
  mouseAL.obj, 
  only.pos = TRUE,   # only positive markers (upregulated in cluster)
  min.pct = 0.25,    # expressed in at least 25% of cells in cluster
  logfc.threshold = 0.25
)

# view top markers per cluster
head(all.markers)


top10 <- all.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 10, wt = avg_log2FC)

# heatmap of top 10 markers per cluster
DoHeatmap(mouseAL.obj, features = top10$gene) + NoLegend()

#Step4
new.cluster.ids <- c("Astrocytes", "Neurons", "Microglia", "Oligodendrocytes")

names(new.cluster.ids) <- levels(mouseAL.obj)
mouseAL.obj <- RenameIdents(mouseAL.obj, new.cluster.ids)

# plot again with labels
DimPlot(mouseAL.obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#Step 5
install.packages('SingleR')
library(SingleR)
install.packages('celldex')
library(celldex)

# use mouse RNA reference
ref <- MouseRNAseqData()
pred <- SingleR(test = GetAssayData(mouseAL.obj, slot = "data"), ref = ref, labels = ref$label.main)

# add predictions to metadata
mouseAL.obj$SingleR <- pred$labels

DimPlot(mouseAL.obj, group.by = "SingleR", label = TRUE)

