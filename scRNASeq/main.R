setwd("C:/Users/lalit/Desktop/scRNA-seq")
library(dplyr)
library(Seurat)
library(patchwork)

pbmc.data <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# Calculate mitochondrial QC metrics with the PercentageFeatureSet() function. Use the set of all genes starting with MT- as a set of mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter cells that have unique feature counts over 2,500 or less than 200, and filter cells that have >5% mitochondrial counts - INFORMED using violin counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Normalize
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features (feature selection), top 2000 most variable features for later PCA
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Perform linear dimensional reduction. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
DimPlot(pbmc, reduction = "pca")

#  DimHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. 
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Implement a resampling test inspired by the JackStraw procedure.
pbmc <- JackStraw(pbmc, num.replicate = 100) # This step may take 3 minutes.
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)

# Cluster the cells using first 10 PCs
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# Finding differentially expressed features (cluster biomarkers

# Find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%  group_by(cluster) %>%  slice_max(n = 2, order_by = avg_log2FC)

# Feature marker genes
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A")) # This is a single-line command.

# Use canonical markers to easily match the unbiased clustering to known cell types

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",                     "NK", "DC", "Platelet") 

names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#examine hetergeneous expression of the top 10 expressed genes in the dataset​

DotPlot(pbmc, features = top10, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

# Find differentially expressed features between CD14+ Monocytes and all other cells, only

# search for positive markers
monocyte.de.markers <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = NULL, only.pos = TRUE)

# view results
head(monocyte.de.markers)
View(monocyte.de.markers)
monocyte_genes <- rownames(monocyte.de.markers)

#Extract gene names of top 20 differentially expressed CD14+ monocytes
monocyte_top20 <- monocyte_genes[1:20]
DotPlot(pbmc, features = monocyte_top20, cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()
