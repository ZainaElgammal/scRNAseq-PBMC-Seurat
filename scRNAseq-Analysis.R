# ==============================================================================
# Single-Cell RNA-seq Analysis of PBMCs using Seurat v5
# Script Name: scRNAseq-Analysis.R
# Author: Zaina Elgammal
# Description: Standard workflow for QC, normalization, clustering, and cell annotation
# ==============================================================================

# Load required packages
library(Seurat)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Load Dataset & Initialization
# ------------------------------------------------------------------------------
# Load PBMC 10x Genomics dataset (Downloaded from 10x Genomics repository)
# Note: Ensure dataset directory exists in local environment
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# 2. Quality Control & Filtering
# ------------------------------------------------------------------------------
# Calculate mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC Visualizations
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter low-quality cells and potential doublets
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# ------------------------------------------------------------------------------
# 3. Data Normalization & Feature Selection
# ------------------------------------------------------------------------------
# Global-scaling normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# ------------------------------------------------------------------------------
# 4. Dimensionality Reduction & Clustering
# ------------------------------------------------------------------------------
# Scale data and run PCA
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Graph-based clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# UMAP Visualization
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", label = TRUE)

# ------------------------------------------------------------------------------
# 5. Cell Type Annotation
# ------------------------------------------------------------------------------
# Assign canonical cell type markers
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", 
                     "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

# Save annotated UMAP plot
pdf("scRNAseq_PBMC_UMAP.pdf", width = 8, height = 6)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()