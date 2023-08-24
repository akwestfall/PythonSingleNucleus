##### python snRNAseq analysis #####
##### I have no idea what I am doing but we are going to try some things

### (0) ENVIRONMENT
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(magrittr)
library(viridis)
library(grid)
library(gridExtra)

setwd("~/Dropbox/CastoeLabFolder/projects/snakephysiolremod/_Manuscripts/__Python_scRNA_Atlas_MS/snRNAseq/")

### (1) TUTORIAL/VIGNETTE #####
##### Data should be gzipped
ar14.data <- Read10X(data.dir = "data/originalRes/AR14/filtered_feature_bc_matrix/")
ar14 <- CreateSeuratObject(counts = ar14.data, project = "ar14", min.cells = 3, min.features = 200)

ar14[["percent.mt"]] <- PercentageFeatureSet(ar14, pattern = "^MT-")
VlnPlot(ar14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", ncol = 3))
FeatureScatter(ar14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

### Normalization
ar14 <- NormalizeData(ar14, normalization.method = "LogNormalize", scale.factor = 10000)

### Identify highly variable features
ar14 <- FindVariableFeatures(ar14, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(ar14), 10)

plot1 <- VariableFeaturePlot(ar14)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

### Scale data
all.genes <- rownames(ar14)
ar14 <- ScaleData(ar14, features = all.genes)

### Linear dimensional reduction
ar14 <- RunPCA(ar14, features = VariableFeatures(object = ar14))
print(ar14[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(ar14, dims = 1:2, reduction = "pca")
DimPlot(ar14, reduction = "pca")
DimHeatmap(ar14, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(ar14, dims = 1:15, cells = 500, balanced = TRUE)

### Determine dimensionality
ar14 <- JackStraw(ar14, num.replicate = 100)
ar14 <- ScoreJackStraw(ar14, dims = 1:20)
JackStrawPlot(ar14, dims = 1:20)
ElbowPlot(ar14) ### elbow seems to occur around 13-15 for AR14

### Cluster cells
ar14 <- FindNeighbors(ar14, dims = 1:13)
ar14 <- FindClusters(ar14, resolution = 0.5)
head(Idents(ar14), 5)

### Non-linear dimensional reduction
ar14 <- RunUMAP(ar14, dims = 1:13)
DimPlot(ar14, reduction = "umap")
saveRDS(ar14, file = "data/AR14/ar14_SeuratTutorial.rds")

### Identify Biomarkers
ar14.markers <- FindAllMarkers(ar14, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- ar14.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
cluster0.markers <- FindMarkers(ar14, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(ar14, features = c("ZEB2", "EMCN"))
VlnPlot(ar14, features = c("APOA4", "PTPRC"), slot = "counts", log = TRUE)
FeaturePlot(ar14, features = c("LGR5", "MSI1", "SOX9", "BMI1", "BMPR1A", "LOC103059990", "SMAD4", "MFSD2A", "APOA4"))
DoHeatmap(ar14, features = top10$gene) + NoLegend()

write.csv(top10, file = "data/AR14/top10_clusterMarkers.csv")

### (2) DATA ######
samples <- c("AQ16", "AR4", "AP22", "AR14")

for (sample in samples) {
  print(sample)
  assign(paste0(sample,".data"),
         Read10X(data.dir = paste0("data/reanalyzeRes/", sample, "/outs/filtered_feature_bc_matrix/")))
  assign(sample,
         CreateSeuratObject(counts = get(paste0(sample, ".data")), project = sample, min.cells = 3, min.features = 200))
}

VlnPlot(AQ16, features = c("nFeature_RNA", "nCount_RNA"))
FeatureScatter(AQ16, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(AR4, features = c("nFeature_RNA", "nCount_RNA"))
FeatureScatter(AR4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(AP22, features = c("nFeature_RNA", "nCount_RNA"))
FeatureScatter(AP22, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(AR14, features = c("nFeature_RNA", "nCount_RNA"))
FeatureScatter(AR14, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

### Merge
si <- merge(x = AQ16, y = list(AR4, AP22, AR14), add.cell.ids = samples, project = "si")
head(colnames(si))
table(si$orig.ident)

saveRDS(si, file = "raw_cells.RDS")

### Filter out MT genes and cells with >9% MT gene expression
si[["percent.mt"]] <- PercentageFeatureSet(si, pattern = "^cds-YP|^rna-NC")

VlnPlot(si, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
FeatureScatter(si, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggplot(si@meta.data, aes(x = `nCount_RNA`, y = `percent.mt`)) +
  geom_point(size = 0.5, alpha = 0.05) +
  geom_hline(yintercept = 9, linetype = "dashed", colour = "red") +
  theme_classic()

si.filt <- subset(x = si,
                  subset = percent.mt < 10)

### Normalize
si.filt <- NormalizeData(si.filt, normalization.method = "LogNormalize", scale.factor = 10000)

### Identify highly variable features
si.filt <- FindVariableFeatures(si.filt, selection.method = "vst",
                                mean.function = ExpMean,
                                dispersion.function = LogVMR,
                                dispersion.cutoff = c(0.0125, 3))

### Scale data by expression and percentage mitochondrial expression
si.filt <- ScaleData(object = si.filt,
                     vars.to.regress = c("nCount_RNA", "percent.mt"))

### PCA
si.filt <- RunPCA(si.filt, features = VariableFeatures(si.filt))
PCAPlot(si.filt)

print(si.filt[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(si.filt, dims = 1:2, reduction = "pca")
DimPlot(si.filt, reduction = "pca")
DimHeatmap(si.filt, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(si.filt, dims = 1:15, cells = 500, balanced = TRUE)

### Determine dimensionality
si.filt <- JackStraw(si.filt, num.replicate = 100, dims = 100) # takes some time
si.filt <- ScoreJackStraw(si.filt, dims = 1:50)
JackStrawPlot(si.filt, dims = 1:50)
ElbowPlot(si.filt) ### elbow seems to occur around 16 for si.filt

### Cluster cells
si.filt <- FindNeighbors(si.filt, dims = 1:50)
si.filt <- FindClusters(si.filt, resolution = 0.5)
head(Idents(si.filt), 5)

### tSNE
si.filt <- RunTSNE(object = si.filt, dims.use = 1:50)
ident <- TSNEPlot(si.filt, group.by = "orig.ident", cols = plasma(4))
clust <- TSNEPlot(si.filt, group.by = "seurat_clusters")
grid.arrange(ident, clust, nrow = 1)

### UMAP
si.filt <- RunUMAP(si.filt, dims = 1:50)
umap1 <- UMAPPlot(si.filt, group.by = "orig.ident", cols = plasma(4))
umap2 <- UMAPPlot(si.filt, group.by = "seurat_clusters")
grid.arrange(umap1, umap2, nrow = 1)

FeaturePlot(object = si.filt, features = c("MUC1", "XBP1", "VIL1", "CHGA"))

saveRDS(si.filt, file = "data/reanalyzeRes/allCells_Filtered_Normalized.rds")
si.filt <- readRDS("~/Dropbox/CastoeLabFolder/projects/snakephysiolremod/_Manuscripts/__Python_scRNA_Atlas_MS/snRNAseq/data/allCells_Filtered_Normalized.rds")

### Get biomarkers
si.filt.markers <- FindAllMarkers(si.filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- si.filt.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
cluster10.markers <- FindMarkers(si.filt, ident.1 = 10, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(si.filt, features = c("SLC13A1"))
VlnPlot(si.filt, features = c("APOA4", "PTPRC"))
FeaturePlot(si.filt, features = c("LGR5", "MSI1", "SOX9", "BMI1", "BMPR1A", "LOC103059990", "SMAD4", "MFSD2A", "APOA4"))
DoHeatmap(si.filt, features = top10$gene) + NoLegend()

write.csv(top10, file = "Seurat_results/allCells_top10_clusterMarkers.csv")
