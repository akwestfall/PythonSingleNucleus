##### python snRNAseq analysis 3 #####
##### Cell sub-clustering 101 -- specifically, enterocytes

### (0) ENVIRONMENT ######
library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(magrittr)
library(viridis)
library(grid)
library(gridExtra)
library(RColorBrewer)

setwd("~/Dropbox/CastoeLabFolder/projects/snakephysiolremod/_Manuscripts/_Python_earlySignaling_MS/snRNAseq/")

si.filt <- readRDS("data/reanalyzeRes/allCells_Filtered_Normalized.rds")

### (1) Separate cell groups
# ENTEROCYTES: 0, 1, 2, 3, 6, 8, 14 (stem), 20

enterocytes <- row.names(subset(x = si.filt,
                           subset = seurat_clusters %in% c(0, 1, 2, 3, 6, 8, 14, 20)))

si <- readRDS("raw_cells.RDS")
si[["percent.mt"]] <- PercentageFeatureSet(si, pattern = "^cds-YP|^rna-NC")

si.filt <- subset(x = si,
                  subset = percent.mt < 10)

si.ent <- si.filt[enterocytes,]
si.ent <- NormalizeData(si.ent, normalization.method = "LogNormalize", scale.factor = 10000)
si.ent <- FindVariableFeatures(si.ent, selection.method = "vst",
                                mean.function = ExpMean,
                                dispersion.function = LogVMR,
                                dispersion.cutoff = c(0.0125, 3))

### Scale data by expression and percentage mitochondrial expression
si.ent <- ScaleData(object = si.ent,
                     vars.to.regress = c("nCount_RNA", "percent.mt"))

### PCA
si.ent <- RunPCA(si.ent, features = VariableFeatures(si.ent))
PCAPlot(si.ent)

print(si.ent[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(si.ent, dims = 1:2, reduction = "pca")
DimPlot(si.ent, reduction = "pca")
DimHeatmap(si.ent, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(si.ent, dims = 1:15, cells = 500, balanced = TRUE)

### Determine dimensionality
si.ent <- JackStraw(si.ent, num.replicate = 100) # takes some time
si.ent <- ScoreJackStraw(si.ent, dims = 1:20)
JackStrawPlot(si.ent, dims = 1:20)
ElbowPlot(si.ent) ### elbow seems to occur around 11 for si.ent

### Cluster cells
si.ent <- FindNeighbors(si.ent, dims = 1:11)
si.ent <- FindClusters(si.ent, resolution = 0.5)
head(Idents(si.ent), 5)

### tSNE
si.ent <- RunTSNE(object = si.ent, dims.use = 1:11)
ident <- TSNEPlot(si.ent, group.by = "orig.ident", cols = plasma(4))
clust <- TSNEPlot(si.ent, group.by = "seurat_clusters")
grid.arrange(ident, clust, nrow = 1)

### UMAP
si.ent <- RunUMAP(si.ent, dims = 1:11)
umap1 <- UMAPPlot(si.ent, group.by = "orig.ident", cols = plasma(4))
umap2 <- UMAPPlot(si.ent, group.by = "seurat_clusters")
grid.arrange(umap1, umap2, nrow = 1)

FeaturePlot(object = si.ent, features = c("MUC1", "XBP1", "BEST4", "CHGA"))

saveRDS(si.ent, file = "data/reanalyzeRes/enterocytes_Filtered_Normalized.rds")
si.ent <- readRDS("data/reanalyzeRes/allCells_Filtered_Normalized.rds")

si.ent.markers <- FindAllMarkers(si.ent, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- si.ent.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)
cluster10.markers <- FindMarkers(si.ent, ident.1 = 10, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)