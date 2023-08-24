##### python snRNAseq analysis 2 #####
##### Now we have our data and clusters, it's time to play with them

### (0) ENVIRONMENT
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

### Get biomarkers
si.filt.markers <- FindAllMarkers(si.filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- si.filt.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)

VlnPlot(si.filt, features = c("ITGA2B"))
VlnPlot(si.filt, features = c("BUD13", "COL3A1"))


FeaturePlot(si.filt, features = c("BEST4", "GUCY2C", "CFTR")) + theme_light()
grid.arrange(umap2, feats)


DoHeatmap(si.filt, features = top10$gene) + NoLegend()

umap <- UMAPPlot(si.filt,
                 group.by = "seurat_clusters",
                 cols = colorRampPalette(brewer.pal(8, "Dark2"))(21)) +
  theme_light() +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2))) +
  ggtitle("UMAP clusters") +
  theme(axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title.y = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", margin=margin(0,0,0,0)),
        legend.margin=margin(t = 0, unit='cm'),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.5, 'lines'))

featureTheme <- theme(axis.title.x = element_blank(),
                      axis.text.x = element_text(size = 6),
                      axis.text.y = element_text(size = 6),
                      axis.title.y = element_text(size = 8),
                      legend.margin=margin(t = 0, unit='cm'),
                      plot.title = element_text(size = 8, face = "bold", margin=margin(0,0,0,0)),
                      legend.position = "right",
                      legend.text = element_text(size = 6),
                      legend.key.size = unit(0.5, 'lines'))

# KIAA1211 - 3/8
kiaa1211 <- FeaturePlot(si.filt, features = c("KIAA1211"), order = TRUE) +
  theme_light() + ggtitle("KIAA1211") + labs(y = "") + featureTheme

# BEST4 -6
best4 <- FeaturePlot(si.filt, features = c("BEST4"), order = TRUE) +
  theme_light() + ggtitle("BEST4") + labs(y = "") + featureTheme

# OLFM4 - 14
olfm4 <- FeaturePlot(si.filt, features = c("OLFM4"), order = TRUE) +
  theme_light() + ggtitle("OLFM4") + labs(y = "") + featureTheme

# NRG1 - 13
nrg1 <- FeaturePlot(si.filt, features = c("NRG1"), order = TRUE) +
  theme_light() + ggtitle("NRG1") + labs(y = "") + featureTheme

grid.arrange(kiaa1211, best4, olfm4, nrg1, nrow = 2)
lay = rbind(c(1, 1),
            c(1, 1),
            c(2, 3),
            c(4, 5))
grid.arrange(umap, kiaa1211, best4, olfm4, nrg1, layout_matrix = lay)

## Cluster pathway analysis
BiocManager::install("ReactomeGSA")
library(ReactomeGSA)
gsva_results <- analyse_sc_clusters(si.filt, verbose = TRUE)
pathway_expression <- pathways(gsva_results)
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))

max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))
max_difference$diff <- max_difference$max - max_difference$min
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

head(max_difference)

plot_gsva_pathway(gsva_results, pathway_id = rownames(max_difference)[15])
plot_gsva_heatmap(gsva_results, max_pathways = 15, margins = c(6,20))
gsva_result <- gsva_results
plot_gsva_pca(gsva_result)

gsva_paths <- gsva_result@results$Seurat$pathways
gsva_pathExp <- gsva_paths[,3:23]
row.names(gsva_pathExp) <- gsva_paths$Name

library(pheatmap)
pheatmap(gsva_pathExp[which(row.names(gsva_pathExp) %in% max_difference$name[1:40]),])

FeaturePlot(si.filt, features = c("OTOP2"), order = TRUE)

### Balloon/bubble chart
si.all.markers <- FindAllMarkers(si.filt, only.pos = FALSE, min.pct = 0, logfc.threshold = 0)
saveRDS(si.all.markers, file = "data/reanalyzeRes/allCells_allMarkers.rds")
top10 <- si.filt.markers %>% dplyr::group_by(cluster) %>% dplyr::top_n(n = 10, wt = avg_log2FC)

