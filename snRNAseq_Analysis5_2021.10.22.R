### Stacked Violin Plot for snRNAseq Analysis with Seurat
library(Seurat)
library(patchwork)
library(ggplot2)

### CODE BELOW FROM: https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
### and modified to suit my needs


## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + theme_light() +
    theme(legend.position = "none", 
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = 6, angle = 0), 
          axis.text.y = element_blank(), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

features <- c("ZPR1", "TREH", 
              "GUCY2C", "BEST4", "PTPRC", "MRC1", "TNXB",
              "MYH11", "OTOP2", "FLT3", "OLFM4", "GABBR2",
              "ANK1", "GCG", "ABCC9", "AMN", "CDH10", "PLEK")
features <- c("BMP1", "BMP2", "BMP3", "BMP4", "BMP5", "BMP6", "BMP7")
DoHeatmap(si.filt, features = features)
markers <- StackedVlnPlot(obj = si.filt, features = features, group.by = "orig.ident")
FeaturePlot(si.filt, features = features)

bmp1 <- FeaturePlot(si.filt, features = "BMP1") + theme_bw() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.25, 'cm'),
        legend.title=element_blank(),
        legend.margin=margin(c(0,1,0,0)),
        plot.margin=margin(c(0,0,0,1)),
        plot.title = element_text(face = "bold"))
bmp2 <- FeaturePlot(si.filt, features = "BMP2") + theme_bw() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.25, 'cm'),
        legend.title=element_blank(),
        legend.margin=margin(c(0,1,0,0)),
        plot.margin=margin(c(0,0,0,1)),
        plot.title = element_text(face = "bold"))
bmp3 <- FeaturePlot(si.filt, features = "BMP3") + theme_bw() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.25, 'cm'),
        legend.title=element_blank(),
        legend.margin=margin(c(0,1,0,0)),
        plot.margin=margin(c(0,0,0,1)),
        plot.title = element_text(face = "bold"))
bmp4 <- FeaturePlot(si.filt, features = "BMP4") + theme_bw() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.25, 'cm'),
        legend.title=element_blank(),
        legend.margin=margin(c(0,1,0,0)),
        plot.margin=margin(c(0,0,0,1)),
        plot.title = element_text(face = "bold"))
bmp5 <- FeaturePlot(si.filt, features = "BMP5") + theme_bw() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.25, 'cm'),
        legend.title=element_blank(),
        legend.margin=margin(c(0,1,0,0)),
        plot.margin=margin(c(0,0,0,1)),
        plot.title = element_text(face = "bold"))
bmp6 <- FeaturePlot(si.filt, features = "BMP6") + theme_bw() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.25, 'cm'),
        legend.title=element_blank(),
        legend.margin=margin(c(0,1,0,0)),
        plot.margin=margin(c(0,0,0,1)),
        plot.title = element_text(face = "bold"))
bmp7 <- FeaturePlot(si.filt, features = "BMP7") + theme_bw() +
  theme(text = element_text(size = 6),
        legend.key.width = unit(0.25, 'cm'),
        legend.title=element_blank(),
        legend.margin=margin(c(0,1,0,0)),
        plot.margin=margin(c(0,0,0,1)),
        plot.title = element_text(face = "bold"))

grid.arrange(bmp1, bmp2, bmp3,
             bmp4, bmp5, bmp6,
             bmp7, StackedVlnPlot(obj = si.filt, features = features),
             layout_matrix = rbind(c(1, 2, 3),
                                   c(4, 5, 6),
                                   c(7, 8, 8)))

si.filt@meta.data$cat <- NA
si.filt@meta.data$cat[which(si.filt@meta.data$seurat_clusters %in% c(0,1,2,3,4,6,12,14,15,17,19),)] <- "epithelial"
si.filt@meta.data$cat[which(si.filt@meta.data$seurat_clusters %in% c(5,7,13,16,21),)] <- "hematopoietic"
si.filt@meta.data$cat[which(si.filt@meta.data$seurat_clusters %in% c(8),)] <- "endothelial"
si.filt@meta.data$cat[which(si.filt@meta.data$seurat_clusters %in% c(9),)] <- "neural"
si.filt@meta.data$cat[which(si.filt@meta.data$seurat_clusters %in% c(10,11,20,18),)] <- "mesenchymal"

library(wesanderson)
umapCol <- UMAPPlot(si.filt, group.by = "cat", cols = wes_palette("Moonrise3", n = 5)) +
  theme_light() +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("cell type categories")

TSNEPlot(si.filt, group.by = "cat", cols = wes_palette("Moonrise3", n = 5)) +
  theme_light() +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle("cell type categories")

ggplot(si.filt@meta.data, aes(fill=seurat_clusters, x=orig.ident)) +
  geom_bar(position = "fill") +
  scale_x_discrete(breaks = c("AQ16", "AR4", "AP22", "AR14")) +
  theme_light() +
  ggtitle("proportion of cell types by sample") +
  xlab("sample")

epithelial <- subset(x = si.filt, subset = cat == "epithelial")
saveRDS(epithelial, file="results/seuratOutputs/epithelialCells.rds")
