### CellChat by timepoint

# (0) Environment
library(CellChat)
library(patchwork)
library(Seurat)
library(cowplot)
options(stringsAsFactors = FALSE)

# (1) Data subsets
Idents(si.filt) <- si.filt$cell_type2

fasted <- subset(si.filt, subset = orig.ident == "AQ16")
x6hrpf <- subset(si.filt, orig.ident == "AR4")
x12hrpf <- subset(si.filt, orig.ident == "AP22")
x1dpf <- subset(si.filt, orig.ident == "AR14")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# (2) Timepoint analysis
### (a) fasted ######
data.input <- fasted@assays$RNA@data
meta <- fasted@meta.data
meta$labels <- meta$cell_type2 # cellchat doesn't like numeric labels
meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
cell.use <- rownames(meta)
data.input <- data.input[, cell.use]
meta <- meta[cell.use, ]
unique(meta$labels) # check the cell labels
cc.fasted <- createCellChat(object = data.input, meta = meta, group.by = "labels")

cc.fasted@DB <- CellChatDB.use

cc.fasted <- subsetData(cc.fasted) # This step is necessary even if using the whole database

cc.fasted <- identifyOverExpressedGenes(cc.fasted, group.by = "labels")
cc.fasted <- identifyOverExpressedInteractions(cc.fasted)
cc.fasted <- projectData(cc.fasted, PPI.human)
cc.fasted <- computeCommunProb(cc.fasted, type = "truncatedMean", trim = 0.1, raw.use = FALSE)
cc.fasted <- filterCommunication(cc.fasted, min.cells = 10)
cc.fasted <- computeCommunProbPathway(cc.fasted)
cc.fasted <- aggregateNet(cc.fasted)

groupSize <- as.numeric(table(cc.fasted@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc.fasted@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.fasted@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

saveRDS(cc.fasted, "~/Desktop/cc.fasted.rds")

### (b) 6hrpf ######
data.input <- x6hrpf@assays$RNA@data
meta <- x6hrpf@meta.data
meta$labels <- meta$cell_type2 # cellchat doesn't like numeric labels
meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
cell.use <- rownames(meta)
data.input <- data.input[, cell.use]
meta <- meta[cell.use, ]
unique(meta$labels) # check the cell labels
cc.6hrpf <- createCellChat(object = data.input, meta = meta, group.by = "labels")

cc.6hrpf@DB <- CellChatDB.use

cc.6hrpf <- subsetData(cc.6hrpf) # This step is necessary even if using the whole database

cc.6hrpf <- identifyOverExpressedGenes(cc.6hrpf, group.by = "labels")
cc.6hrpf <- identifyOverExpressedInteractions(cc.6hrpf)
cc.6hrpf <- projectData(cc.6hrpf, PPI.human)
cc.6hrpf <- computeCommunProb(cc.6hrpf, type = "truncatedMean", trim = 0.1, raw.use = FALSE)
cc.6hrpf <- filterCommunication(cc.6hrpf, min.cells = 10)
cc.6hrpf <- computeCommunProbPathway(cc.6hrpf)
cc.6hrpf <- aggregateNet(cc.6hrpf)

groupSize <- as.numeric(table(cc.6hrpf@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc.6hrpf@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.6hrpf@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

saveRDS(cc.6hrpf, "~/Desktop/cc.6hrpf.rds")

### (c) 12hrpf ######
data.input <- x12hrpf@assays$RNA@data
meta <- x12hrpf@meta.data
meta$labels <- meta$cell_type2 # cellchat doesn't like numeric labels
meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
cell.use <- rownames(meta)
data.input <- data.input[, cell.use]
meta <- meta[cell.use, ]
unique(meta$labels) # check the cell labels
cc.12hrpf <- createCellChat(object = data.input, meta = meta, group.by = "labels")

cc.12hrpf@DB <- CellChatDB.use

cc.12hrpf <- subsetData(cc.12hrpf) # This step is necessary even if using the whole database

cc.12hrpf <- identifyOverExpressedGenes(cc.12hrpf, group.by = "labels")
cc.12hrpf <- identifyOverExpressedInteractions(cc.12hrpf)
cc.12hrpf <- projectData(cc.12hrpf, PPI.human)
cc.12hrpf <- computeCommunProb(cc.12hrpf, type = "truncatedMean", trim = 0.1, raw.use = FALSE)
cc.12hrpf <- filterCommunication(cc.12hrpf, min.cells = 10)
cc.12hrpf <- computeCommunProbPathway(cc.12hrpf)
cc.12hrpf <- aggregateNet(cc.12hrpf)

groupSize <- as.numeric(table(cc.12hrpf@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc.12hrpf@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.12hrpf@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

saveRDS(cc.12hrpf, "~/Desktop/cc.12hrpf.rds")

### (d) 1dpf ######
data.input <- x1dpf@assays$RNA@data
meta <- x1dpf@meta.data
meta$labels <- meta$cell_type2 # cellchat doesn't like numeric labels
meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))
cell.use <- rownames(meta)
data.input <- data.input[, cell.use]
meta <- meta[cell.use, ]
unique(meta$labels) # check the cell labels
cc.1dpf <- createCellChat(object = data.input, meta = meta, group.by = "labels")

cc.1dpf@DB <- CellChatDB.use

cc.1dpf <- subsetData(cc.1dpf) # This step is necessary even if using the whole database

cc.1dpf <- identifyOverExpressedGenes(cc.1dpf, group.by = "labels")
cc.1dpf <- identifyOverExpressedInteractions(cc.1dpf)
cc.1dpf <- projectData(cc.1dpf, PPI.human)
cc.1dpf <- computeCommunProb(cc.1dpf, type = "truncatedMean", trim = 0.1, raw.use = FALSE)
cc.1dpf <- filterCommunication(cc.1dpf, min.cells = 10)
cc.1dpf <- computeCommunProbPathway(cc.1dpf)
cc.1dpf <- aggregateNet(cc.1dpf)

groupSize <- as.numeric(table(cc.1dpf@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc.1dpf@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.1dpf@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

saveRDS(cc.1dpf, "~/Desktop/cc.1dpf.rds")

cc.fasted <- readRDS("/Users/aundreawestfall/Desktop/Regeneration_Python/snRNAseq_python/python_cellCommunication/cc.fasted.rds")
cc.6hrpf <- readRDS("/Users/aundreawestfall/Desktop/Regeneration_Python/snRNAseq_python/python_cellCommunication/cc.6hrpf.rds")
cc.12hrpf <- readRDS("/Users/aundreawestfall/Desktop/Regeneration_Python/snRNAseq_python/python_cellCommunication/cc.12hrpf.rds")
cc.1dpf <- readRDS("/Users/aundreawestfall/Desktop/Regeneration_Python/snRNAseq_python/python_cellCommunication/cc.1dpf.rds")

### Centrality stuff
cc.fasted <- netAnalysis_computeCentrality(cc.fasted, slot.name = "netP")
cc.6hrpf <- netAnalysis_computeCentrality(cc.6hrpf, slot.name = "netP")
cc.12hrpf <- netAnalysis_computeCentrality(cc.12hrpf, slot.name = "netP")
cc.1dpf <- netAnalysis_computeCentrality(cc.1dpf, slot.name = "netP")



# (2) Visualization
groupSize00 <- as.numeric(table(cc.fasted@idents))
groupSize06 <- as.numeric(table(cc.6hrpf@idents))
groupSize12 <- as.numeric(table(cc.12hrpf@idents))
groupSize24 <- as.numeric(table(cc.1dpf@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cc.1dpf@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.1dpf@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

lacteal0 <- netVisual_bubble(cc.fasted, sources.use = c("BEST4+/SPIB+"), targets.use = c("lymphatic endothelial", "vascular endothelial", "myeloid", "lymphoid"), signaling = "NOTCH",
                 remove.isolate = TRUE) + theme(axis.text = element_text(size = 5))
lacteal6 <- netVisual_bubble(cc.6hrpf, sources.use = c("BEST4+/SPIB+"), targets.use = c("lymphatic endothelial", "vascular endothelial", "myeloid", "lymphoid"), signaling = "NOTCH",
                 remove.isolate = TRUE) + theme(axis.text = element_text(size = 5))
lacteal12 <- netVisual_bubble(cc.12hrpf, sources.use = c("BEST4+/SPIB+"), targets.use = c("lymphatic endothelial", "vascular endothelial", "myeloid", "lymphoid"), signaling = "NOTCH",
                 remove.isolate = TRUE) + theme(axis.text = element_text(size = 5))
lacteal24 <- netVisual_bubble(cc.1dpf, sources.use = c("BEST4+/SPIB+"), targets.use = c("lymphatic endothelial", "vascular endothelial", "myeloid", "lymphoid"), signaling = "NOTCH",
                 remove.isolate = TRUE) + theme(axis.text = element_text(size = 5))

lacteal0 <- lacteal0$data %>% mutate(sample = "fasted")
lacteal6 <- lacteal6$data %>% mutate(sample = "6hrpf")
lacteal12 <- lacteal12$data %>% mutate(sample = "12hrpf")
lacteal24 <- lacteal24$data %>% mutate(sample = "1dpf")

lacteal <- rbind(lacteal0, lacteal6) %>% rbind(lacteal12) %>% rbind(lacteal24) %>%
  mutate(sample = factor(sample, levels = c("fasted", "6hrpf", "12hrpf", "1dpf")))

lacteal %>% ggplot(aes(x = target, y = interaction_name_2, color = prob, group = sample)) + geom_point(size = 3) + facet_grid(cols = vars(sample)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.y = element_blank()) +
  scale_color_continuous(low = "gray85", high = "green3") +
  xlab("target cell") +
  ylab("ligand-receptor pair")
  



netVisual_bubble(cc.fasted, sources.use = c("BEST4+/SPIB+"), targets.use = c("stromal", "fibroblast", "smooth muscle", "mesothelial", "enterocyte"),
                 remove.isolate = TRUE) + theme(axis.text = element_text(size = 5)) + theme(legend.position = "none") +
  netVisual_bubble(cc.6hrpf, sources.use = c("BEST4+/SPIB+"), targets.use = c("stromal", "fibroblast", "smooth muscle", "mesothelial", "enterocyte"),
                   remove.isolate = TRUE) + theme(axis.text = element_text(size = 5)) + theme(legend.position = "none") +
  netVisual_bubble(cc.12hrpf, sources.use = c("BEST4+/SPIB+"), targets.use = c("stromal", "fibroblast", "smooth muscle", "mesothelial", "enterocyte"), 
                   remove.isolate = TRUE) + theme(axis.text = element_text(size = 5)) + theme(legend.position = "none") +
  netVisual_bubble(cc.1dpf, sources.use = c("BEST4+/SPIB+"), targets.use = c("stromal", "fibroblast", "smooth muscle", "mesothelial", "enterocyte"),
                   remove.isolate = TRUE) + theme(axis.text = element_text(size = 5)) + plot_layout(nrow = 1)
netVisual_chord_gene(cc.fasted, sources.use = c("stromal"), targets.use = c("enterocyte"),
                     slot.name = "netP", legend.pos.x = 10)

FeaturePlot(si.filt, features = c("LAMA1", "LAMA2", "LAMA3", "LAMA4", "LOC103055002"), order = TRUE)

par(mfrow = c(2,4), xpd=TRUE)
netVisual_circle(cc.fasted@net$count, vertex.weight = groupSize00, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.fasted@net$weight, vertex.weight = groupSize00, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cc.6hrpf@net$count, vertex.weight = groupSize06, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.6hrpf@net$weight, vertex.weight = groupSize06, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cc.12hrpf@net$count, vertex.weight = groupSize12, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.12hrpf@net$weight, vertex.weight = groupSize12, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_circle(cc.1dpf@net$count, vertex.weight = groupSize24, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cc.1dpf@net$weight, vertex.weight = groupSize24, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_chord_gene(cc.fasted, sources.use = c("BEST4+/SPIB+", "enterocyte"), targets.use = c("enterocyte", "stromal", "smooth muscle", "mesothelial", "fibroblast", "BEST4+/SPIB+"),
                     slot.name = "netP", legend.pos.x = 10, lab.cex = 0.5)
netVisual_chord_gene(cc.fasted, sources.use = c("BEST4+/SPIB+"), targets.use = c("enterocyte", "stromal", "smooth muscle", "mesothelial", "fibroblast", "BEST4+/SPIB+"),
                     slot.name = "netP", legend.pos.x = 10, lab.cex = 0.5)

notch1 <- netAnalysis_signalingRole_scatter(cc.12hrpf, signaling = c("NOTCH"))
gg1 <- netAnalysis_signalingRole_scatter(cc.12hrpf)
netAnalysis_signalingRole_network(cc.1dpf, signaling = "NOTCH", font.size = 10)
plot_grid(notch1, notch2)

### Centrality stuff
cc.fasted <- netAnalysis_computeCentrality(cc.fasted, slot.name = "netP")
cc.6hrpf <- netAnalysis_computeCentrality(cc.6hrpf, slot.name = "netP")
cc.12hrpf <- netAnalysis_computeCentrality(cc.12hrpf, slot.name = "netP")
cc.1dpf <- netAnalysis_computeCentrality(cc.1dpf, slot.name = "netP")


library(NMF)
library(ggalluvial)
selectK(cc.fasted, pattern = "incoming")

### outgoing n =
## fasted = 3
##  6hrpf = 3
## 12hrpf = 3
##   1dpf = 2

### incoming n =
## fasted = 3
##  6hrpf = 3
## 12hrpf = 1
##   1dpf = 3

cc.fasted <- identifyCommunicationPatterns(cc.fasted, pattern = "outgoing", k = 3, height = 12)
cc.fasted <- identifyCommunicationPatterns(cc.fasted, pattern = "incoming", k = 3, height = 12)

cc.6hrpf <- identifyCommunicationPatterns(cc.6hrpf, pattern = "outgoing", k = 3, height = 12)
cc.6hrpf <- identifyCommunicationPatterns(cc.6hrpf, pattern = "incoming", k = 4, height = 12)

cc.12hrpf <- identifyCommunicationPatterns(cc.12hrpf, pattern = "outgoing", k = 4, height = 12)
cc.12hrpf <- identifyCommunicationPatterns(cc.12hrpf, pattern = "incoming", k = 3, height = 12)

cc.1dpf <- identifyCommunicationPatterns(cc.1dpf, pattern = "outgoing", k = 5, height = 12)
cc.1dpf <- identifyCommunicationPatterns(cc.1dpf, pattern = "incoming", k = 4, height = 12)


netAnalysis_river(cc.6hrpf, pattern = "outgoing")
netAnalysis_river(cc.6hrpf, pattern = "incoming")



ht1 <- netAnalysis_signalingRole_heatmap(cc.fasted, pattern = "outgoing", height = 12)
ht2 <- netAnalysis_signalingRole_heatmap(cc.fasted, pattern = "incoming", height = 12)
ht1 + ht2

netAnalysis_river(cc.fasted, slot.name = "netP", pattern = "outgoing")
netAnalysis_dot(cc.fasted, pattern = "outgoing")


cellchat <- identifyCommunicationPatterns(cc.fasted, pattern = "incoming", k = 3)
identifyCommunicationPatterns(cc.fasted, pattern = "incoming", k = 3, height = 12)
netAnalysis_river(cellchat, pattern = "incoming")

gg00 <- netAnalysis_signalingRole_scatter(cc.fasted, signaling = c("BMP")) + ggtitle("fasted") + theme_bw()
gg06 <- netAnalysis_signalingRole_scatter(cc.6hrpf, signaling = c("BMP")) + ggtitle("6hrpf") + theme_bw()
gg12 <- netAnalysis_signalingRole_scatter(cc.12hrpf, signaling = c("BMP")) + ggtitle("12hrpf") + theme_bw()
gg24 <- netAnalysis_signalingRole_scatter(cc.1dpf, signaling = c("BMP")) + ggtitle("1dpf") + theme_bw()
gg00 + gg06 + gg12 + gg24



##### COMPARISON ANALYSIS VERSION 1.6.1
cc.fasted <- updateCellChat(cc.fasted)
cc.6hrpf <- updateCellChat(cc.6hrpf)
cc.12hrpf <- updateCellChat(cc.12hrpf)
cc.1dpf <- updateCellChat(cc.1dpf)
object.list <- list(fasted = cc.fasted, x6hrpf = cc.6hrpf)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
interactionWeight <- netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", vertex.label.cex = 0.5, edge.width.max = 3, title.name = "")

gg1 <- netVisual_heatmap(cellchat, font.size = 6, font.size.title = 7)
gg2 <- netVisual_heatmap(cellchat, measure = "weight", font.size = 6, font.size.title = 7)
gg1 + gg2 #### use this one!!!
row1 <- plot_grid(NULL, grid.grabExpr(draw(gg1, heatmap_legend_side = "bottom")), grid.grabExpr(draw(gg2, heatmap_legend_side = "bottom")), nrow = 1, rel_widths = c(1, 1.5, 1.5))

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "BEST4+/SPIB+")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "enterocyte")
patchwork::wrap_plots(plots = list(gg1,gg2))

cellchat <- computeNetSimilarityPairwise(cellchat, slot.name = "netP", type = c("functional"))
cellchat <- netEmbedding(cellchat, slot.name = "netP", type = "functional", umap.method = "uwot", pathway.remove = NULL)
cellchat <- netClustering(cellchat, slot.name = "netP", type = "functional", do.parallel = FALSE, k = 6)

functional <- cellchat@netP$similarity$functional$dr$`1-2` %>% as.data.frame
functional$group <- cellchat@netP$similarity$functional$group$`1-2`

functional <- functional %>% rownames_to_column() %>% as_tibble() %>%
  separate(rowname, c("pathway","sample"), sep = "--")

funcAll <- ggplot(functional, aes(x = UMAP1, y = UMAP2, color = factor(group), shape = sample, label = pathway)) + geom_point() +
  theme_bw() +
  ggtitle("Functional similarity of communication networks") +
  guides(shape = guide_legend("time point", keyheight = 0.4, title.position = "top", direction = "horizontal"),
         color = guide_legend("cluster", keyheight = 0.4, title.position = "top", direction = "horizontal")) +
  scale_shape_discrete(labels = c("fasted", "6hrpf")) +
  theme(legend.position = "top",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 6),
        legend.title = element_text(face = "bold")) 
  
pathways.show <- c("COLLAGEN") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}

#> guides(size = guide_legend("proportion cells", keyheight = 0.4, title.position = "top", direction = "horizontal"),
#> color = guide_colorbar("avg. expression", barwidth = 3, barheight = 0.5, title.position = "top", direction = "horizontal")) +
#> scale_x_discrete(labels = marker.names) +
#> scale_y_discrete(position = "right") + 
#>  theme(text = element_text(size = 6),
#>        plot.title = element_blank(),
#>        axis.title.y = element_blank(),
#>        axis.title.x = element_text(size = 7, face = "bold", hjust = 0),
#>        axis.text.x = element_text(size = 6),
#>        strip.text = element_text(size = 6, face = "bold", color = "white"),
#>        legend.position = "top", legend.justification = "right", legend.margin = margin(5, 5, -15, 5), legend.title = element_text(face = "bold"),
#>        plot.margin = margin(5, 5, 5, 5))

funcFac <- ggplot(functional, aes(x = UMAP1, y = UMAP2, color = factor(group), shape = sample, label = pathway)) +
  geom_point() +
  geom_label_repel(aes(fill = factor(group), segment.color = factor(group)), color = "white",
                   size = 2.14, box.padding = 0.1, label.padding = 0.1, 
                   max.overlaps = 4, force = 2, force_pull = 0.75, min.segment.length = 0.2) +
  theme_bw() +
  guides(shape = guide_legend("time point", keyheight = 0.4, title.position = "top", direction = "horizontal"),
         color = guide_legend("cluster", keyheight = 0.4, title.position = "top", direction = "horizontal")) +
  scale_shape_discrete(labels = c("fasted", "6hrpf")) +
  theme(legend.position = "none",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(size = 6),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  facet_wrap(vars(factor(group)), scales = "free", ncol = 3)

row2 <- plot_grid(funcAll, funcFac, ncol = 2, rel_widths = c(1, 1.5), align = "h", axis = "tb")

plot_grid(row1, row2, nrow = 2)

cellchat <- computeNetSimilarityPairwise(cellchat, slot.name = "netP", type = c("structural"))
cellchat <- netEmbedding(cellchat, slot.name = "netP", type = "structural", umap.method = "uwot", pathway.remove = NULL)
cellchat <- netClustering(cellchat, slot.name = "netP", type = "structural", do.parallel = FALSE, k = 6)

structural <- cellchat@netP$similarity$structural$dr$`1-2` %>% as.data.frame
structural$group <- cellchat@netP$similarity$structural$group$`1-2`
structural <- structural %>% rownames_to_column() %>% as_tibble() %>%
  separate(rowname, c("pathway","sample"), sep = "--")
ggplot(structural, aes(x = UMAP1, y = UMAP2, color = factor(group), shape = sample)) + geom_point() + geom_label_repel(aes(label = pathway)) + theme_bw()

netVisual_embeddingPairwise(cellchat, type = "structural", label.size = 3.5, top.label = 10)
netVisual_embeddingPairwiseZoomIn(cellchat, type = "structural", nCol = 2)

rankSimilarity(cellchat, type = "functional")

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2


i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 17)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 17)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 17, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 17, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 17, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 17, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

plot_grid(netVisual_bubble(cellchat, comparison = c(1, 2), angle.x = 45, signaling = "PDGF"), FeaturePlot(si.filt, "PDGFA", order = TRUE), nrow = 1, rel_widths = c(3,1))

netVisual_individual(cc.6hrpf, signaling = c("CNTN"))
                     
                     
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "fasted"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in fasted
net.up <- subsetCommunication(cellchat, net = net, datasets = "fasted",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in 6hrpf, i.e.,downregulated in fasted
net.down <- subsetCommunication(cellchat, net = net, datasets = "fasted",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("fasted", "x6hrpf")) # set factor level
plotGeneExpression(cellchat, signaling = "TENASCIN", split.by = "datasets", colors.ggplot = T)

#### GO THROUGH FIG2 TO GET CLUSTUMAP AND SI.META THEN RESUME HERE!!!
clust.cols <- c("#A2B627", "#F3A546", "#17BECF", "#A39FC9", "#F1788D", 
                         "#4E79A7", "#A0CBE8", "#59A14F", "#8CD17D", "#B6992D",
                         "#F1CE63", "#86BCB6", "#79706E", "#BAB0AC", "#D37295",
                         "#FABFD2", "#B07AA1", "#D4A6C8", "#9D7660", "#D7B5A6")
                         names(clust.cols) <- c(0:19)
                         
                         clust.labels <- data.frame(labels = names(clust.cols),
                                                    bg = clust.cols)
                         row.names(clust.labels) <- clust.labels$labels
                         
                         clust.labels$x <- NA
                         clust.labels$y <- NA
                         
dotsUMAP <- ggplot(si.meta, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters)) +
  rasterise(geom_point(size = 0.1, alpha = 0.5, shape = 20, color = "gray70"), dpi = 300, scale = 0.5) +
  theme_void() + coord_equal() +
  #xlab(paste0("UMAP1", sprintf('\u2192'))) + ylab(paste0("UMAP2", sprintf('\u2192'))) +
  scale_color_manual(values = clust.cols) +
  theme(legend.position = "none",
        plot.margin = margin(1,1,1,1),
        axis.title.x = element_text(size = 6, hjust = 0.01),
        axis.title.y = element_text(size = 6, hjust = 0.01, angle = 90)) +
  geom_point(data = clust.labels, aes(x = x, y = y, color = labels), size = 6, alpha = 0.5) +
  geom_text(data = clust.labels, aes(x = x, y = y, label = labels), color = "black", fontface = "bold", size = 2.5)

plot_grid(dotsUMAP, dotsUMAP, dotsUMAP, dotsUMAP, nrow = 1, labels = "e  BEST4 cell intracellular communication", label_size = 8, hjust = 0, label_x = 0.05)

#par(mfrow = c(1,4), xpd=TRUE)
dev.off()
b400h <- netVisual_chord_gene(cc.fasted, sources.use = c("BEST4+/SPIB+"), slot.name = "netP")
b406h <- netVisual_chord_gene(cc.6hrpf, sources.use = c("BEST4+/SPIB+"), slot.name = "netP")
b412h <- netVisual_chord_gene(cc.12hrpf, sources.use = c("BEST4+/SPIB+"), slot.name = "netP")
b424h <- netVisual_chord_gene(cc.1dpf, sources.use = c("BEST4+/SPIB+"), slot.name = "netP")

netVisual_bubble(cc.fasted, sources.use = "BEST4+/SPIB+", remove.isolate = TRUE)
netAnalysis_signalingRole_network(cc.12hrpf)


gg00 <- netAnalysis_signalingRole_scatter(cc.fasted)
gg06 <- netAnalysis_signalingRole_scatter(cc.6hrpf)
gg12 <- netAnalysis_signalingRole_scatter(cc.12hrpf)
gg24 <- netAnalysis_signalingRole_scatter(cc.1dpf)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cc.fasted, signaling = c("FGF", "NOTCH"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg00 + gg06 + gg12 + gg24 + plot_layout(nrow = 1)

par(mfrow = c(2,4), xpd=TRUE)
netVisual_aggregate(cc.fasted, signaling = "FGF", layout = "chord")
netVisual_aggregate(cc.6hrpf, signaling = "FGF", layout = "chord")
netVisual_aggregate(cc.12hrpf, signaling = "FGF", layout = "chord")
netVisual_aggregate(cc.1dpf, signaling = "FGF", layout = "chord")

#par(mfrow = c(2,2), xpd=TRUE)
#netVisual_aggregate(cc.fasted, signaling = "NOTCH", layout = "chord")
#netVisual_aggregate(cc.6hrpf, signaling = "NOTCH", layout = "chord")
#netVisual_aggregate(cc.12hrpf, signaling = "NOTCH", layout = "chord")
#netVisual_aggregate(cc.1dpf, signaling = "NOTCH", layout = "chord")

netAnalysis_signalingRole_network(cc.fasted, signaling = c("NOTCH"))
netAnalysis_signalingRole_network(cc.6hrpf, signaling = c("NOTCH"))
netAnalysis_signalingRole_network(cc.12hrpf, signaling = c("NOTCH"))
netAnalysis_signalingRole_network(cc.1dpf, signaling = c("NOTCH"))

netAnalysis_signalingRole_network(cc.fasted, signaling = c("FGF"))
netAnalysis_signalingRole_network(cc.6hrpf, signaling = c("FGF"))
netAnalysis_signalingRole_network(cc.12hrpf, signaling = c("FGF"))
netAnalysis_signalingRole_network(cc.1dpf, signaling = c("FGF"))



fgfFasted <- netAnalysis_contribution(cc.fasted, signaling = c("FGF"))
fgfFasted <- fgfFasted$data %>% mutate(time = "fasted")

notchFasted <- netAnalysis_contribution(cc.fasted, signaling = c("NOTCH"))
notchFasted <- notchFasted$data %>% mutate(time = "fasted")

fgf6hr <- netAnalysis_contribution(cc.6hrpf, signaling = c("FGF"))
fgf6hr <- fgf6hr$data %>% mutate(time = "6hrpf")

notch6hr <- netAnalysis_contribution(cc.6hrpf, signaling = c("NOTCH"))
notch6hr <- notch6hr$data %>% mutate(time = "6hrpf")

fgf12hr <- netAnalysis_contribution(cc.12hrpf, signaling = c("FGF"))
fgf12hr <- fgf12hr$data %>% mutate(time = "12hrpf")

notch12hr <- netAnalysis_contribution(cc.12hrpf, signaling = c("NOTCH"))
notch12hr <- notch12hr$data %>% mutate(time = "12hrpf")

fgf1d <- netAnalysis_contribution(cc.1dpf, signaling = c("FGF"))
fgf1d <- fgf1d$data %>% mutate(time = "1dpf")

notch1d <- netAnalysis_contribution(cc.1dpf, signaling = c("NOTCH")) 
notch1d <- notch1d$data %>% mutate(time = "1dpf")

fgfLR <- rbind(fgfFasted, fgf6hr) %>% rbind(fgf12hr) %>% rbind(fgf1d) %>%
  mutate(time = factor(time, levels = c("fasted", "6hrpf", "12hrpf", "1dpf")))

notchLR <- rbind(notchFasted, notch6hr) %>% rbind(notch12hr) %>% rbind(notch1d) %>%
  mutate(time = factor(time, levels = c("fasted", "6hrpf", "12hrpf", "1dpf"))) %>%
  filter(contribution > 0)

fgfLR_plot <- fgfLR %>% ggplot(aes(x = reorder_within(name, contribution, time), y = contribution)) +
  geom_col() +
  ggtitle("FGF L-R pairs") +
  facet_grid(rows = vars(time), scales = "free", space = "free") + coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = expansion(add = 0.01), name = "relative contribution") +
  scale_x_discrete(labels = fgfLR$name) +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", size = 8),
        text = element_text(size = 7),
        axis.title.x = element_text(size = 7, face = "bold"))

notchLR_plot <- notchLR %>% ggplot(aes(x = reorder_within(name, contribution, time), y = contribution)) +
  geom_col() +
  ggtitle("NOTCH L-R pairs") +
  facet_grid(rows = vars(time), scales = "free", space = "free") + coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = expansion(add = 0.01), name = "relative contribution") +
  scale_x_discrete(labels = notchLR$name) +
  theme(axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", size = 8),
        text = element_text(size = 7),
        axis.title.x = element_text(size = 7, face = "bold"))


##### recreating a hierarchy plot because. the CellChat function is busted
cells.level <- levels(cc.fasted@idents)
color.use <- scPalette(nrow(net))
names(color.use) <- cells.level


getSignalProbs <- function(object, pathway, thresh, time) {
  object = object
  signaling = pathway
  thresh = thresh
  time = time
  
  pairLR <- searchPair(signaling = signaling, pairLR.use = object@LR$LRsig, key = "pathway_name", matching.exact = T, pair.only = T)
  net <- object@net
  
  pairLR.use.name <- dimnames(net$prob)[[3]]
  pairLR.name <- intersect(rownames(pairLR), pairLR.use.name)
  pairLR <- pairLR[pairLR.name, ]
  prob <- net$prob
  pval <- net$pval
  
  prob[pval > thresh] <- 0
  
  if (length(pairLR.name) > 1) { pairLR.name.use <- pairLR.name[apply(prob[,,pairLR.name], 3, sum) != 0] } else { pairLR.name.use <- pairLR.name[sum(prob[,,pairLR.name]) != 0] }
  
  pairLR <- pairLR[pairLR.name.use,]
  prob <- prob[,,pairLR.name.use]
  pval <- pval[,,pairLR.name.use]
  
  prob.sum <- apply(prob, c(1,2), sum)
  
  df.net <- reshape2::melt(prob.sum, value.name = "value")
  colnames(df.net)[1:2] <- c("source","target")
  
  df.net$source <- factor(df.net$source, levels = cells.level)
  df.net$target <- factor(df.net$target, levels = cells.level)
  df.net$value[is.na(df.net$value)] <- 0
  net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  
  df.net <- df.net %>% mutate(start = "source",
                              end = "target",
                              source.size = as.vector(table(object@idents)[source]),
                              target.size = as.vector(table(object@idents)[target]),
                              source = factor(source, levels = rev(cells.level)),
                              target = factor(target, levels = rev(cells.level)),
                              time = time) 
  
  return(df.net)
}

df.Fas <- getSignalProbs(cc.fasted, "NOTCH", 0.05, "fasted")
df.6hr <- getSignalProbs(cc.6hrpf, "NOTCH", 0.05, "6hrpf")
df.12hr <- getSignalProbs(cc.12hrpf, "NOTCH", 0.05, "12hrpf")
df.24hr <- getSignalProbs(cc.1dpf, "NOTCH", 0.05, "1dpf") %>% filter(source != "other mesenchyme", target != "other mesenchyme")

df.notch <- rbind(df.Fas, df.6hr) %>% rbind(df.12hr) %>% rbind(df.24hr)

hierarch <- df.notch  %>% mutate(time = factor(time, levels = c("fasted", "6hrpf", "12hrpf", "1dpf"))) %>%
  ggplot(aes(x = start, y = source, xend = end, yend = target)) +
  geom_segment(aes(linewidth = value, color = source), alpha = 0.5) +
  geom_point(aes(x = start, y = source, color = source, size = source.size), shape = 16) +
  geom_point(aes(x = end, y = target, color = target, size = target.size), shape = 21, fill = "white") +
  scale_linewidth_continuous(range = c(0, 3)) +
  scale_radius() +
  scale_color_manual(values = color.use) + 
  scale_x_discrete(expand = expansion(add = 0.2), position = "top") +
  scale_y_discrete(drop = FALSE) +
  facet_grid(cols = vars(time)) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin = margin(3, 3, 3, 3),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6, hjust = 1),
        strip.placement = "outside")

bars <- notchLR %>% ggplot(aes(x = reorder_within(name, -contribution, time), y = contribution)) +
  geom_col() +
  ggtitle("NOTCH L-R pairs") +
  facet_grid(cols = vars(time), scales = "free", space = "free") +
  theme_classic() +
  scale_y_continuous(expand = expansion(add = 0.01), name = "relative contribution") +
  #scale_x_discrete(labels = notchLR$name) +
  theme(axis.title.x = element_blank(),
        strip.background =  element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(face = "bold", size = 8),
        text = element_text(size = 7),
        axis.line = element_line(size = 0.5),
        axis.title.y = element_text(size = 7, face = "bold"),
        axis.text.x = element_text(angle = 90))

ggplot_build(bars)$layout$panel_params[[1]]$x$get_labels() %>% gsub(pattern = "___.*", replacement = "", perl = TRUE)
ggplot_build(bars)$layout$panel_params[[2]]$x$get_labels() %>% gsub(pattern = "___.*", replacement = "", perl = TRUE)
ggplot_build(bars)$layout$panel_params[[3]]$x$get_labels() %>% gsub(pattern = "___.*", replacement = "", perl = TRUE)
ggplot_build(bars)$layout$panel_params[[4]]$x$get_labels() %>% gsub(pattern = "___.*", replacement = "", perl = TRUE)


fgfBars <- fgfLR %>% ggplot(aes(x = reorder_within(name, -contribution, time), y = contribution)) +
  geom_col() +
  ggtitle("FGF L-R pairs") +
  facet_grid(cols = vars(time), scales = "free", space = "free") + 
  theme_classic() +
  scale_y_continuous(expand = expansion(add = 0.01), name = "relative contribution") +
  #scale_x_discrete(labels = notchLR$name) +
  theme(axis.title.x = element_blank(),
        strip.background =  element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(face = "bold", size = 8),
        text = element_text(size = 7),
        axis.line = element_line(size = 0.5),
        axis.title.y = element_text(size = 7, face = "bold"),
        axis.text.x = element_text(angle = 90))

netVisual_chord_gene(cc.1dpf, sources.use = c("BEST4+/SPIB+", "OLFM4+/SPDEF+"), signaling = "FGF", lab.cex = 0.5,legend.pos.y = 30)


######
## Cluster 1
cl1_6hr <- c("COLLAGEN", "LAMININ", "THBS", "FN1",
             "FGF", 'NOTCH', "EPHB", "PDGF", "IGF",
             "GAS", "L1CAM", "TGFb", "CSF", "ACTIVIN",
             "GDNF", "KIT", "NT", "EPO", "CNTN",
             "NGF", "FLT3", "BAFF")
cl16hr_plot <- netAnalysis_signalingRole_scatter(cc.6hrpf, signaling = cl1_6hr, do.label = F, dot.size = c(1,4)) +
  ggtitle("cluster 1 (6hrpf)") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

## Cluster 2
cl2_fasted <- c("PECAM1","CD45","ANGPTL","NCAM","PROS","SEMA5","GDNF","ESAM")
cl2_6hr <- c("PECAM1","CD45","NCAM","APP","GCG","SEMA7","CSPG4","SEMA5","PERIOSTIN",
             "GRN","ESAM","NRXN","AGRN","NECTIN","ANGPTL","CD34","CD99","SELPLG","OCLN")

cl2fas_plot <- netAnalysis_signalingRole_scatter(cc.fasted, signaling = cl2_fasted, do.label = F, dot.size = c(1,4))+
  ggtitle("cluster 2 (fasted)") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
cl26hr_plot <- netAnalysis_signalingRole_scatter(cc.6hrpf, signaling = cl2_6hr, do.label = F, dot.size = c(1,4))+
  ggtitle("cluster 2 (6hrpf)") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))


## Cluster 3
cl3_fasted <- c("CDH1","EPHA","CDH","CADM","GRN","ACTIVIN","FLT3","CSF",
                "EPO", "NT", "OCLN", "BAFF", "COMPLEMENT")
cl3_6hr <- c("EPHA","CADM")

cl3fas_plot <- netAnalysis_signalingRole_scatter(cc.fasted, signaling = cl3_fasted, do.label = F, dot.size = c(1,4))+
  ggtitle("cluster 3 (fasted)") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
cl36hr_plot <- netAnalysis_signalingRole_scatter(cc.6hrpf, signaling = cl3_6hr, do.label = F, dot.size = c(1,4))+
  ggtitle("cluster 3 (6hrpf)") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

## Cluster 4
cl4_fasted <- c("NRG","GCG","PERIOSTIN","BMP","VWF","KIT","VTN","NGF")
cl4_6hr <- c("TENASCIN","NRG","EGF","VWF","CDH1","RELN","BMP")
cl4fas_plot <- netAnalysis_signalingRole_scatter(cc.fasted, signaling = cl4_fasted, do.label = F, dot.size = c(1,4))+
  ggtitle("cluster 4 (fasted)") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
cl46hr_plot <- netAnalysis_signalingRole_scatter(cc.6hrpf, signaling = cl4_6hr, do.label = F, dot.size = c(1,4))+
  ggtitle("cluster 4 (6hrpf)") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

## Cluster 5
cl5_fasted <- c("COLLAGEN","LAMININ","PDGF","FN1","THBS","IGF","L1CAM","VEGF","CDH5","TGFb")
cl5_6hr <- c("VEGF","CDH5","SEMA3","SEMA6","CDH","PARs","VTN")
cl5fas_plot <- netAnalysis_signalingRole_scatter(cc.fasted, signaling = cl5_fasted, do.label = F, dot.size = c(1,4))+
  ggtitle("cluster 5 (fasted)") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
cl56hr_plot <- netAnalysis_signalingRole_scatter(cc.6hrpf, signaling = cl5_6hr, do.label = F, dot.size = c(1,4))+
  ggtitle("cluster 5 (6hrpf)") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

## Cluster 6
cl6_fasted <- c("FGF","NOTCH","LIGHT","SEMA6","JAM","EPHB","SEMA3","GAS","ANGPT")
cl6_6hr <- c("HSPG","BTLA","ALCAM","CD6","LIGHT","ANGPT","VISFATIN")
cl6fas_plot <- netAnalysis_signalingRole_scatter(cc.fasted, signaling = cl6_fasted, do.label = F, dot.size = c(1,4))+
  xlab("outgoing interaction strength") + ylab("incoming interaction strength") +
  ggtitle("cluster 6 (6hrpf)") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))
cl66hr_plot <- netAnalysis_signalingRole_scatter(cc.6hrpf, signaling = cl6_6hr, do.label = F, dot.size = c(1,4))+
  ggtitle("cluster 6 (6hrpf)") +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 7),
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

plot_grid(NULL, cl16hr_plot, cl2fas_plot, cl26hr_plot, cl3fas_plot, cl36hr_plot,
          cl4fas_plot, cl46hr_plot, cl5fas_plot, cl56hr_plot, cl6fas_plot, cl66hr_plot, 
          nrow = 2)

plot_grid(NULL, cl2fas_plot, cl3fas_plot, cl4fas_plot, cl5fas_plot, cl6fas_plot, 
          cl16hr_plot, cl26hr_plot, cl36hr_plot, cl46hr_plot, cl56hr_plot, cl66hr_plot, 
          nrow = 2, align = "hv", axis = "tblr")


netVisual_heatmap(cc.fasted, signaling = c("VEGF"), remove.isolate = FALSE, color.heatmap = "Reds")
netVisual_heatmap(cc.6hrpf, signaling = c("VEGF"), remove.isolate = FALSE, color.heatmap = "Reds")
netVisual_heatmap(cc.12hrpf, signaling = c("VEGF"), remove.isolate = FALSE, color.heatmap = "Reds")
netVisual_heatmap(cc.1dpf, signaling = c("VEGF"), remove.isolate = FALSE, color.heatmap = "Reds")
netVisual_bubble(cc.6hrpf, sources.use = "BEST4+/SPIB+", signaling = c("FGF","NOTCH","OCLN"), remove.isolate = FALSE)
netVisual_bubble(cc.12hrpf, sources.use = "BEST4+/SPIB+", signaling = c("FGF","NOTCH","OCLN"), remove.isolate = FALSE)
netVisual_bubble(cc.1dpf, sources.use = "BEST4+/SPIB+", signaling = c("FGF","NOTCH","OCLN"), remove.isolate = FALSE)

plot_grid(
netAnalysis_signalingRole_scatter(cc.fasted, signaling = c("FGF")) + ggtitle("a. Fasted FGF signaling") + theme(legend.position = "none", plot.title = element_text(hjust = 0, face = "bold")),
netAnalysis_signalingRole_scatter(cc.6hrpf, signaling = c("FGF")) + ggtitle("b. 6hrpf FGF signaling") + theme(legend.position = "none", plot.title = element_text(hjust = 0, face = "bold")),
netAnalysis_signalingRole_scatter(cc.12hrpf, signaling = c("FGF")) + ggtitle("c. 12hrpf FGF signaling") + theme(legend.position = "none", plot.title = element_text(hjust = 0, face = "bold")),
netAnalysis_signalingRole_scatter(cc.1dpf, signaling = c("FGF")) + ggtitle("d. 1dpf FGF signaling") + theme(legend.position = "none", plot.title = element_text(hjust = 0, face = "bold")),
netAnalysis_signalingRole_scatter(cc.fasted, signaling = c("NOTCH")) + ggtitle("e. Fasted NOTCH signaling") + theme(legend.position = "none", plot.title = element_text(hjust = 0, face = "bold")),
netAnalysis_signalingRole_scatter(cc.6hrpf, signaling = c("NOTCH")) + ggtitle("f. 6hrpf NOTCH signaling") + theme(legend.position = "none", plot.title = element_text(hjust = 0, face = "bold")),
netAnalysis_signalingRole_scatter(cc.12hrpf, signaling = c("NOTCH")) + ggtitle("g. 12hrpf NOTCH signaling") + theme(legend.position = "none", plot.title = element_text(hjust = 0, face = "bold")),
netAnalysis_signalingRole_scatter(cc.1dpf, signaling = c("NOTCH")) + ggtitle("h. 1dpf NOTCH signaling") + theme(legend.position = "none", plot.title = element_text(hjust = 0, face = "bold")),
ncol = 2)

ht1 <- netAnalysis_signalingRole_heatmap(cc.fasted, pattern = "outgoing", height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(cc.fasted, pattern = "incoming", height = 15)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cc.6hrpf, pattern = "outgoing", height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(cc.6hrpf, pattern = "incoming", height = 15)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cc.12hrpf, pattern = "outgoing", height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(cc.12hrpf, pattern = "incoming", height = 15)
ht1 + ht2

ht1 <- netAnalysis_signalingRole_heatmap(cc.1dpf, pattern = "outgoing", height = 15)
ht2 <- netAnalysis_signalingRole_heatmap(cc.1dpf, pattern = "incoming", height = 15)
ht1 + ht2

### Crypt genes
FeaturePlot(si.filt, features = c("MKI67", "CCNB1", "CCND1", "MCM2", "PCNA", "OLFM4"), order = TRUE)

### Bottom villus genes
FeaturePlot(si.filt, features = c("NLRP6", "LYPD8", "IL18", "REG1", "REG3A"))

### Mid-villus genes
FeaturePlot(si.filt, features = c("SLC5A1", "SLC2A5", "SLC2A2", "SLC7A7", "SLC7A8", "SLC7A9"), order = TRUE)

### Villus tip genes
FeaturePlot(si.filt, features = c("ADA", "NT5E", "SLC28A2", "CREB3L3", "APOA1", "APOB"), order = TRUE)

FeaturePlot(si.filt, features = c("MKI67", "SLC7A7"), blend = TRUE, order = TRUE)
FeaturePlot(si.filt, features = c("CCND1", "RAB18"), blend = TRUE, order = TRUE)




