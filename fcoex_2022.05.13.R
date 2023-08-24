library(fcoex)
library(DESeq2)
library(Seurat)

setwd("~/Desktop/snRNAseq_python")

## Whole sample:
si.filt <- readRDS("reanalysis.RDS")  
sce <- as.SingleCellExperiment(si.filt, assay = "RNA")

## One sample at a time:
sub <- subset(si.filt, subset = orig.ident == "AR14")
sce <- as.SingleCellExperiment(sub, assay = "RNA")

target <- colData(sce)
target <- target$seurat_clusters

# Get normalized table from the pre-processing
exprs <- as.data.frame(assay(sce, 'logcounts'))

# Create fcoex object
fc <- new_fcoex(data.frame(exprs),target)

fc <- discretize(fc, number_of_bins = 16)

fc <- find_cbf_modules(fc,n_genes_selected_in_first_step = 200, verbose = FALSE, is_parallel = TRUE)
fc <- get_nets(fc)
## AQ16: 19 modules
## AR4: 11 modules
## AP22: 
## AR14: 11 modules

saveRDS(fc, file = "fc_AR14_08.03.RDS")
fc <- readRDS("fc_allCells_07.25.RDS")
fc.fasted <- readRDS("fc_AQ16_07.25.RDS")
fc.6hrpf <- readRDS("fc_AR4_07.25.RDS")
fc.12hrpf <- readRDS("fc_AP22_07.25.RDS")

network_plots <- show_net(fc)
### ALL:


### BEST4 cells:
# 14 --> NUPR1 (increases stemness, critical regulator of antioxidants)
# 16 --> LOC103064373, uncharacterized loci

### 11: CALCRL
##### LOC103050314 MEIS1-like
##### LOC103060783 GYPC
##### LOC103068020 rho GTPase-activating protein 20-like
##### LOC103057419 AKAP2/PALM2

### AQ16 (fasted): ######
## 1. FLNA
## 2. FNBP1
## 3. SLC26A9
## 4. MEIS2-like
## 5. SYNE1 NULL
## 6. DYSF
## 7. ADGRL1
## 8. ABCC4 NULL
## 9. TRPM6
## 10. LOC112541228 uncharacterized ncRNA -- cell adhesion, breast cancer
## 11. LOC103064376 = DOCK2
## 12. ZEB2
## 13. AKT3
## 14. PBX3
## 15. CCDC88A
## 16. ZNF521
## 17. MAGI2
## 18. FAM129A NULL
## 19. LOC103058238 NULL

### AP22 (12hrpf): ######
## 1. LOC103061637 = MFSD2A-like
## 2. SPTBN5 NULL
## 3. FNBP1
## 4. FLNA
## 5. LOC103058653 = CLEC2D-like
## 6. EGFL7
## 7. ZEB2
## 8. PLTP
## 9. LOC103066265 NULL
## 10. CDHR2 NULL
## 11. PTPRC NULL
## 12. TOX
## 13. LOC103062401 NULL

### AR14 (1dpf): ######
## 1. TSN1
## 2. LOC103054953 (serine protease 27-like)/PRSS27-like
## 3. FNBP1
## 4. EGFL7 (inhibits PDFG-BB-induced smooth muscle cell migration and promotes adhesion and angiogenesis)
## 5. ZEB2 (key role in TGFbeta signaling during early development, interacts with SMADs; destabilizes intestinal barrier --> immune infiltration --> cancer)
## 6. COL6A3
## 7.  KCNMA1
## 8. FAP
## 9. LOC103067860 null
## 10. LOC112540173 null
## 11. MEIS1 null
######

### AQ16 (fasted): ######
######

plot_grid(network_plots[[1]], network_plots[[3]],
             network_plots[[6]], network_plots[[7]],
             network_plots[[13]], network_plots[[17]])

a <- network_plots[[1]]
b <- network_plots[[2]]
c <- network_plots[[6]]
d <- network_plots[[7]]
e <- network_plots[[13]]
f <- network_plots[[16]]

## network 6: CSF1R:
#### LOC103048864 - Factor H
#### LOC103057442 - DMBT1-like (development, angiogenesis)
#### LOC103055338 - ALOX5-like (inflammation, immune response)
#### LOC103049533 - STAB1-like (angiogenesis)
#### LOC103064031 - NOX2 (forms antioxidants)
#### LOC103058156 - STAB1
#### LOC103064043 - CPA3 -- mucosal mast cell subtype

## LOC103056554 == MEIS2-like

plot_grid(a, b, c, d, e, f, nrow = 3)

a <- network_plots[[1]]
b <- network_plots[[2]]
c <- network_plots[[3]]
d <- network_plots[[4]]
e <- network_plots[[5]]
f <- network_plots[[6]]
k <- network_plots[[11]]
l <- network_plots[[12]]

nets3 <- plot_grid(c, f, e, nrow = 1)
nets2 <- plot_grid(d, k, l, nrow = 3)

netsA <- plot_grid(b, nets2, rel_widths = c(1, 0.5))
plot_grid(netsA, nets3, rel_heights = c(1, 0.5), nrow = 2)

r1c2 <- plot_grid(c, d, ncol = 1)
r1 <- plot_grid(b, r1c2, nrow = 1, rel_widths = c(2, 1))

r2c1 <- plot_grid(e, f, ncol = 1)
r2 <- plot_grid(r2c1, a, nrow = 1, rel_heights = c(1, 2.5))

plot_grid(r1, r2, nrow = 2)

fc <- recluster(fc)

gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- pathwayPCA::read_gmt(gmt_fname)
#fc <- mod_ora(fc, gmt_in)

si.fc <- si.filt
Idents(si.fc) <- idents(fc)[[1]]

for (i in names(module_genes(fc))) {
  Idents(si.fc) <-   fc@mod_idents[[i]]
  
  print(paste("Checking for anticorrelation in module", i))
  
  # Identify markers only for module genes:
  module_genes_in_clusters <-
    FindAllMarkers(
      si.fc,
      logfc.threshold = 1,
      only.pos = TRUE,
      features = fc@module_list[[i]],
      verbose = FALSE
    )
  
  # If there are markers of the HN cluster, it means that we have anticorrelation
  if ("HN" %in% module_genes_in_clusters$cluster) {
    module_genes_in_clusters$module = i
    message(paste0("anticorrelated genes found for module ", i))
  }
}

Idents(si.fc) <-   fc@mod_idents[["APOA4"]]

markers_for_cluster_APOA4 <-
  FindAllMarkers(
    si.fc,
    logfc.threshold = 1,
    only.pos = TRUE,
    features = fc@module_list[["APOA4"]],
    verbose = FALSE
  )

head(markers_for_cluster_APOA4)
APOA4.HN <- c("OTOP3", "KIAA1211", "LOC103048842", "SLC43A2", "APOA4")
FeaturePlot(si.filt, APOA4.HN, order = TRUE)

OTOP3 <- FeaturePlot(si.filt, features = "OTOP3", order = TRUE) +
  ggtitle("OTOP3") +
  scale_color_gradient(low = "#d3d3d3", high = "#CE6091") +
  theme_void() + coord_cartesian(expand = FALSE) +
  xlab(paste0("UMAP1", sprintf('\u2192'))) + ylab(paste0("UMAP2", sprintf('\u2192'))) +
  guides(colour = guide_colorbar(barheight = 2, barwidth = 0.5, title.position = "top", direction = "vertical")) +
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6, hjust = 0.01),
        axis.title.y = element_text(size = 6, hjust = 0.01, angle = 90),
        legend.position = "right")
OTOP3$layers[[1]]$aes_params$size <- 0.01
OTOP3$layers[[1]]$aes_params$alpha <- 0.2
OTOP3$layers[[1]]$aes_params$shape <- 20

KIAA1211 <- FeaturePlot(si.filt, features = "KIAA1211", order = TRUE) +
  ggtitle("KIAA1211") +
  scale_color_gradient(low = "#d3d3d3", high = "#CE6091") +
  theme_void() + coord_cartesian(expand = FALSE) +
  xlab(paste0("UMAP1", sprintf('\u2192'))) + ylab(paste0("UMAP2", sprintf('\u2192'))) +
  guides(colour = guide_colorbar(barheight = 2, barwidth = 0.5, title.position = "top", direction = "vertical")) +
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6, hjust = 0.01),
        axis.title.y = element_text(size = 6, hjust = 0.01, angle = 90),
        legend.position = "right")
KIAA1211$layers[[1]]$aes_params$size <- 0.01
KIAA1211$layers[[1]]$aes_params$alpha <- 0.4
KIAA1211$layers[[1]]$aes_params$shape <- 20

LOC103048842 <- FeaturePlot(si.filt, features = "LOC103048842", order = TRUE) +
  ggtitle("SLC7A9") +
  scale_color_gradient(low = "#d3d3d3", high = "#CE6091") +
  theme_void() + coord_cartesian(expand = FALSE) +
  xlab(paste0("UMAP1", sprintf('\u2192'))) + ylab(paste0("UMAP2", sprintf('\u2192'))) +
  guides(colour = guide_colorbar(barheight = 2, barwidth = 0.5, title.position = "top", direction = "vertical")) +
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6, hjust = 0.01),
        axis.title.y = element_text(size = 6, hjust = 0.01, angle = 90),
        legend.position = "right")
LOC103048842$layers[[1]]$aes_params$size <- 0.01
LOC103048842$layers[[1]]$aes_params$alpha <- 0.4
LOC103048842$layers[[1]]$aes_params$shape <- 20

SLC43A2 <- FeaturePlot(si.filt, features = "SLC43A2", order = TRUE) +
  ggtitle("SLC43A2") +
  scale_color_gradient(low = "#d3d3d3", high = "#CE6091") +
  theme_void() + coord_cartesian(expand = FALSE) +
  xlab(paste0("UMAP1", sprintf('\u2192'))) + ylab(paste0("UMAP2", sprintf('\u2192'))) +
  guides(colour = guide_colorbar(barheight = 2, barwidth = 0.5, title.position = "top", direction = "vertical")) +
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6, hjust = 0.01),
        axis.title.y = element_text(size = 6, hjust = 0.01, angle = 90),
        legend.position = "right")
SLC43A2$layers[[1]]$aes_params$size <- 0.01
SLC43A2$layers[[1]]$aes_params$alpha <- 0.4
SLC43A2$layers[[1]]$aes_params$shape <- 20


apoa4.mod <- UMAPPlot(si.fc) +
  scale_color_manual(values = c("#d3d3d3", "#CE6091"), labels = c("negative", "positive")) + 
  ggtitle("APOA4 module expression") +
  theme_void() + coord_cartesian(expand = FALSE) +
  xlab(paste0("UMAP1", sprintf('\u2192'))) + ylab(paste0("UMAP2", sprintf('\u2192'))) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2),
                               keywidth = 0.5, direction = "horizontal")) +
  theme(text = element_text(size = 6),
        axis.title.x = element_text(size = 6, hjust = 0.01),
        axis.title.y = element_text(size = 6, hjust = 0.01, angle = 90),
        legend.position = "bottom")
apoa4.mod$layers[[1]]$aes_params$size <- 0.01
apoa4.mod$layers[[1]]$aes_params$alpha <- 0.5

umap.sample + ggtitle("Sample") + scale_color_manual(values = plasma(4), labels = c("fasted", "6hrpf", "12hrpf", "1dpf")) + guides(color = guide_legend(nrow = 1, override.aes = list(size = 2)))

col2 <- plot_grid(apoa4.mod,
                  umap.sample + ggtitle("Time point") + scale_color_manual(values = plasma(4), labels = c("fasted", "6hrpf", "12hrpf", "1dpf")) + guides(color = guide_legend(nrow = 1, override.aes = list(size = 2))),
                  ncol = 1)

title_gg <- ggplot() + 
  labs(title = "Negatively correlated genes") + 
  umap.mini_theme
col3 <- plot_grid(title_gg, OTOP3, KIAA1211, LOC103048842, SLC43A2, ncol = 1,
                  rel_heights = c(0.1, 1, 1, 1, 1))

plot_grid(col1, col2, col3, ncol = 3)

fc.fasted <- readRDS("fc_AQ16_07.25.RDS")
network_plots.fasted <- show_net(fc.fasted)
fastedA <- network_plots.fasted[[1]]
fastedB <- network_plots.fasted[[2]]
fastedC <- network_plots.fasted[[3]]
fastedD <- network_plots.fasted[[4]]
fastedE <- network_plots.fasted[[6]]
fastedF <- network_plots.fasted[[7]]
fastedG <- network_plots.fasted[[9]]
fastedH <- network_plots.fasted[[10]]
fastedI <- network_plots.fasted[[11]]
fastedJ <- network_plots.fasted[[12]]
fastedK <- network_plots.fasted[[13]]
fastedL <- network_plots.fasted[[14]]
fastedM <- network_plots.fasted[[15]]
fastedN  <- network_plots.fasted[[16]]

plot_grid(fastedA, fastedB, fastedC,
          fastedD, fastedE, fastedF,
          fastedG, fastedH, fastedI,
          fastedJ, fastedK, fastedL,
          fastedM, fastedN)

fc.6hrpf <- readRDS("fc_AR4_07.25.RDS")
network_plots.6hrpf <- show_net(fc.6hrpf)
a6hrpfA <- network_plots.6hrpf[[1]]
a6hrpfB <- network_plots.6hrpf[[2]]
a6hrpfC <- network_plots.6hrpf[[3]]
a6hrpfD <- network_plots.6hrpf[[4]]
a6hrpfE <- network_plots.6hrpf[[6]]
a6hrpfF <- network_plots.6hrpf[[7]]
a6hrpfG <- network_plots.6hrpf[[9]]
a6hrpfH <- network_plots.6hrpf[[10]]
a6hrpfI <- network_plots.6hrpf[[11]]

plot_grid(a6hrpfA, a6hrpfB, a6hrpfC,
          a6hrpfD, a6hrpfE, a6hrpfF,
          a6hrpfG, a6hrpfH, a6hrpfI)

fc.12hrpf <- readRDS("fc_AP22_07.25.RDS")
network_plots.12hrpf <- show_net(fc.12hrpf)

FeaturePlot(si.filt, c("LOC112541014", "LOC112540815"), order = TRUE)

LOC112541014 <- FeaturePlot(si.filt, features = "LOC112541014", order = TRUE) +
  ggtitle("LOC112541014") +
  scale_color_gradient(low = "#d3d3d3", high = "#C44326") + umap.mini_theme +
  xlab("UMAP 1") + ylab("UMAP 2") +
  coord_cartesian(xlim=c(-10, 14), ylim=c(-16, 9), expand = FALSE, clip = "off") +
  annotate(x=-10, xend=-10, y=-16, yend=-10, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  annotate(x=-10, xend=-3, y=-16, yend=-16, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt')))
LOC112541014$layers[[1]]$aes_params$size <- 0.01
LOC112541014$layers[[1]]$aes_params$alpha <- 0.9
LOC112541014$layers[[1]]$aes_params$shape <- 20

LOC112540815 <- FeaturePlot(si.filt, features = "LOC112540815", order = TRUE) +
  ggtitle("LOC112540815") +
  scale_color_gradient(low = "#d3d3d3", high = "#C44326") + umap.mini_theme +
  xlab("UMAP 1") + ylab("UMAP 2") +
  coord_cartesian(xlim=c(-10, 14), ylim=c(-16, 9), expand = FALSE, clip = "off") +
  annotate(x=-10, xend=-10, y=-16, yend=-10, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt'))) +
  annotate(x=-10, xend=-3, y=-16, yend=-16, colour="black", lwd=0.25, geom="segment",
           arrow = arrow(type='closed', length = unit(4,'pt')))
LOC112540815$layers[[1]]$aes_params$size <- 0.01
LOC112540815$layers[[1]]$aes_params$alpha <- 0.9
LOC112540815$layers[[1]]$aes_params$shape <- 20

grid.arrange(LOC112540815, LOC112541014, nrow = 1)
