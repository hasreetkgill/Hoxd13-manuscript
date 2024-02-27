# scRNA-seq analysis for project: 'Hox genes modulate physical forces..."
# Hasreet Gill, November 2023

# Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(glmGamPoi)
library(reticulate)
library(leidenAlg)
library(igraph)
library(SingleCellExperiment)
library(scater)
library(cowplot)
library(edgeR)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(pheatmap)
library(apeglm)
library(png)
library(DESeq2)
library(RColorBrewer)
library(scran)
library(ggrepel)
library(SeuratData)
# library(SeuratWrappers)
# library(Azimuth)
library(ggplot2)
library(scCustomize)
library(viridis)
library("wesanderson")
options(future.globals.maxSize = 1e9)
options(Seurat.object.assay.version = "v5")

# Set folder containing folders with pre-processed Seurat files as working directory
setwd("C:/Users/lalig/OneDrive/Documents/Tabin Lab/Hoxd13 project/scRNA-seq/hgl5")

# Load pre-processed Seurat objects & combine into one singlecell object for filtering/processing
# singlecell.filter <- readRDS( file = "midhind_filter_processed_annotated.rds" ) # Read in saved filtered Seurat objects
# hgl6 <- readRDS(file = "HGL6/20231030_HGL_reannotated_HGL6.rds") # pPIDNB midgut control
# hgl7 <- readRDS(file = "HGL7/20231030_HGL_reannotated_HGL7.rds") # pPINDB/pPIDNB-cHoxd13co hindgut
# hgl8 <- readRDS(file = "HGL8/20231030_HGL_reannotated_HGL8.rds") # pPIDNB-cHoxd13co midgut
# hgl5m <- readRDS(file = "HGL8/20231030_HGL_reannotated_HGL8.rds") # e14 midgut
# hgl5h <- readRDS(file = "HGL8/20231030_HGL_reannotated_HGL8.rds") # e14 midgut

# Figure 4 and S5
singlecell.filter <- tribble(
  ~seurat, ~name,
  hgl6, "HGL6",
  hgl7, "HGL7",
  hgl8, "HGL8",
)

# Figure S7
singlecell.filter <- tribble(
  ~seurat, ~name,
  hgl5m, "HGL5M",
  hgl5h, "HGL5H",
)


## Filter cells NOTE: this has already been performed for publicly available Seurat objects!
# Remove cells with high mitochondrial content, too high or too low RNA features, red blood cells (RBCs), and doublets called from genotype analysis
# singlecell.filter <-
#  singlecell %>%
#  mutate(
#    seurat2 = map2( name, seurat,
#                    function(l, s) {
#                      if (grepl("smallb", l)) {
#                        return(s)
#                      }
#                      subset(
#                        s,
#                       subset = percent.mito < 0.15 & 
#                          nFeature_RNA > 1500 & 
#                          nFeature_RNA < 5000 & 
#                          percent.rbc < 0.0003 &
#                          donor_id.stringent != "doublet"
#                      )
#                    }
#    )
#  )

# Save table with before and after filtering data
# singlecell.filter %>%
#  mutate(
#    before.count = map( seurat,  ~length(colnames(.x)) ),
#    after.count  = map( seurat2, ~length(colnames(.x)) ),
#  ) %>%
#  dplyr::select( name, before.count, after.count ) %>%
#  unnest( cols = c(before.count, after.count)) %>%
#  mutate(
#    fraction = after.count / before.count
#  ) %>%
#  # saveRDS( file=paste0( intermediate.prefix, "02_filter_table.rds"))

# Replace original seurat objects with filtered seurat objects, check file size & save
# singlecell.filter$seurat <- singlecell.filter$seurat2
# singlecell.filter$seurat2 <- NULL
# format(object.size(singlecell.filter), units = "Gb", standard = "auto", digits = 1L )
# saveRDS( singlecell.filter, file=paste0( intermediate.prefix, "02_filter.rds" ) )
# singlecell.filter <- readRDS( file = "20231126_filter_processed_annot.rds" ) ######## Read in saved filtered data


## SCtransform and cluster (basic)
defaultW <- getOption("warn")
options(warn=-1)

# Cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

singlecell.filter <-
  singlecell.filter %>%
  mutate(
    seurat = map( seurat,
                  function(s) {
                    DefaultAssay( s ) <- "RNA"
                    s <- NormalizeData( object = s, verbose = T )
                    s <- FindVariableFeatures( object = s, verbose = T )
                    # s <- ScaleData( object = s, vars.to.regress = c("nCount_RNA", "percent.mito") )
                    s <- CellCycleScoring(s, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)                    
                    # Now run SCT
                    s <- SCTransform(
                      object = s,
                      return.only.var.genes = FALSE,
                      method="glmGamPoi",
                      # vars.to.regress = c("S.Score", "G2M.Score"),
                      # vars.to.regress = "CC.Difference",
                      min_cells = 1 # to make sure all genes are returned.
                    )
                    s <- RunPCA( object = s, npcs = 50 )
                    s <- RunUMAP( object = s, dims = 1:35, verbose = FALSE )
                    s <- FindNeighbors( object = s, dims = 1:35, verbose = FALSE )
                    s <- FindClusters(
                      object = s,
                      algorithm = 4,
                      resolution = 0.4,  # seems less number of clusters in the end
                      verbose = FALSE )
                    return(s)
                  }
    )
  )

# Identities for e14 midgut and hindgut, Figure S7
mid.idents <- c("Epithelium", "SmoothMuscle1", "Subepithelial1", "Subepithelial2", "Inner Circ",
                "Mature Circ", "ICC1", "Mesenchyme", "ICC2", "Endothelium",
                "Neural2", "PDGFRAhi", "Neural3", "Cycling", "SmoothMuscle2", "ICC3", "Mesothelium")
hind.idents <- c("Inner Circ+Mature Circ", "Subepithelial", "Mesenchyme", "GREM1+BMPER+", "Epithelium",
                 "Cycling1", "Endothelium", "ICC1", "Cycling2", "ICC2",
                 "Neural1", "PDGFRAhi", "SmoothMuscle1", "Neural2", "G2M Circ",
                 "Epithelium", "SmoothMuscle2", "Neural3", "Immune", "Mesothelium")
names(mid.idents) <- levels(singlecell.filter$seurat[[1]])
singlecell.filter$seurat[[1]] <- RenameIdents(singlecell.filter$seurat[[1]], mid.idents)
names(hind.idents) <- levels(singlecell.filter$seurat[[2]])
singlecell.filter$seurat[[2]] <- RenameIdents(singlecell.filter$seurat[[2]], hind.idents)

# Identities for e15 pPIDNB midgut, pPIDNB-HOXD13 midgut, hindgut, Figures 4 & S5
mid.idents <- c("NeuralMesenchyme", "SmoothMuscle3", "Mesenchyme2", "SmoothMuscle1", "Neural1",
                "SmoothMuscle2", "Mesenchyme3", "Epithelium1", "Mesenchyme4", "Neural2",
                "Mesenchyme5", "Mesenchyme6", "Endothelium", "ICC", "Neural3", "SmoothMuscle4", "Neural4",
                "Immune")
hind.idents <- c("NeuralMesenchyme", "Mesenchyme4", "Epithelium1", "SmoothMuscle1", "Neural1",
                 "SmoothMuscle2", "ICC", "Mesenchyme7", "Endothelium", "Neural2",
                 "Neural3", "Mesenchyme5", "SmoothMuscle3", "Mesenchyme6", "Mesenchyme8",
                 "Neural4", "SmoothMuscle4", "Immune", "Mesenchyme9")
hox.idents <- c("NeuralMesenchyme", "Mesenchyme4", "SmoothMuscle1", "Epithelium1", "SmoothMuscle3",
                "Mesenchyme3", "SmoothMuscle2", "Mesenchyme5", "Neural3", "Endothelium",
                "Neural2", "ICC", "SmoothMuscle4", "Mesenchyme6", "Mesenchyme8",
                "Epithelium2", "Immune")

## Plotting
options( repr.plot.width = 7, repr.plot.height = 7 )

# UMAP for each sample
map(
  singlecell.filter$seurat,
  ~DimPlot( 
    .x, 
    reduction = "umap",
    label = T
  ) + theme(text = element_text(size = 10)
  ) + plot_annotation(title = .x@project.name)
)

# Hoxd13 expression endogenous Fig 4A
FeaturePlot( 
  singlecell.filter$seurat[[1]], 
  features = "HOXD13",
  order = T
)

# Just isolate and re-cluster mesenchyme clusters Figs 4B, E, S5D, S7F
hind.mes <- singlecell.filter$seurat[[2]] 
DefaultAssay(hind.mes) <- "RNA"
hind.mes <- subset(hind.mes, idents = c("Mesenchyme4", "SmoothMuscle2", "Mesenchyme5","ICC")) # replace with metadata idents of interest
# Figure S5
hind.mes <- NormalizeData( object = hind.mes, verbose = T )
hind.mes <- FindVariableFeatures( object = hind.mes, verbose = T )
hind.mes <- SCTransform(
  object = hind.mes,
  return.only.var.genes = FALSE,
  method="glmGamPoi",
  # vars.to.regress = c("S.Score", "G2M.Score"),
  # vars.to.regress = "CC.Difference",
  min_cells = 1 # to make sure all genes are returned.
)
hind.mes <- RunPCA( object = hind.mes, npcs = 50 )
hind.mes <- RunUMAP( object = hind.mes, dims = 1:30, verbose = FALSE )
hind.mes <- FindNeighbors( object = hind.mes, dims = 1:30, verbose = FALSE )
hind.mes <- FindClusters(
  object = hind.mes,
  algorithm = 4,
  resolution = 0.1,  # seems less number of clusters in the end
  verbose = FALSE )

# Plot expression of a gene just for subset of clusters
DimPlot( hind.mes, reduction = "umap",
         label = T)
FeaturePlot( 
  hind.mes, 
  features = "ACTA2",
  order = T
)

# Dotplot Fig 4B
pal <- viridis(n = 10, option = "mako")
DotPlot_scCustom(
  subset(hind, idents = c(
    "Mesenchyme4", "Mesenchyme7", "SmoothMuscle2", "ICC")
  ),
  scale.by = "size",
  features = rev(c("HOXD13", 
                   "LINGO2", "FAM132A", "ANO1",
                   "ACTG2", "TAGLN", "MYH11",
                   "KIT", "ADAMTSL1", "LMO3",
                   "COL1A1", "COL14A1", "COL1A2",
                   "SPON2", "SHISA2", "SERPIND1",
                   "CCDC102B", "SPRY3", "GUCY1A3",
                   "MASTL", "GTSE1", "PRC1", "CENPE", "NDC80", "NUSAP1")),
  colors_use = wes_palette("Zissou1", 100, type = "continuous"),
  remove_axis_titles = TRUE,
  x_lab_rotate = TRUE,
  flip_axes = TRUE
)


# Averaged expression for a gene set Fig S7F
gene.set <- c("INHBA","TGFB2")

# Get mean expression of genes of interest per cell
mean.exp <- colMeans(x = hind.mes@assays[["SCT"]]@data[gene.set, ], na.rm = TRUE)

# Add mean expression values in 'object@meta.data$gene.set.score'
if (all(names(x = mean.exp) == rownames(x = hind.mes@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  hind.mes@meta.data$gene.set.score <- mean.exp
}

# Plot mean expression using Seurat::FeaturePlot()
FeaturePlot_scCustom(seurat_object = hind.mes, features = "gene.set.score",
                     colors_use = rev(pal), na_cutoff = NULL)


## Differential gene expression
# Cluster markers
cluster_markers <- map(
  singlecell.filter$seurat,
  ~FindAllMarkers( 
    .x, 
    only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
  ) 
)

# # Save .csv for each set of cluster markers
# write.csv(cluster_markers[[1]],"mid_cluster_markers.csv")
# write.csv(cluster_markers[[2]],"hind_cluster_markers.csv")
# write.csv(cluster_markers[[3]],"hox_cluster_markers.csv")

# # Look at top markers per cluster
# top10.markers.hox <-
#   cluster_markers[[3]] %>%
#     group_by(cluster) %>%
#     slice_max(n = 10, order_by = avg_log2FC)


## Pseudobulk analysis

# Add classification for transgene expression 
# Annotate by mNeon, Hox transgene, and endogenous Hox expression NOTE: already present in publicly available Seurat objects!
cell.class <- map(
  singlecell.filter$seurat, 
  ~ FetchData(.x, 
              vars = c("cHoxd13co", "mNeonGreen.P2A.TetOn", "HOXD13", "orig.ident")
  )
)
ctrl <- map(
  cell.class, 
  ~ gsub("TRUE", "Ctrl", gsub( "FALSE", "",
                               if(any(.x[["rna_HOXD13"]]==0)) {
                                 .x[["mNeonGreen.P2A.TetOn"]]>0 & .x[["cHoxd13co"]]==0
                               } 
                               else {
                                 .x[["HOXD13"]]==0
                               }
  ))
)
hox <- map(
  cell.class, 
  ~ gsub("TRUE", "Hox+", gsub("FALSE", "",
                              if(any(.x[["rna_HOXD13"]]==0 & .x[["orig.ident"]]=="HGL8")) {
                                .x[["cHoxd13co"]]>0
                              } 
                              else if(any(.x[["rna_HOXD13"]]==0 & .x[["orig.ident"]]=="HGL6")) {
                                .x[["orig.ident"]]!="HGL6"
                              } 
                              else {
                                .x[["HOXD13"]]>0
                              }
  ))
)


# Merge cluster-based cell annotations with Hox expression-based annotations
names(mid.idents) <- levels(singlecell.filter$seurat[[1]])
singlecell.filter$seurat[[1]] <- RenameIdents(singlecell.filter$seurat[[1]], mid.idents)
singlecell.filter$seurat[[1]]$hox.exp <- paste0(hox[[1]], ctrl[[1]])
singlecell.filter$seurat[[1]]$cell.type <- paste0(singlecell.filter$seurat[[1]]$hox.exp, 
                                                  as.character(Idents(
                                                    singlecell.filter$seurat[[1]])))
singlecell.filter$seurat[[1]]$gut.cond <- rep("Mid", 4410)

names(hind.idents) <- levels(singlecell.filter$seurat[[2]])
singlecell.filter$seurat[[2]] <- RenameIdents(singlecell.filter$seurat[[2]], hind.idents)
singlecell.filter$seurat[[2]]$hox.exp <- paste0(hox[[2]], ctrl[[2]])
singlecell.filter$seurat[[2]]$cell.type <- paste0(singlecell.filter$seurat[[2]]$hox.exp, 
                                                  as.character(Idents(
                                                    singlecell.filter$seurat[[2]])))
singlecell.filter$seurat[[2]]$gut.cond <- rep("Hind", 9022)

names(hox.idents) <- levels(singlecell.filter$seurat[[3]])
singlecell.filter$seurat[[3]] <- RenameIdents(singlecell.filter$seurat[[3]], hox.idents)
singlecell.filter$seurat[[3]]$hox.exp <- paste0(hox[[3]], ctrl[[3]])
singlecell.filter$seurat[[3]]$cell.type <- paste0(singlecell.filter$seurat[[3]]$hox.exp,
                                                  as.character(Idents(
                                                    singlecell.filter$seurat[[3]])))
singlecell.filter$seurat[[3]]$gut.cond <- rep("Mid+Hox", 5217)


# Add metadata annotations for sex NOTE: already present in publicly available Seurat objects!
singlecell.filter$seurat[[1]]$sex <- singlecell.filter$seurat[[1]]$donor_id.stringent
singlecell.filter$seurat[[1]]$sex <- gsub("donor0", "M", singlecell.filter$seurat[[1]]$sex)
singlecell.filter$seurat[[1]]$sex <- gsub("donor1", "F", singlecell.filter$seurat[[1]]$sex)

singlecell.filter$seurat[[2]]$sex <- singlecell.filter$seurat[[2]]$donor_id.stringent
singlecell.filter$seurat[[2]]$sex <- gsub("donor0", "F", singlecell.filter$seurat[[2]]$sex)
singlecell.filter$seurat[[2]]$sex <- gsub("donor1", "M", singlecell.filter$seurat[[2]]$sex)

singlecell.filter$seurat[[3]]$sex <- singlecell.filter$seurat[[3]]$donor_id.stringent
singlecell.filter$seurat[[3]]$sex <- gsub("donor0", "M", singlecell.filter$seurat[[3]]$sex)
singlecell.filter$seurat[[3]]$sex <- gsub("donor1", "M", singlecell.filter$seurat[[3]]$sex)
singlecell.filter$seurat[[3]]$sex <- gsub("donor2", "F", singlecell.filter$seurat[[3]]$sex)


# Set up metadata as desired for aggregation and DE analysis
# Create single cell experiment object
counts <- map(
  singlecell.filter$seurat, 
  ~.x[["RNA"]]$counts
)
metadata <- map(
  singlecell.filter$seurat,
  ~.x[[]]
)

mid.sce <- SingleCellExperiment(
  assays = list(counts = counts[[1]]),
  colData = metadata[[1]])
hind.sce <- SingleCellExperiment(
  assays = list(counts = counts[[2]]),
  colData = metadata[[2]])
hox.sce <- SingleCellExperiment(
  assays = list(counts = counts[[3]]),
  colData = metadata[[3]])
sce <- cbind(mid.sce, hind.sce, hox.sce)

# Lump identities for pseudobulk DE analysis
# mid.sce$cell.type <- gsub("Mesenchyme2", "Mesenchyme4", mid.sce$cell.type)
# mid.sce$cell.type <- gsub("CtrlMesenchyme2", "CtrlMesenchyme4", mid.sce$cell.type)
# hind.sce$cell.type <- gsub("[<>+]", "", hind.sce$cell.type)
# hind.sce$cell.type <- gsub("HoxMesenchyme7", "HoxMesenchyme4", hind.sce$cell.type)
# hind.sce$cell.type <- gsub("CtrlMesenchyme7", "CtrlMesenchyme4", hind.sce$cell.type)

rm(mid.sce, hind.sce, hox.sce)

# Identify groups for aggregation of counts
colData(sce)

# Aggregate counts according to specified groups
reduced_sce <- pseudobulk(
  sce, group_by = vars(gut.cond, donor_id.stringent, cell.type, Phase, sex), 
  n_cells = n())
colData(reduced_sce)

# Filter entries that are not useful
reduced_sce.filter <- reduced_sce[, 
                                  reduced_sce$n_cells > 10 &
                                    reduced_sce$donor_id.stringent!="unassigned"] 
colData(reduced_sce.mid)

## Code for further subsetting
reduced_sce.mid <- reduced_sce.filter[,
                                      reduced_sce.filter$gut.cond=="Mid" &
                                        reduced_sce.filter$cell.type=="ICC"]
reduced_sce.hind <- reduced_sce.filter[,
                                       reduced_sce.filter$gut.cond=="Hind" &
                                         reduced_sce.filter$cell.type=="HoxICC"]
reduced_sce.hox <- reduced_sce.filter[,
                                      reduced_sce.filter$gut.cond=="Mid+Hox" &
                                        # reduced_sce.filter$donor_id.stringent != "donor0" &
                                        reduced_sce.filter$cell.type=="Hox+ICC"]
reduced_sce.filter <- cbind(reduced_sce.mid, reduced_sce.hox)
rm(reduced_sce.mid, reduced_sce.hind, reduced_sce.hox)
colData(reduced_sce.filter)

## Run Gamma Poisson model pseudobulk analysis
fit <- glm_gp(
  reduced_sce.filter, 
  design = ~ cell.type, # ~ cell.type + gut.cond
  size_factors = "deconvolution", verbose = TRUE)
colnames(fit$Beta)

# Find differentially expressed genes for a given comparison
res.midhind <- test_de(
  fit, 
  contrast = cond(
    cell.type = "Hox+ICC") - cond(
      cell.type = "ICC")
)

res.midhind.up <- res.midhind[rev(order(res.midhind$lfc)),]
res.midhox.up <- res.midhox[rev(order(res.midhox$lfc)),]
res.midhind.down <- res.midhind[order(res.midhind$lfc),]
res.midhox.down <- res.midhox[order(res.midhox$lfc),]


## Plotting pseudobulk results
# Volcano plot
# Make column in results object with genes to label
gene.labels.midhind.up = res.midhind.up$name[res.midhind.up$adj_pval < 0.01 &
                                               res.midhind.up$lfc>1]
gene.labels.midhox.up = res.midhox.up$name[res.midhox.up$pval < 0.05 &
                                             res.midhox.up$lfc>1]
heatmap.genes.up <- intersect(gene.labels.midhind.up, gene.labels.midhox.up)
exclude = grep("LOC", heatmap.genes.up, value = T)
heatmap.genes.up <- heatmap.genes.up[!heatmap.genes.up %in% exclude]

gene.labels.midhind.down = res.midhind.down$name[res.midhind.down$adj_pval < 0.01 &
                                                   res.midhind.down$lfc< (-1)]
gene.labels.midhox.down = res.midhox.down$name[res.midhox.down$pval < 0.05 &
                                                 res.midhox.down$lfc< (-1)]
heatmap.genes.down <- intersect(gene.labels.midhind.down, gene.labels.midhox.down)
exclude = grep("LOC", heatmap.genes.down, value = T)
heatmap.genes.down <- heatmap.genes.down[!heatmap.genes.down %in% exclude]

# Plot results
library(ggplot2, warn.conflicts = FALSE)
p <- ggplot(res, aes(x = lfc, y = -log2(adj_pval), label = name)) +
  geom_point(aes(color = pval < 0.01), size = 1)
p+geom_text_repel(aes(label=label), size=3, nudge_y = 0.5, max.overlaps = 6) #+ 
# ggtitle("pPIDNB-cHoxd13co Cluster 5 vs. pPIDNB Cluster 2") +
# theme( 
#   plot.title = element_text(size=14, face="bold"),
# )

EnhancedVolcano(
  res.midhind,
  res.midhind$name,
  x="lfc",
  y="adj_pval",
  # selectLab = c('FGL2','GREM1','LVRN','HOXC6','HOXC8'),
  # selectLab = c('DCLK1','NELL2','FGL2','CHRM3', 'FAM20C','TENM4',
  #               'LAMA3','CHRM2','EPHA6','KIF26A'),
  xlim = c(-10,10),
  # ylim = c(0, 6),
  xlab = bquote(~log[2] ~ "fold change"),
  ylab = bquote(~-log[10] ~ italic(p) ~ "-adj."),
  axisLabSize = 18,
  title = NULL,
  subtitle = NULL,
  # title = "Volcano plot",
  # subtitle = bquote(italic(EnhancedVolcano)),
  # caption = paste0("total = ", nrow(toptable), " variables"),
  # titleLabSize = 18,
  # subtitleLabSize = 14,
  captionLabSize = 14,
  pCutoff = 1e-02,
  FCcutoff = 1,
  cutoffLineType = "longdash",
  cutoffLineCol = "gray",
  cutoffLineWidth = 0.4,
  pointSize = 2,
  labSize = 0.5,
  labCol = "black",
  labFace = "plain",
  # labhjust = 0.5,
  # labvjust = 1.5,
  boxedLabels = FALSE,
  shape = 19,
  shapeCustom = NULL,
  col = c("cornsilk3", "cadetblue3", "royalblue", "coral"),
  colCustom = NULL,
  colAlpha = 1/2,
  colGradient = NULL,
  colGradientBreaks = c(pCutoff, 1),
  colGradientLabels = c("0", "1.0"),
  colGradientLimits = c(0, 1),
  # legendLabels = c("NS", expression(log[2] ~ FC), "adj. p-value", expression(p - value ~ and
  #                                                                       ~ log[2] ~ FC)),
  legendPosition = "none",
  legendLabels = NULL,
  # legendPosition = "top",
  # legendLabSize = 14,
  # legendIconSize = 5,
  # legendDropLevels = TRUE,
  encircle = NULL,
  encircleCol = "black",
  encircleFill = "pink",
  encircleAlpha = 3/4,
  encircleSize = 2.5,
  shade = NULL,
  shadeFill = "grey",
  shadeAlpha = 1/2,
  shadeSize = 0.01,
  shadeBins = 2,
  drawConnectors = FALSE,
  widthConnectors = 0.5,
  typeConnectors = "closed",
  endsConnectors = "first",
  lengthConnectors = unit(0.01, "npc"),
  colConnectors = "grey10",
  arrowheads = TRUE,
  hline = NULL,
  hlineType = "longdash",
  hlineCol = "black",
  hlineWidth = 0.4,
  vline = NULL,
  vlineType = "longdash",
  vlineCol = "black",
  vlineWidth = 0.4,
  gridlines.major = FALSE,
  gridlines.minor = TRUE,
  border = "partial",
  borderWidth = 0.8,
  borderColour = "black",
  raster = FALSE
)

# Heatmaps
# subset cells of interest
cells <- c(
  # WhichCells(subset(seurat_integrated, subset = cell.type == "Hox+SmoothMuscle1" &
  #                     orig.ident == "HGL7")),
  # # WhichCells(subset(seurat_integrated, subset = cell.type == "Hox+SmoothMuscle4" &
  # #          orig.ident == "HGL7")),
  WhichCells(subset(seurat_integrated, subset = cell.type == "CtrlICC" &
                      orig.ident == "HGL6")),
  # WhichCells(subset(seurat_integrated, subset = cell.type == "CtrlMesenchyme4" &
  #          orig.ident == "HGL6")),
  WhichCells(subset(seurat_integrated, subset = cell.type == "Hox+ICC" &
                      orig.ident == "HGL8"))),
#          #   ,
#   WhichCells(subset(seurat_integrated, subset = cell.type == "Hox+SmoothMuscle4" &
# orig.ident == "HGL8")))
# ,
# WhichCells(subset(seurat_integrated, subset = cell.type == "Hox+SmoothMuscle3" &
#                     orig.ident == "HGL7")),
WhichCells(subset(seurat_integrated, subset = cell.type == "CtrlSmoothMuscle3" &
                    orig.ident == "HGL6")),
WhichCells(subset(seurat_integrated, subset = cell.type == "Hox+SmoothMuscle3" &
                    orig.ident == "HGL8")))
# WhichCells(subset(seurat_integrated, subset = cell.type == "Hox+Mesenchyme4" &
#                     orig.ident == "HGL7")),
# WhichCells(subset(seurat_integrated, subset = cell.type == "Hox+Mesenchyme8" &
#                     orig.ident == "HGL7")),
# WhichCells(subset(seurat_integrated, subset = cell.type == "Mesenchyme2" &
#                     orig.ident == "HGL6")),
# WhichCells(subset(seurat_integrated, subset = cell.type == "Mesenchyme4" &
#                     orig.ident == "HGL6")),
# WhichCells(subset(seurat_integrated, subset = cell.type == "Hox+Mesenchyme4" &
#                     orig.ident == "HGL8")),
# WhichCells(subset(seurat_integrated, subset = cell.type == "Hox+Mesenchyme8" &
#                     orig.ident == "HGL8"))
)

seurat_integrated$heatmap <- paste0(seurat_integrated$cell.type, seurat_integrated$orig.ident)
seurat_integrated$heatmap <- gsub("[<>+]", "", seurat_integrated$heatmap)
seurat_integrated$heatmap <- gsub("CtrlSmoothMuscle1HGL6", "F", seurat_integrated$heatmap)
# seurat_integrated$heatmap <- gsub("SmoothMuscle4HGL6", "D", seurat_integrated$heatmap)
seurat_integrated$heatmap <- gsub("HoxSmoothMuscle1HGL8", "E", seurat_integrated$heatmap)
# seurat_integrated$heatmap <- gsub("HoxSmoothMuscle4HGL8", "E", seurat_integrated$heatmap)
seurat_integrated$heatmap <- gsub("HoxSmoothMuscle1HGL7", "D", seurat_integrated$heatmap)
# seurat_integrated$heatmap <- gsub("HoxSmoothMuscle4HGL7", "F", seurat_integrated$heatmap)
seurat_integrated$heatmap <- gsub("CtrlMesenchyme4HGL6", "C", seurat_integrated$heatmap)
seurat_integrated$heatmap <- gsub("CtrlICCHGL6", "B", seurat_integrated$heatmap)
seurat_integrated$heatmap <- gsub("HoxICCHGL8", "A", seurat_integrated$heatmap)
seurat_integrated$heatmap <- gsub("Mesenchyme2HGL6", "J", seurat_integrated$heatmap)
seurat_integrated$heatmap <- gsub("Mesenchyme4HGL6", "K", seurat_integrated$heatmap)
seurat_integrated$heatmap <- gsub("HoxMesenchyme4HGL8", "L", seurat_integrated$heatmap)
seurat_integrated$heatmap <- gsub("HoxMesenchyme8HGL8", "M", seurat_integrated$heatmap)
seurat_integrated$heatmap <- gsub("HoxMesenchyme4HGL7", "N", seurat_integrated$heatmap)
seurat_integrated$heatmap <- gsub("HoxMesenchyme8HGL7", "O", seurat_integrated$heatmap)


# heatmap.genes <- c(heatmap.genes.up.SM1, heatmap.genes.down.SM1) 
# heatmap.genes <- c(heatmap.genes.up.SM1, heatmap.genes.down.SM1,
#                    heatmap.genes.up.SM3, heatmap.genes.down.SM3,
#                    heatmap.genes.up.M4, heatmap.genes.down.M4)
heatmap.genes <- c(heatmap.genes.up, heatmap.genes.down)
saveRDS(heatmap.genes, file="heatmap_genes.rds" )

heatmap.genes.sm1.reordered <- c("INHBA", "STARD9", "PDE5A","CHRM3", "AQP3", "GRID2", "RTN1")
heatmap.genes.sm3.reordered <- c("NELL2","DCLK1","FGL2","CHRM3",  "LAMA3", "EPHA6") 
heatmap.genes.sm1.reordered.short <- c("INHBA", "CHRM3", "PDE5A", 
                                       "RTN1", "GFRA1", "GRID2")
heatmap.genes.sm3.reordered.short <- c("FGL2","ENPP2","DCLK1","LAMA3", "CHRM2","EPHA6") 

library("wesanderson")
pal <- wes_palette("Chevalier1")
view(pal)
pal <- c("#FF0000", "#FD6467", "#C6CDF7", "#46ACC8")

DoHeatmap(
  seurat_integrated,
  group.by = "heatmap",
  assay = "SCT",
  # slot = "counts",
  features = heatmap.genes,
  cells = cells,
  size = 3
) +
  scale_fill_gradientn(colours = rev(pal))


# Dataset integration
split_seurat <- singlecell.filter$seurat
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
seurat_integrated <- RunPCA(object = seurat_integrated)
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP for integrated data                             
pal <- wes_palette("Moonrise3")
pal <- c(pal, "#FD6467")
DimPlot(seurat_integrated,
        reduction = "umap",
        group.by = "final.calls.rescued",
        # group.by = "gut.cond",
        label=T)

# Violin plot with integrated data Fig 4D
seurat_integrated$violin <- paste0(seurat_integrated$cell.type, seurat_integrated$orig.ident)
seurat_integrated$violin <- gsub("[<>+]", "", seurat_integrated$violin)
seurat_integrated$violin <- gsub("HoxSmoothMuscle1HGL7", "SM", seurat_integrated$violin)
seurat_integrated$violin <- gsub("HoxSmoothMuscle1HGL8", "SM", seurat_integrated$violin)
seurat_integrated$violin <- gsub("CtrlSmoothMuscle1HGL6", "SM", seurat_integrated$violin)


pal <- wes_palette("Moonrise3")
subset <- subset(seurat_integrated, subset = violin == "SM")
Stacked_VlnPlot(seurat_object = subset, 
                features = c("INHBA", "STARD9","PDE5A", "CHRM3", "AQP3", "RTN1"),
                assay = "RNA",
                split.by = "gut.cond",
                plot_legend = TRUE,
                pt.size = 0.3,
                colors_use = pal)


# Violin plot cluster markers Fig S7B
# color palettes
moonrise <- wes_palette("Moonrise3")
grandbuda <- wes_palette("GrandBudapest1")
darj <- wes_palette("Darjeeling1")
pal.mid <- c(moonrise[1:2], grandbuda[2], darj[3])
pal.hind <- c(moonrise[1:3], darj[3], grandbuda[2] )

# Group clusters according to broad identity
hind_vln <- subset(hind, idents=c("Subepithelial","PDGFRAhi","GREM1+BMPER+","Mesenchyme","Cycling1","Cycling2"))
mid_vln <- subset(mid, idents=c("Subepithelial1","Subepithelial2", "PDGFRAhi","Mesenchyme","Cycling"))
mid.idents <- c("Subepithelial","Subepithelial","Mesenchyme","PDGFRAhi","Cycling")
hind.idents <- c("Subepithelial","Mesenchyme","GREM1+BMPER+","Cycling","Cycling","PDGFRAhi")
names(mid.idents) <- levels(mid_vln)
mid_vln <- RenameIdents(mid_vln, mid.idents)
names(hind.idents) <- levels(hind_vln)
hind_vln <- RenameIdents(hind_vln, hind.idents)

# Plot results
Stacked_VlnPlot(seurat_object = hind_vln,
                features = c("GREM1","BMPER","RSPO3"),
                assay = "RNA",
                # split.by = "cell.type",
                plot_legend = FALSE,
                pt.size = 0.2,
                colors_use=pal.hind,
                y.max=4,
                same.y.lims=TRUE)
Stacked_VlnPlot(seurat_object = mid_vln, 
                features = c("GREM1","BMPER","RSPO3"),
                assay = "RNA",
                # split.by = "cell.type",
                plot_legend = FALSE,
                pt.size = 0.2,
                colors_use=pal.mid,
                y.max=5,
                same.y.lims=TRUE)
