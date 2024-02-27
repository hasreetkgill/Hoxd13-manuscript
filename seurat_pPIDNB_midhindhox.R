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
# intermediate.prefix <- "20231126_"

# setwd("C:/Users/lalig/OneDrive/Documents/Tabin Lab/Hoxd13 project/scRNA-seq/hgl5")

singlecell.filter <- readRDS( file = "midhind_filter_processed_annotated.rds" ) ######## Read in saved filtered data
mid.markers <- readRDS( file = "midmarkers_122023.rds" ) ######## Read in saved filtered data
hind.markers <- readRDS( file = "hindmarkers_122023.rds" ) ######## Read in saved filtered data
seurat_integrated <- readRDS( file = "20231204_midhind_integrated.rds" ) ######## Read in saved filtered data
hind.mus <- readRDS(file = "hindmesclustersonly.rds")

# Load pre-processed Seurat objects & combine into one singlecell object for filtering/processing
hgl6 <- readRDS(file = "HGL6/20231030_HGL_reannotated_HGL6.rds") # pPIDNB midgut control
hgl7 <- readRDS(file = "HGL7/20231030_HGL_reannotated_HGL7.rds") # pPINDB/pPIDNB-cHoxd13co hindgut
hgl8 <- readRDS(file = "HGL8/20231030_HGL_reannotated_HGL8.rds") # pPIDNB-cHoxd13co midgut

singlecell <- tribble(
  ~seurat, ~name,
  hgl6, "HGL6",
  hgl7, "HGL7",
  hgl8, "HGL8",
)
mesmus.2 <- tribble(
  ~seurat, ~name,
  mid.mus, "midmesmus",
  hind.mus, "hindmesmus",
)

rm(hgl6, hgl7, hgl8)


##### Filter cells
# Remove cells with high mitochondrial content, too high or too low RNA features, red blood cells (RBCs), and doublets called from genotype analysis
singlecell.filter <-
  singlecell %>%
  mutate(
    seurat2 = map2( name, seurat,
                    function(l, s) {
                      if (grepl("smallb", l)) {
                        return(s)
                      }
                      subset(
                        s,
                        subset = percent.mito < 0.15 & 
                          nFeature_RNA > 1500 & 
                          nFeature_RNA < 5000 & 
                          percent.rbc < 0.0003 &
                          donor_id.stringent != "doublet"
                      )
                    }
    )
  )

# Save table with before and after filtering data
singlecell.filter %>%
  mutate(
    before.count = map( seurat,  ~length(colnames(.x)) ),
    after.count  = map( seurat2, ~length(colnames(.x)) ),
  ) %>%
  dplyr::select( name, before.count, after.count ) %>%
  unnest( cols = c(before.count, after.count)) %>%
  mutate(
    fraction = after.count / before.count
  ) %>%
  # saveRDS( file=paste0( intermediate.prefix, "02_filter_table.rds"))

# Replace original seurat objects with filtered seurat objects, check file size & save
singlecell.filter$seurat <- singlecell.filter$seurat2
singlecell.filter$seurat2 <- NULL
format(object.size(singlecell.filter), units = "Gb", standard = "auto", digits = 1L )
# saveRDS( singlecell.filter, file=paste0( intermediate.prefix, "02_filter.rds" ) )
singlecell.filter <- readRDS( file = "20231126_filter_processed_annot.rds" ) ######## Read in saved filtered data





##### SCtransform and cluster (basic)
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

hind.mes <- singlecell.filter$seurat[[3]]
DefaultAssay(hind.mes) <- "RNA"
hind.mes <- subset(hind.mes, idents = c("Mesenchyme4", "SmoothMuscle2", "Mesenchyme5","ICC"))
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
hind.mes.markers <- FindAllMarkers(hind.mes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
hox.mes.icc<-hind.mes
hox.mes.icc.markers <- hind.mes.markers

DimPlot( mid.mes.icc, reduction = "umap",
         label = T)
DimPlot( hind.mes.icc, reduction = "umap",
         label = T)
DimPlot( hox.mes.icc, reduction = "umap",
         label = T)
FeaturePlot( 
  hox.mes., 
  features = "cHoxd13co",
  order = T
)

hind.mus <- singlecell.filter$seurat[[3]]
DefaultAssay(hind.mus) <- "RNA"
hind.mus <- subset(hind.mus, idents = c("SmoothMuscle1", "SmoothMuscle3", "SmoothMuscle4"))
hind.mus <- NormalizeData( object = hind.mus, verbose = T )
hind.mus <- FindVariableFeatures( object = hind.mus, verbose = T )
hind.mus <- SCTransform(
  object = hind.mus,
  return.only.var.genes = FALSE,
  method="glmGamPoi",
  # vars.to.regress = c("S.Score", "G2M.Score"),
  # vars.to.regress = "CC.Difference",
  min_cells = 1 # to make sure all genes are returned.
)

# hind.mus <- subset(hind.mus, idents = c("SmoothMuscle1", "SmoothMuscle3", "SmoothMuscle4"))
hind.mus <- RunPCA( object = hind.mus, npcs = 50 )
hind.mus <- RunUMAP( object = hind.mus, dims = 1:15, verbose = FALSE )
hind.mus <- FindNeighbors( object = hind.mus, dims = 1:15, verbose = FALSE )
hind.mus <- FindClusters(
  object = hind.mus,
  algorithm = 4,
  resolution = 0.1,  # seems less number of clusters in the end
  verbose = FALSE )


ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
color_list <- ggplotColours(n=16)
DimPlot( hox.mes.icc, 
         reduction = "umap",
         label = T,
         cols = color_list[c(1,3,7,9,7,11)], 
         pt.size = 1)
# 5,8,10
# 1,4,12,2
mid <- singlecell.filter$seurat[[1]]
DimPlot( hox, 
         reduction = "umap",
         label = T)


pal <- viridis(n = 10, option = "plasma")
DefaultAssay(hox.mes.icc)<-"SCT"
FeaturePlot_scCustom(seurat_object = hox, features = "COL1A1", colors_use = rev(pal), pt.size=2)
FeaturePlot( 
  hox.mes.icc, 
  features = "ACTA2"
  # order = T
)


hind.mus.markers <- FindAllMarkers(hind.mus, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
hind.mus.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
hind.mes.markers <- FindAllMarkers(hind.mes, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
hind.mes.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

mid<- singlecell.filter$seurat[[1]]
hind.markers <- FindAllMarkers(hind, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
hind.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


levels(hind) <- c("SmoothMuscle3", "SmoothMuscle1", "SmoothMuscle4", "ICC",  "Mesenchyme4", "Mesenchyme7", "Mesenchyme5", "Mesenchyme8", "SmoothMuscle2","Mesenchyme6",
                  "Mesenchyme9", "NeuralMesenchyme", "Neural1", "Neural2", "Neural3", "Neural4", "Immune", "Epithelium1", "Endothelium")
levels(hind) <- rev(levels(hind))


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

DotPlot_scCustom(
  hind.mus,
  scale.by = "size",
  features = rev(c("ZNF536", 
                   "DACH1", "NRXN1", "KCNQ3",
                   "DLG2", "NPAS3", "COL14A1",
                   "SLIT2", "EBF1", "PDE3A",
                   "GREM1", "BMPER", "NEXN",
                   "PLAC9", "RSPO3", "TUBA1A1",
                   "HMGB2", "SMC2", "TOP2A",
                   "HMGB1", "DIAPH3", "MXD3", 
                   "CIT", "NDC80", "MKI67",
                   "KCNH1", "NDST4", "NRG1", 
                   "SLC35F3", "TRPA1")),
  colors_use = wes_palette("Zissou1", 100, type = "continuous"),
  remove_axis_titles = TRUE,
  x_lab_rotate = TRUE,
  flip_axes = TRUE
)

DotPlot_scCustom(
  hind.mus,
  scale.by = "size",
  features = rev(c("FBLN1", 
                   "THBS2", "VCAN", "TAGLN",
                   "MYL9", "NPAS3", "COL14A1",
                   "SLIT2", "EBF1", "PDE3A",
                   "GREM1", "BMPER", "NEXN",
                   "PLAC9", "RSPO3", "TUBA1A1",
                   "HMGB2", "SMC2", "TOP2A",
                   "HMGB1", "DIAPH3", "MXD3", 
                   "CIT", "NDC80", "MKI67",
                   "KCNH1", "NDST4", "NRG1", 
                   "SLC35F3", "TRPA1")),
  colors_use = wes_palette("Zissou1", 100, type = "continuous"),
  remove_axis_titles = TRUE,
  x_lab_rotate = TRUE,
  flip_axes = TRUE
)

DotPlot(object = subset(hind, idents = c("SmoothMuscle3", "SmoothMuscle1", "SmoothMuscle4", "ICC",  "Mesenchyme4", "Mesenchyme7", "Mesenchyme5", "Mesenchyme8")), 
        scale.by = "size",
        features = c("HOXD13", 
                     "LINGO2", "FAM132A", "ANO1",
                     "ACTG2", "TAGLN", "MYH11",
                     "KIT", "ADAMTSL1", "LMO3",
                     "COL1A1", "NPAS3", "COL1A2",
                     "SPON2", "SHISA2", "SERPIND1",
                     "CCDC102B", "SPRY3", "GUCY1A3",
                     "MASTL", "GTSE1", "PRC1", "CENPE", "NDC80", "NUSAP1"))

options(warn=defaultW)

# Save filtered & processed data
saveRDS( singlecell.filter, file = "midhind_filter_processed.rds")
singlecell.filter <- readRDS(file = "20231106_02_filter_processed.rds") ######## Read in saved filtered & processed data





##### Plotting
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

gene.set.muspooled <- names.mes.mus.sorted.muspooled
gene.set.mespooled <- names.mes.mus.sorted.mespooled
gene.set.mes <- names.mes.sorted

write.csv(as.data.frame(gene.set.muspooled), "geneset_mesvsmuspooled.csv")
write.csv(as.data.frame(gene.set.mespooled), "geneset_mespooledvsmus.csv")
write.csv(as.data.frame(gene.set.mes), "geneset_mesvscircmus.csv")

saveRDS( gene.set.muspooled, file = "geneset_mesvsmuspooled.rds")
saveRDS( gene.set.mespooled, file = "geneset_mespooledvsmus.rds")
saveRDS( gene.set.mes, file = "geneset_mesvscircmus.rds")


saveRDS( seurat_integrated.mus, file = "integrated_mus.rds")


# gene.set <- c("CDH1","TCF12","KRT7","KLF5",
#               "EEF2")
# gene.set <- c("TAGLN","FBLN1","THBS2","MYL9",
#               "VCAN","FHL3","TUBB","LRIG1")
gene.set <- c("INHBA","TGFB2")
  
# Get mean expression of genes of interest per cell
mean.exp <- colMeans(x = hind.mus@assays[["SCT"]]@data[gene.set, ], na.rm = TRUE)

# Add mean expression values in 'object@meta.data$gene.set.score'
if (all(names(x = mean.exp) == rownames(x = hind.mus@meta.data))) {
  cat("Cell names order match in 'mean.exp' and 'object@meta.data':\n", 
      "adding gene set mean expression values in 'object@meta.data$gene.set.score'")
  hind.mus@meta.data$gene.set.score <- mean.exp
}

# Plot mean expression using Seurat::FeaturePlot()
FeaturePlot_scCustom(seurat_object = hind.mus, features = "gene.set.score",
                     colors_use = rev(pal), na_cutoff = NULL)
FeaturePlot(object = hind, features = "ACTA2")


####################################################
mid.idents <- c("Epithelium1", "SmoothMuscle1", "Neural/Mesenchyme", "Neural1", "Inner Circ",
                "Mature Circ", "ICC1", "Mesenchyme5", "ICC2", "Endothelium",
                "Neural2", "Mesenchyme6", "Neural3", "G2/M+Immune", "SmoothMuscle2", "ICC3", "Mesothelium")
hind.idents <- c("Inner Circ+Mature Circ", "Neural/Mesenchyme", "Mesenchyme5", "Mesenchyme7", "Epithelium1",
                 "G2/M1", "Endothelium", "ICC1", "G2/M2", "ICC2",
                 "Neural2", "Mesenchyme6", "SmoothMuscle2", "Neural3", "G2M Circ",
                 "Epithelium2", "SmoothMuscle1", "Neural2", "Immune", "Mesothelium")
names(mid.idents) <- levels(singlecell.filter$seurat[[1]])
singlecell.filter$seurat[[1]] <- RenameIdents(singlecell.filter$seurat[[1]], mid.idents)
names(hind.idents) <- levels(singlecell.filter$seurat[[2]])
singlecell.filter$seurat[[2]] <- RenameIdents(singlecell.filter$seurat[[2]], hind.idents)



hind.mus <- singlecell.filter$seurat[[2]]
# DefaultAssay(hind.mus) <- "RNA"
# hind.mus <- subset(hind.mus, idents = c("Neural/Mesenchyme", "Mesenchyme5", "Mesenchyme7","Mesenchyme6", "SmoothMuscle2","SmoothMuscle1","Inner Circ+Mature Circ"))
hind.mus <- subset(hind.mus, idents = c("Neural/Mesenchyme", "Mesenchyme5", "Mesenchyme7","Mesenchyme6", "SmoothMuscle2"))
hind.mus <- NormalizeData( object = hind.mus, verbose = T )
hind.mus <- FindVariableFeatures( object = hind.mus, verbose = T )
hind.mus <- SCTransform(
  object = hind.mus,
  return.only.var.genes = FALSE,
  method="glmGamPoi",
  # vars.to.regress = c("S.Score", "G2M.Score"),
  # vars.to.regress = "CC.Difference",
  min_cells = 1 # to make sure all genes are returned.
)

# hind.mus <- subset(hind.mus, idents = c("SmoothMuscle1", "SmoothMuscle3", "SmoothMuscle4"))
hind.mus <- RunPCA( object = hind.mus, npcs = 50 )
hind.mus <- RunUMAP( object = hind.mus, dims = 1:15, verbose = FALSE )
hind.mus <- FindNeighbors( object = hind.mus, dims = 1:15, verbose = FALSE )
hind.mus <- FindClusters(
  object = hind.mus,
  algorithm = 4,
  resolution = 0.25,  # seems less number of clusters in the end
  verbose = FALSE )

DimPlot( hind,
         reduction = "umap",
         label = T)


mid.mus <- singlecell.filter$seurat[[1]]
# DefaultAssay(mid.mus) <- "RNA"
# mid.mus <- subset(mid.mus, idents = c("Neural/Mesenchyme", "Neural1", "Mesenchyme5","Mesenchyme6", "SmoothMuscle2","SmoothMuscle1","Inner Circ","Mature Circ"))
mid.mus <- subset(mid.mus, idents = c("Neural/Mesenchyme", "Neural1", "Mesenchyme5","Mesenchyme6", "SmoothMuscle2"))
mid.mus <- NormalizeData( object = mid.mus, verbose = T )
mid.mus <- FindVariableFeatures( object = mid.mus, verbose = T )
mid.mus <- SCTransform(
  object = mid.mus,
  return.only.var.genes = FALSE,
  method="glmGamPoi",
  # vars.to.regress = c("S.Score", "G2M.Score"),
  # vars.to.regress = "CC.Difference",
  min_cells = 1 # to make sure all genes are returned.
)

# hind.mus <- subset(hind.mus, idents = c("SmoothMuscle1", "SmoothMuscle3", "SmoothMuscle4"))
mid.mus <- RunPCA( object =mid.mus, npcs = 50 )
mid.mus <- RunUMAP( object = mid.mus, dims = 1:15, verbose = FALSE )
mid.mus <- FindNeighbors( object = mid.mus, dims = 1:15, verbose = FALSE )
mid.mus <- FindClusters(
  object = mid.mus,
  algorithm = 4,
  resolution = 0.25,  # seems less number of clusters in the end
  verbose = FALSE )

mus <- tribble(
  ~seurat, ~name,
  mid.mus, "Mid",
  hind.mus, "Hind",
)

FeaturePlot( 
  hind, 
  features = c("KLF5", "RPLP0", "KRT7", "RPL5", "RPSAP58",
               "RPL4", "EEF2", "NTN3", "HLX", "NKX3-2", "HSPA2",
               "ETV6"),
  order = T)

# Feature plot for same gene across different samples
map(
  singlecell.filter$seurat,
  ~FeaturePlot( 
    .x, 
    features = c( grep( "ACTG2", rownames(.x), value = T )),
    order = T
  ) + plot_annotation(title = .x@project.name)
)

# Feature plot for genes in one sample
mid <- singlecell.filter$seurat[[1]]
hind <- singlecell.filter$seurat[[2]]
hox <- singlecell.filter$seurat[[3]]
hox.donor0 <-hox[, 
                 hox$donor_id.stringent=="donor0"]
hox.donor1 <- hox[, 
                 hox$donor_id.stringent=="donor1"]
hox.donor2 <- hox[, 
                  hox$donor_id.stringent=="donor2"]
mid.donor0 <-mid[, 
                 mid$donor_id.stringent=="donor0"]
mid.donor1 <- mid[, 
                  mid$donor_id.stringent=="donor1"]
FeaturePlot(object = mid.donor1, features = "HINTW", combine = FALSE)
lapply(X = fp, FUN = function(x) x + theme(text = element_text(size = 10)))
CombinePlots(plots = fp)




##### Differential gene expression
# Cluster markers
cluster_markers <- map(
  singlecell.filter$seurat,
  ~FindAllMarkers( 
    .x, 
    only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
  ) 
)


hind.mus.markers <- FindAllMarkers(hind.mus, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# # Save .csv for each set of cluster markers
# write.csv(cluster_markers[[1]],"mid_cluster_markers.csv")
# write.csv(cluster_markers[[2]],"hind_cluster_markers.csv")
# write.csv(cluster_markers[[3]],"hox_cluster_markers.csv")

# # Look at top markers per cluster
# top10.markers.hox <-
#   cluster_markers[[3]] %>%
#     group_by(cluster) %>%
#     slice_max(n = 10, order_by = avg_log2FC)
# 


# Pseudobulk analysis
# Add classification for transgene expression
# Annotate by mNeon, Hox transgene, and endogenous Hox expression
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

# Cell types reasoned from cluster markers
# mid.idents <- c("NeuralMesenchyme", "Mesenchyme1", "Mesenchyme2", "SmoothMuscle1", "Neural",
#                            "SmoothMuscle2", "Mesenchyme3", "Epithelium1", "Mesenchyme4", "DividingNeural",
#                            "Mesenchyme5", "Mesenchyme6", "Endothelium", "ICC", "Neural2", "DividingMesenchyme1",
#                            "Neural3", "Immune")
# hind.idents <- c("NeuralMesenchyme", "Mesenchyme4", "Epithelium1", "SmoothMuscle1", "Neural",
#                            "SmoothMuscle2", "ICC", "Mesenchyme7", "Endothelium", "DividingNeural",
#                            "Neural2", "Mesenchyme5", "Mesenchyme1", "Mesenchyme6", "DividingMesenchyme1",
#                            "Neural3", "DividingMuscle", "Immune", "DividingMesenchyme2")
# hox.idents <- c("NeuralMesenchyme", "Mesenchyme4", "SmoothMuscle1", "Epithelium1", "Mesenchyme1",
#                             "Mesenchyme3", "SmoothMuscle2", "Mesenchyme5", "Neural2", "Endothelium",
#                             "DividingNeural", "ICC", "DividingMuscle", "Mesenchyme6", "DividingMesenchyme1",
#                             "Epithelium2", "Immune")

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

# ??????????????
# hox <- singlecell.filter$seurat[[3]]
# hox <- subset(hox, subset = donor_id.stringent != "donor2")
# sm1.hox.cell.names <- row.names(subset(hox, subset = cell.type == "Hox+SmoothMuscle1")@meta.data)
# levels(hox@active.ident) = c(levels(hox@active.ident), 18)
# levels(hox@meta.data$seurat_clusters) = c(levels(hox@meta.data$seurat_clusters), 18)
# hox@active.ident[which(row.names(hox@meta.data) %in%  sm1.hox.cell.names)] = 18
# hox@meta.data$seurat_clusters[which(row.names(hox@meta.data) %in%  sm1.hox.cell.names)] = 18
# 
# x <- FindMarkers(hox, ident.1 = 18)
# DotPlot(object = hind, features = "INHBA")
# 
# 
# 
# new.hox.idents <- paste0(singlecell.filter$seurat[[3]]$hox.exp,
#                      as.character(Idents(
#                        singlecell.filter$seurat[[3]])))
# new.hox.idents <- levels(singlecell.filter$seurat[[3]])
# singlecell.filter$seurat[[3]] <- RenameIdents(singlecell.filter$seurat[[3]], new.hox.idents)
# test <- FindMarkers(singlecell.filter$seurat[[3]], hox.exp = "Hox+")








rm(hox,ctrl, cell.class)

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



singlecell.filter <- readRDS(file = "20231126_filter_processed_annot.rds") ######## Read in saved filtered & processed data
mid <- singlecell.filter$seurat[[1]]
hind <- singlecell.filter$seurat[[2]]
hox <- singlecell.filter$seurat[[3]]

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

mid.sce$cell.type <- gsub("Mesenchyme2", "Mesenchyme4", mid.sce$cell.type)
mid.sce$cell.type <- gsub("CtrlMesenchyme2", "CtrlMesenchyme4", mid.sce$cell.type)


hind.sce <- SingleCellExperiment(
  assays = list(counts = counts[[2]]),
  colData = metadata[[2]])

# hind.sce$cell.type <- gsub("[<>+]", "", hind.sce$cell.type)
# 
# hind.sce$cell.type <- gsub("HoxMesenchyme7", "HoxMesenchyme4", hind.sce$cell.type)
# hind.sce$cell.type <- gsub("CtrlMesenchyme7", "CtrlMesenchyme4", hind.sce$cell.type)


hox.sce <- SingleCellExperiment(
  assays = list(counts = counts[[3]]),
  colData = metadata[[3]])
sce <- cbind(mid.sce, hind.sce, hox.sce)

# sce$cell.type <- gsub("Mesenchyme8", "Mesenchyme4", sce$cell.type)
# sce$cell.type <- gsub("Hox+Mesenchyme8", "Hox+Mesenchyme4", sce$cell.type)
# sce$cell.type <- gsub("CtrlMesenchyme8", "CtrlMesenchyme4", sce$cell.type)

# sce$cell.type <- gsub("SmoothMuscle4", "SmoothMuscle1", sce$cell.type)
# sce$cell.type <- gsub("Hox+SmoothMuscle4", "Hox+SmoothMuscle1", sce$cell.type)
# sce$cell.type <- gsub("CtrlSmoothMuscle4", "CtrlSmoothMuscle1", sce$cell.type)

rm(mid.sce, hind.sce, hox.sce)

# Identify groups for aggregation of counts
colData(sce)

# sce$cell.type <- gsub("Mesenchyme4", "Mesenchyme", sce$cell.type)
# sce$cell.type <- gsub("Mesenchyme8", "Mesenchyme", sce$cell.type)
# sce$cell.type <- gsub("Hox+Mesenchyme4", "Hox+Mesenchyme", sce$cell.type)
# sce$cell.type <- gsub("Hox+Mesenchyme8", "Hox+Mesenchyme", sce$cell.type)

# Aggregate counts according to specified groups
reduced_sce <- pseudobulk(
  sce, group_by = vars(gut.cond, donor_id.stringent, cell.type, Phase, sex), 
  n_cells = n())
colData(reduced_sce)

# Filter entries that are not useful
reduced_sce.filter <- reduced_sce[, 
                             reduced_sce$n_cells > 10 &
                             reduced_sce$donor_id.stringent!="unassigned"] 
# colData(reduced_sce.mid)
                             
############ Code for further subsetting

#midgut mesenchyme2=hindgut mesenchyme4?
#midgut mesenchyme4=hindgut mesenchyme7?
reduced_sce.mid <- reduced_sce.filter[,
                               reduced_sce.filter$gut.cond=="Mid" &
                                 reduced_sce.filter$cell.type=="ICC"]
# reduced_sce.mid2 <- reduced_sce.filter[,
#                                        reduced_sce.filter$gut.cond=="Mid" &
#                                          reduced_sce.filter$cell.type=="Mesenchyme"]
reduced_sce.hind <- reduced_sce.filter[,
                               reduced_sce.filter$gut.cond=="Hind" &
                               reduced_sce.filter$cell.type=="HoxICC"]
# # reduced_sce.hind2 <- reduced_sce.filter[,
#                                         reduced_sce.filter$gut.cond=="Hind" &
#                                           reduced_sce.filter$cell.type=="Hox+Mesenchyme"]
reduced_sce.hox <- reduced_sce.filter[,
                               reduced_sce.filter$gut.cond=="Mid+Hox" &
                                 # reduced_sce.filter$donor_id.stringent != "donor0" &
                                 reduced_sce.filter$cell.type=="Hox+ICC"]

reduced_sce.filter <- cbind(reduced_sce.mid, reduced_sce.hox)
rm(reduced_sce.mid, reduced_sce.hind, reduced_sce.hox)
# reduced_sce.filter@colData@listData[["donor_id.stringent"]] <- gsub("donor2","donor1",
                                                                    # reduced_sce.filter@colData@listData[["donor_id.stringent"]])
colData(reduced_sce.filter)

# colData(reduced_sce_filter)
############



# Run Gamma Poisson model pseudobulk analysis
# group <- factor(
#   paste0(
#     reduced_sce.filter$gut.cond, 
#     reduced_sce.filter$cell.type)
# )

fit <- glm_gp(
  reduced_sce.filter, 
  design = ~ cell.type, # ~ cell.type + gut.cond
  size_factors = "deconvolution", verbose = TRUE)
colnames(fit$Beta)

# saveRDS( fit, file=paste0( intermediate.prefix, "fit_MidMid+Hox_CtrlHox+_SmoothMuscle1.rds" ) )
# fit <- readRDS(file = "20231113_fit_celltype_SmoothMuscle1.rds") ######## Read in fit

# Find differentially expressed genes for a given comparison
res.midhox.ICC.v2 <- test_de(
  fit, 
  contrast = cond(
    cell.type = "Hox+ICC") - cond(
      cell.type = "ICC")
    )
    # cond(
    # cell.type = "Hox+SmoothMuscle1", gut.cond = "Mid+Hox") - cond(
    #   cell.type = "CtrlSmoothMuscle1", gut.cond = "Mid"))

saveRDS( res.midhox.ICC.v2, file="res_MidvHox_ICC_Phase.rds" )
# write.csv(res,"SM1_Hox+Ctrl_MidvMidHox_DEGs.csv")
# # res <- readRDS(file = "res_Hox+CtrlSM1_MidvMidHox.rds") ######## Read in DE results







res.midhind<-res.midhind.ICC
res.midhox<-res.midhox.ICC

# res.midhind <- readRDS(file = "res_MidvHind_SM3_Phase.rds") ######## Read in DE results
# res.midhox <- readRDS(file = "res_MidvHox_SM3_Phase_noDonor0.rds") ######## Read in DE results

res.midhind.up <- res.midhind[rev(order(res.midhind$lfc)),]
res.midhox.up <- res.midhox[rev(order(res.midhox$lfc)),]
res.midhind.down <- res.midhind[order(res.midhind$lfc),]
res.midhox.down <- res.midhox[order(res.midhox$lfc),]


##### Plotting pseudobulk results
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

# saveRDS( heatmap.genes.M4.2, file="heatmap_genes_MidMesClus2&4_HindHoxClus4.rds" )


# shared.SM3.M4 <- intersect(heatmap.genes.M4, heatmap.genes.SM3)
# shared.SM3.M4 <- intersect(heatmap.genes.SM3, heatmap.genes.M4)
# saveRDS( shared.SM1.SM3, file="heatmap_sharedgenes_SMClus1&4SMClus3.rds" )

# res$label = if_else(res$name %in% gene.labels, 
#                     res$name, "")
# Plot results
library(ggplot2, warn.conflicts = FALSE)
p <- ggplot(res, aes(x = lfc, y = -log2(adj_pval), label = name)) +
  geom_point(aes(color = pval < 0.01), size = 1)
p+geom_text_repel(aes(label=label), size=3, nudge_y = 0.5, max.overlaps = 6) #+ 
  # ggtitle("pPIDNB-cHoxd13co Cluster 5 vs. pPIDNB Cluster 2") +
 # theme( 
 #   plot.title = element_text(size=14, face="bold"),
 # )


mid <- singlecell.filter$seurat[[1]]
hind <- singlecell.filter$seurat[[2]]
hox <- singlecell.filter$seurat[[3]]
singlecell.filter.all <- merge(mid, y = c(hind,hox), add.cell.ids = c("M", "H", "X"))
singlecell.filter.all[["RNA"]] <- split(singlecell.filter.all[["RNA"]], f = singlecell.filter.all$gut.cond)

res <- read.csv("InnerMes_HindvsMid_byPadj.csv", header = TRUE)

EnhancedVolcano(
  res.midhind.ICC,
  res.midhind.ICC$name,
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


EnhancedVolcano(
  res.midhox,
  res.midhox$name,
  selectLab = c('DCLK1','NKX3-2','NELL2','FGL2','ENPP2','ISL2','CHRM3','FAM20C','TENM4',
                'LAMA3','CHRM2','EPHA6','KIF26A'),
  x="lfc",
  y="pval",
  # selectLab = NULL,
  xlim = c(-10,10),
  ylim = c(0, 8),
  xlab = bquote(~log[2] ~ "fold change"),
  ylab = bquote(~-log[10] ~ italic(p)),
  axisLabSize = 18,
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
  labSize = 2,
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
  legendLabels = c("NS", expression(log[2] ~ FC), "adj. p-value", expression(p - value ~ and
                                                                             ~ log[2] ~ FC)),
  legendPosition = "top",
  legendLabSize = 14,
  legendIconSize = 5,
  legendDropLevels = TRUE,
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
# saveRDS( cells, file="heatmap_cells.rds" )
saveRDS( cells, file="heatmap_cells_midnoMesCtrls.rds" )



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

SM3 <- DoHeatmap(
  seurat_integrated,
  group.by = "heatmap",
 assay = "SCT",
 # slot = "counts",
  features = heatmap.genes.sm1.reordered,
  cells = cells,
  size = 3
) +
  scale_fill_gradientn(colours = rev(pal))


saveRDS( singlecell.filter, file=paste0( "20231126_filter_processed_annot.rds" ) )
split_seurat <- mesmus.2$seurat
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
seurat_integrated.mesmus.2 <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

saveRDS( seurat_integrated, file="20240122_midhind_integrated_mesenchyme.rds" )
seurat_integrated <- readRDS( file="20231126_seurat_integrated.rds" )



seurat_integrated <- RunPCA(object = seurat_integrated.mesmus)

# Plot PCA
# PCAPlot(seurat_integrated,
#         split.by = "gut.cond")  


set.seed(123456)

# Run UMAP
seurat_integrated <- RunUMAP(seurat_integrated, 
                             dims = 1:40,
                             reduction = "pca")

# Plot UMAP                             
pal <- wes_palette("Moonrise3")
pal <- c(pal, "#FD6467")
DimPlot(seurat_integrated,
        reduction = "umap",
        group.by = "final.calls.rescued",
        # group.by = "gut.cond",
        label=T)
        
 
DimPlot(seurat_integrated.mus,
        reduction = "umap",
        # group.by = "final.calls.rescued",
        # group.by = "ident",
        label = T)
DimPlot(hind.mus,
        reduction = "umap",
        # group.by = "final.calls.rescued",
        # group.by = "ident",
        label = T) 
DefaultAssay(mid.mus) <- "SCT"
FeaturePlot_scCustom( 
  seurat_object =hind.mus,
  features = c("FBLN1"),
  label = T,
  # colors_use = rev(pal),
  pt.size=2)

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

# Cell cycle scoring
# s.genes <- cc.genes$s.genes
# g2m.genes <- cc.genes$g2m.genes
# seurat <- singlecell.filter$seurat[[3]]
# 
# seurat <- NormalizeData(seurat)
# seurat <- FindVariableFeatures(seurat, selection.method = "vst")
# seurat <- ScaleData(seurat, features = rownames(seurat))
# seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# head(seurat[[]])
# seurat <- RunPCA(seurat, features = c(s.genes, g2m.genes))
# DimPlot(seurat)
# 
# marrow <- ScaleData(marrow, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(marrow))
# marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)
# DimPlot(marrow)
# 
# ## ALT METHOD
# marrow$CC.Difference <- marrow$S.Score - marrow$G2M.Score
# marrow <- ScaleData(marrow, vars.to.regress = "CC.Difference", features = rownames(marrow))
# # cell cycle effects strongly mitigated in PCA
# marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)
# marrow <- RunPCA(marrow, features = c(s.genes, g2m.genes))
# DimPlot(marrow)

# #### if making Seurat object
# meta.stats <-
#   meta.dir2 %>%
#   mutate(
#     path       = map2_chr( dir, name, function(x, y) { paste0( root.dir, x, "/", y, "/outs/", "metrics_summary.csv" ) } ),
#     library    = name
#   ) %>%
#   mutate(
#     metrics = map( path,
#                    ~read_csv( . ) %>%as.data.frame() )
#   ) %>%
#   unnest( cols = "metrics")
# 
# # colnames( meta.stats )
# 
# meta.stats %>%
#   mutate_if( is.integer,
#              scales::comma_format()
#   )  %>%
#   select( -path )
# 
# singlecell <-
#   meta.dir2 %>%
#   mutate(
#     seurat = map2(dir, name, function( d, l ) {
#         warning(paste0(l, "\n"))  # log
#         m <- Read10X(glue::glue("{root.dir}{d}/{l}/outs/filtered_feature_bc_matrix/"))
#         s <- CreateSeuratObject(counts = m, project = l)
#         return(s)
#       }
#     )
# )

# singlecell <-
#   singlecell %>%
#   mutate(
#     seurat = pmap( list(seurat, dir, name),
#                    function( s, d, l  ) {
#                      s@project.name <- l
#                      s$library <- l
#                      s$barcode <- rownames( s@meta.data )
#                      s$barcode2 <- paste0( d, "-", l, "-", s$barcode )
#                      return(s)
#                    }
#     )
#   )
# 
# 
# mito.genes.chick <- c("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
# gene.list.chick <- rownames(singlecell$seurat[[1]])
# 
# read_tsv("FC_08283/uniprotkb_ribosomal_protein_AND_gallus_gallus_2023_10_18.tsv.gz") %>% head()
# 
# read_tsv("FC_08283/uniprotkb_ribosomal_protein_AND_gallus_gallus_2023_10_18.tsv.gz") %>% 
#   dplyr::count( `Reviewed` )
# 
# 
# # This makes the entire gene list
# provisional  <-
#   bind_rows(
#     read_tsv("FC_08283/uniprotkb_ribosomal_protein_AND_gallus_gallus_2023_10_18.tsv.gz", show_col_types = F),
#   ) %>%
#   dplyr::select(`Gene Names`, `Organism` ) %>%
#   dplyr::filter( !is.na(`Gene Names`) ) %>%
#   separate_rows( `Gene Names`, convert = FALSE, sep = " " ) %>%
#   pull(`Gene Names` )
# #dplyr::filter( grepl("14e22", `Gene Names` ) )
# 
# 
# provisional %>% length() # How many distinct names do you find in the UniProt list?
# 
# intersect( gene.list.chick, provisional ) %>% length() # How many gene names can you retrieve by simple intersection?
# # Potential genes that might be missed
# potential <- setdiff( provisional, gene.list.chick ) # Left over names that do not have a match
# 
# potential
# 
# grep("RPL", gene.list.chick, value = T) %>% sort()
# 
# ribo.genes.chick <- intersect( gene.list.chick, c(provisional, grep("RPL", gene.list.chick, value = T)) )
# ribo.genes.chick %>% length()
# 
# 
# grep("^HB", gene.list.chick, value = T)
# 
# rbc.genes.chick <- c("HBE", "HBBA", "HBBR", "HBZ", "HBM")
# 
# 
# singlecell <- 
#   singlecell %>%
#   mutate(
#     seurat = map2( seurat, name,
#                    function( s, l ) {
#                      warning( paste0(l, "\n") )
# 
#                        s$percent.mito <- Matrix::colSums(GetAssayData(object = s, slot = "counts")[mito.genes.chick, ])/Matrix::colSums(GetAssayData(object = s, slot = "counts"))
#                        s$percent.rbc <- Matrix::colSums(GetAssayData(object = s, slot = "counts")[rbc.genes.chick, ])/Matrix::colSums(GetAssayData(object = s, slot = "counts"))
#                        s$percent.ribo <- Matrix::colSums(GetAssayData(object = s, slot = "counts")[ribo.genes.chick, ])/Matrix::colSums(GetAssayData(object = s, slot = "counts"))
# 
#                      return(s)
#                    }
#     )
#   )
# 
# 
# 
# # The following is repeated for some reason?
# # singlecell <-
# #   singlecell %>%
# #   mutate(
# #     seurat = pmap( list(seurat, dir, name),
# #                    function( s, d, l  ) {
# #                      s@project.name <- l
# #                      s$library <- l
# #                      s$barcode <- rownames( s@meta.data )
# #                      s$barcode2 <- paste0( d, "-", l, "-", s$barcode )
# #                      return(s)
# #                    }
# #     )
# #   )
# 
# intermediate.prefix <- "20231102_HGL8"
saveRDS(hind.mus, file= "hindmesclustersonly.rds" )
# # data <- readRDS(file = "29231102_.rds")
# 
# data <- singlecell$seurat[[1]]


