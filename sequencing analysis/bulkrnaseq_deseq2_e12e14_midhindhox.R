### Code for analysis of bulk RNA-seq data 
## Written by John C. Lawlor and Hasreet K. Gill

## Setup
# Note: Any step that has three asterisks "***" right after the hashmark requires the paths or filenames to be modified to fit how the data is stored/named.
# *** set working directory
setwd("Dropbox/My PC (MCB-HASREETGILL-1)/Documents/Tabin Lab/Molecular biology/rnaseq/rnaseq output")

# Bioconductor/CRAN libraries
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
library(AnnotationHub)
library(ensembldb)
library(apeglm)
library(tibble)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Gg.eg.db)
library(org.Hs.eg.db)
library(enrichplot)
library(ggnewscale)
library(ashr)
library(readxl)
library("genefilter")
library(ggplot2)
library(grid)
library(gtable)
library(viridis)
library(gridExtra)
library(ggdendro)
library(reshape2)
library("gplots")



###


## Create the annotation
# Access chick data from Annotation Hub, subset w/up to date annotations
ah <- AnnotationHub()
chick_ens <- query(ah, c("Gallus gallus", "EnsDb"))
chick_ens <- chick_ens[["AH89419"]]

# Create a gene-level dataframe 
annotations_ahb <- genes(chick_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, gene_name, entrezid, gene_biotype)

# Determine the indices for the non-duplicated genes, keep only non-duplicates
non_duplicates_idx <- which(duplicated(annotations_ahb$gene_name) == FALSE)
annotations_ahb <- annotations_ahb[non_duplicates_idx, ]

# Keep only first identifier in multiple mapping cases
annotations_ahb$entrezid <- map(annotations_ahb$entrezid,1) %>%  unlist()

# Create a transcript dataframe
txdb <- transcripts(chick_ens, return.type = "data.frame") %>%
  dplyr::select(tx_id, gene_id)
txdb <- txdb[grep("ENSGALT", txdb$tx_id),]

# Create a gene-level dataframe
genedb <- genes(chick_ens, return.type = "data.frame")  %>%
  dplyr::select(gene_id, gene_name)

# Merge the two dataframes together to create the tx2gene 
annotations <- inner_join(txdb, genedb)


###


## Making DESEQ2 dataset from salmon output files
# *** List all directories containing data, create vector of filenames with path. Name each quant file uniquely.
samples <- list.files(path = "C:/Users/hag070/e12_salmon_output/", full.names = T, pattern="salmon$")
# samples <- list.files(path = "C:/Users/hag070/e14_salmon_output/", full.names = T, pattern="salmon$")
files <- file.path(samples, "quant.sf")
names(files) <- str_replace(samples, "./e12_salmon_output/", "") %>%
# names(files) <- str_replace(samples, "./e14_salmon_output/", "") %>%
  str_replace(".salmon", "")

# Run tximport, write counts to an object. Note that ignoreTxVersion = TRUE in case Ensembl has decimaled transcript versions.
txi <- tximport(files, type="salmon", tx2gene=annotations[, c("tx_id", "gene_id")], countsFromAbundance="lengthScaledTPM", ignoreTxVersion = TRUE)
data <- txi$counts %>%
  round() %>%
  data.frame()

# *** Create a sampletable/metadata (necessary for DESEQ2)
sampletype <- factor(c(rep("GFP_mid_meso_E12",4), rep("GFP_hind_meso_E12", 4), rep("HOXD13_mid_meso_E12", 4)))
# sampletype <- factor(c(rep("GFP_mid_meso_E14",3), rep("GFP_hind_meso_E14", 3), rep("HOXD13_mid_meso_E14", 3)))

# samplenum <- (c(1:12))
# batchnum <- (c(1,2,3,4,1,2,3,3,1,1,2,3))
# meta <- data.frame(samplenum, batchnum, sampletype, row.names = colnames(txi$counts))
sampnames <- (c("mgfp_mg_1","mgfp_mg_2","mgfp_mg_3","mgfp_mg_4","mgfp_hg_1","mgfp_hg_2","mgfp_hg_3","mgfp_hg_4","hoxd13_mg_1","hoxd13_mg_2","hoxd13_mg_3","hoxd13_mg_4"))
# sampnames <- (c("mgfp_mg_1","mgfp_mg_2","mgfp_mg_3","mgfp_hg_1","mgfp_hg_2","mgfp_hg_3","hoxd13_mg_1","hoxd13_mg_2","hoxd13_mg_3"))
meta <- data.frame(sampletype, row.names = sampnames)
# meta <- data.frame(sampletype, row.names = colnames(txi$counts))


# Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)

## Replicate removal due to possible sample swap (informed by PCA analysis)
# 4, 6, 12, can be changed out to fit whichever replicates are desired to be removed
# first, sweep the environment then run everything preceding the creation of txi
# the rest of e12 and e14 Rscripts can be run as normal, now with the desired replicates removed

# files <- files[-c(4, 6, 12)]
# txi <- tximport(files, type="salmon", tx2gene=annotations[, c("tx_id", "gene_id")], countsFromAbundance="lengthScaledTPM", ignoreTxVersion = TRUE)
# sampletype <- factor(c(rep("GFP_mid_meso_E14",4), rep("GFP_hind_meso_E14", 4), rep("HOXD13_mid_meso_E14", 4)))
# sampletype <- sampletype[-c(4, 6, 12)]
# meta <- data.frame(sampletype, row.names = colnames(txi$counts))
# dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)

# PCA and further QC could occur here as well



###


## Differential expression analysis
# Run DESEQ
dds <- DESeq(dds)

# *** Define contrasts for RCAS-hoxd13 expression midgut (MODIFY THIS BASED ON DESIRED CONTRAST)
contrastGHo <- c("sampletype", "HOXD13_mid_meso_E12", "GFP_mid_meso_E12")
contrastGG <- c("sampletype", "GFP_hind_meso_E12", "GFP_mid_meso_E12")
# contrastGHo <- c("sampletype", "HOXD13_mid_meso_E14", "GFP_mid_meso_E14")
# contrastGG <- c("sampletype", "GFP_hind_meso_E14", "GFP_mid_meso_E14")

# Extract results for RCAS-hoxd13 expression midgut vs control midgut
res_table_GHo <- results(dds, contrast=contrastGHo, alpha = 0.05) # mid vs hox
res_table_GG <- results(dds, contrast=contrastGG, alpha = 0.05) # mid vs hind

#"Shrink" the log2 fold changes. First save the unshrunken results to compare
res_table_GHo_unshrunken <- res_table_GHo
res_table_GG_unshrunken <- res_table_GG
res_table_GHo <- lfcShrink(dds, contrast=contrastGHo, type="ashr")
res_table_GG <- lfcShrink(dds, contrast=contrastGG, type="ashr")

# Set significance thresholds
padj.cutoff <- 0.05
# log2fold.min <- 0.50

# Create a tibble of results
res_table_tb_GHo <- res_table_GHo %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_table_tb_GG <- res_table_GG %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset the tibble to keep only significant genes. Sig shows the log2 fold expression differences for the desired contrast!
sig_GHo <- res_table_tb_GHo %>%
  dplyr::filter(padj < padj.cutoff) 
  # dplyr::filter(abs(log2FoldChange) > log2fold.min)

sig_GG <- res_table_tb_GG %>%
  dplyr::filter(padj < padj.cutoff)
  # dplyr::filter(abs(log2FoldChange) > log2fold.min)

# Note that running this annotation thing multiple times screws up the columns in sig
sig_annot <- annotations %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::distinct()
sig_GHo <- merge(sig_GHo, sig_annot, by.x="gene", by.y="gene_id") %>%
  arrange(padj)
sig_GG <- merge(sig_GG, sig_annot, by.x="gene", by.y="gene_id") %>%
  arrange(padj)

sig_GGHo <- inner_join(sig_GHo, sig_GG, by="gene") %>%
  dplyr::filter(log2FoldChange.x*log2FoldChange.y > 0) %>%
  arrange(log2FoldChange.x)

sig_GGHo_e12 <- sig_GGHo
sig_GHo_e12 <- sig_GHo
sig_GG_e12 <- sig_GG

sig_GGHo_e14 <- sig_GGHo
sig_GHo_e14 <- sig_GHo
sig_GG_e14 <- sig_GG

sig_GGHo_e12_e14 <- inner_join(sig_GGHo_e12, sig_GGHo_e14, by="gene")
sig_GGHo_e14_unique <- sig_GGHo_e14 %>%
  dplyr::filter(!gene %in% sig_GGHo_e12_e14$gene)
sig_GGHo_e12_unique <- sig_GGHo_e12 %>%
  dplyr::filter(!gene %in% sig_GGHo_e12_e14$gene)


###

## Plotting
# Enhanced volcano plot 
sig_GGHo_e12_e14 <- read_excel("e12e14_venn_sig_midhindhoxd13.xlsx")
tgfbeta_ecm_genes <- read_excel("Tabin Lab/Molecular biology/rnaseq/tgfbeta_ecm_genes.xlsx") # tgfbeta ecm gene list extracted directly from Enrichr gene set "TGF-beta regulation of extracellular matrix"
tgfbeta_ecm_genes <- inner_join(tgfbeta_ecm_genes, sig_GGHo_e12_e14, by="gene_name.x.x") %>%
  dplyr::filter(log2FoldChange.x.x*log2FoldChange.y.x > 0) %>%
  arrange(-log2FoldChange.x.x)

tgf_annot <- annotations %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::distinct()
tgfbeta_ecm_genes <- merge(tgfbeta_ecm_genes, tgf_annot, by.x="gene_name.x.x", by.y="gene_name")%>%
  arrange(-log2FoldChange.x.x)
# tgfbeta_ecm_genes <- merge(tgfbeta_ecm_genes, tgf_annot, by.x="gene_name.x", by.y="gene_name")%>%
#   arrange(-log2FoldChange.x)


if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
res_table_tb_GG_e12 <- res_table_tb_GG
res_table_tb_GHo_e12 <- res_table_tb_GHo
# res_table_tb_GG_e14 <- res_table_tb_GG
# res_table_tb_GHo_e14 <- res_table_tb_GHo

# res_table_tb_volc <- res_table_tb_GG_e12 %>%
res_table_tb_volc <- res_table_tb_GG_e14 %>%
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.50)
galgal6annot <- annotations %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::distinct()
res_table_tb_volc <- bind_cols(res_table_tb_volc, gene_name=galgal6annot$gene_name[match(res_table_tb_volc$gene, galgal6annot$gene_id)])
res = subset(res_table_tb_volc, select = c(gene_name,log2FoldChange,padj)) 

hindhox <- setdiff(sig_GGHo_e14$gene_name.x,sig_GGHo_e12_e14$gene_name.x.x)
e12e14 <- sig_GGHo_e12_e14$gene_name.x.x
tgfe12e14 <- setdiff(tgfbeta_ecm_genes$gene_name.x.x,sig_GGHo_e12_e14$gene_name.x.x)
tgfe14unique <- tgfbeta_ecm_genes$gene_name.x

keyvals.colour <- ifelse(
  res$padj > 0.05, 'ivory2',
  ifelse(res$gene_name %in% e12e14, 'black',
         ifelse(res$gene_name %in% tgfe12e14, 'red1', 'azure4')))
# keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'ivory2'] <- 'hindhox'
names(keyvals.colour)[keyvals.colour == 'red1'] <- 'sigjusthind'
names(keyvals.colour)[keyvals.colour == 'azure4'] <- 'notsig'


 
  EnhancedVolcano(res,
                      lab = res$gene_name,
                      x = 'log2FoldChange',
                      y = 'padj',
                      selectLab = res$gene_name[which(names(keyvals.colour) %in% c('hindhox'))],
                      xlab = bquote(~Log[2]~ 'fold change'),
                      title = 'Custom shape & colour over-ride',
                      pCutoff = 10e-14,
                      FCcutoff = 1.0,
                      pointSize = 1,
                      labSize = 2,
                      # shapeCustom = keyvals.shape,
                      colCustom = keyvals.colour,
                      colAlpha = 1,
                      legendPosition = 'right',
                      legendLabSize = 15,
                      legendIconSize = 5.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      colConnectors = 'grey50',
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      border = 'full',
                      borderWidth = 1.0,
                      borderColour = 'black')



# Plot expression of only one gene
gene1_meta <- meta %>%
  rownames_to_column(var="samplename") %>%
  as_tibble()

# Convert normalized_counts to a data frame and transfer the row names to a new column called "gene"
normalized_counts <- counts(dds, normalized=T) %>%
  data.frame() %>%
  rownames_to_column(var="gene")

# Next, merge together (ensembl IDs) the normalized counts data frame with a subset of the annotations in the tx2gene data frame (only the columns for ensembl gene IDs and gene symbols)
galgal6annot <- annotations %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::distinct()

## This will bring in a column of gene symbols
normalized_counts <- merge(normalized_counts, galgal6annot, by.x="gene", by.y="gene_id")

# Now create a tibble for the normalized counts
normalized_counts <- normalized_counts %>%
  as_tibble()

# Find the Ensembl ID
galgal6annot[galgal6annot$gene_name == "GDF3", "gene_id"]

# Save plotcounts to a data frame object
d <- plotCounts(dds, gene="ENSGALG00000003161", intgroup="sampletype", returnData=TRUE)
write.csv(d, "e12_GDF3.csv")

# Plot the normalized counts, using the samplenames (rownames(d) as labels)
ggplot(d, aes(x = sampletype, y = count, color = sampletype)) +
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) +
  theme_bw() +
  ggtitle("HOXD13 E12") +
  theme(plot.title = element_text(hjust = 0.5))

#or (seems to work better):
plotCounts(dds, gene="ENSGALG00000034616", intgroup="sampletype")



# Gene expression heatmap
# Transform count data using the variance stablilizing transform
deseq2VST <- vst(dds)
# deseq2ResDF <- as.data.frame(res_table)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)

# Keep only the significantly differentiated genes where the fold-change was at least 3
# sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 0.2,])
# sigGenes <- clusters_heatmap$genes
deseq2VST <- deseq2VST[deseq2VST$Gene %in% tgfbeta_ecm_genes$gene_id,]

deseq2VST$Gene <- NULL
deseq2VSTMatrix <- data.matrix(deseq2VST, rownames.force = NA)
deseq2VSTMatrix <- deseq2VSTMatrix[order(factor(rownames(deseq2VSTMatrix), levels = tgfbeta_ecm_genes$gene)),]
par(mar = c(2,2,2,2))
tgf_annot <- annotations %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::distinct()
genenames <- data.frame(rownames(deseq2VSTMatrix))
tgfgenenames <- merge(genenames, tgf_annot, by.x="rownames.deseq2VSTMatrix.", by.y="gene_id", sort=FALSE)
heatmap <- heatmap.2(deseq2VSTMatrix,col=brewer.pal(11,"PiYG"), scale="row",
          key=TRUE, symkey=FALSE, lhei = c(2,7), density.info="none", trace="none", cexRow=.9, dendrogram = c("both","row","column","none"), Colv=FALSE, Rowv=FALSE, labRow=tgfgenenames$gene_name, labCol=FALSE, ColSideColors=rep(c("green","orange","purple"), each=4))

deseq2VST_dend <- melt(deseq2VST, id.vars=c("Gene"))
deseq2VSTMatrix_dend <- dcast(deseq2VST_dend, Gene ~ variable)
rownames(deseq2VSTMatrix_dend) <- deseq2VSTMatrix_dend$Gene
deseq2VSTMatrix_dend$Gene <- NULL
distanceGene <- dist(deseq2VSTMatrix_dend)
distanceSample <- dist(t(deseq2VSTMatrix_dend))

# Cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="average")
clusterSample <- hclust(distanceSample, method="average")

# Construct a dendogram for samples
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + theme_dendro()

geneModel <- as.dendrogram(clusterGene)
geneDendrogramData <- segment(dendro_data(geneModel, type = "rectangle"))

# construct the dendrogram in ggplot
geneDendrogram <- ggplot(geneDendrogramData) + geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + coord_flip() + scale_y_reverse(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) + theme_dendro()
