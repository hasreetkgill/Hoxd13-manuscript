## Setup
#Note: Any step that has three asterisks "***" right after the hashmark requires the paths or filenames to be modified to fit how the data is stored/named.
#set working directory
setwd("Dropbox/My PC (MCB-HASREETGILL-1)/Documents/Tabin Lab/Molecular biology/rnaseq/rnaseq output")

#Bioconductor/CRAN libraries
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


## Making DESEQ2 dataset
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


# 
# # *** Create a sampletable/metadata (necessary for DESEQ2)
sampletype <- factor(c(rep("GFP_mid_meso_E12",4), rep("GFP_hind_meso_E12", 4), rep("HOXD13_mid_meso_E12", 4)))
# sampletype <- factor(c(rep("GFP_mid_meso_E14",3), rep("GFP_hind_meso_E14", 3), rep("HOXD13_mid_meso_E14", 3)))

#samplenum <- (c(1:12))
#batchnum <- (c(1,2,3,4,1,2,3,3,1,1,2,3))
# meta <- data.frame(samplenum, batchnum, sampletype, row.names = colnames(txi$counts))
sampnames <- (c("mgfp_mg_1","mgfp_mg_2","mgfp_mg_3","mgfp_mg_4","mgfp_hg_1","mgfp_hg_2","mgfp_hg_3","mgfp_hg_4","hoxd13_mg_1","hoxd13_mg_2","hoxd13_mg_3","hoxd13_mg_4"))
# sampnames <- (c("mgfp_mg_1","mgfp_mg_2","mgfp_mg_3","mgfp_hg_1","mgfp_hg_2","mgfp_hg_3","hoxd13_mg_1","hoxd13_mg_2","hoxd13_mg_3"))
meta <- data.frame(sampletype, row.names = sampnames)
# meta <- data.frame(sampletype, row.names = colnames(txi$counts))


# # 
# # Create DESeq2Dataset object
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)

## Replicate removal
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
#Run DESEQ
dds <- DESeq(dds)

# *** Define contrasts for RCAS-hoxd13 expression midgut (MODIFY THIS BASED ON DESIRED CONTRAST)
contrastGHo <- c("sampletype", "HOXD13_mid_meso_E12", "GFP_mid_meso_E12")
contrastGG <- c("sampletype", "GFP_hind_meso_E12", "GFP_mid_meso_E12")
# contrastGHo <- c("sampletype", "HOXD13_mid_meso_E14", "GFP_mid_meso_E14")
# contrastGG <- c("sampletype", "GFP_hind_meso_E14", "GFP_mid_meso_E14")

# Extract results for RCAS-hoxd13 expression midgut vs control midgut
#res_table <- results(dds, contrast=contrast_h, alpha = 0.05)
res_table_GHo <- results(dds, contrast=contrastGHo, alpha = 0.05)
res_table_GG <- results(dds, contrast=contrastGG, alpha = 0.05)

#"Shrink" the log2 fold changes. First save the unshrunken results to compare
res_table_GHo_unshrunken <- res_table_GHo
res_table_GG_unshrunken <- res_table_GG

#there is some problem with this. coef can only found by resultsNames(dds), but that won't allow comparison of the desired conditions (Hoxd13 mid versus GFP mid)
#because the contrast is different than the coef this has no effect right now. Will need to be corrected
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

# #Note that running this annotation thing multiple times screws up the columns in sig
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
# 
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
# ###

# count_data <- data
# count_data <- count_data
# colnames(count_data) <- sampnames
# count_data$gene <- rownames(count_data)

# count_data_e12 <- data
# count_data_e12 <- count_data
# colnames(count_data_e12) <- sampnames
# count_data_e12$gene <- rownames(count_data_e12)

# count_data_e14 <- data
# count_data_e14 <- count_data
# colnames(count_data_e14) <- sampnames
# count_data_e14$gene <- rownames(count_data_e14)

# count_data_GHo_e12 <- count_data %>%
#   dplyr::filter(gene %in% sig_GHo$gene) %>%
#   dplyr::select(-gene)
# count_data_GG_e12 <- count_data %>%
#   dplyr::filter(gene %in% sig_GG$gene) %>%
#   dplyr::select(-gene)
# count_data_GGHo_e12 <- count_data %>%
#   dplyr::filter(gene %in% sig_GGHo$gene) %>%
#   dplyr::select(-gene)
# count_data_e12_GGHo_e12_e14 <- count_data %>%
#   dplyr::filter(gene %in% sig_GGHo_e12_e14$gene) %>%
#   dplyr::select(-gene)

# count_data_GHo_e14 <- count_data %>%
#   dplyr::filter(gene %in% sig_GHo$gene) %>%
#   dplyr::select(-gene)
# count_data_GG_e14 <- count_data %>%
#   dplyr::filter(gene %in% sig_GG$gene) %>%
#   dplyr::select(-gene)
# count_data_GGHo_e14 <- count_data %>%
#   dplyr::filter(gene %in% sig_GGHo$gene) %>%
#   dplyr::select(-gene)
# count_data_e14_GGHo_e12_e14 <- count_data %>%
#   dplyr::filter(gene %in% sig_GGHo_e12_e14$gene) %>%
#   dplyr::select(-gene)

# write.csv(count_data, "e12_counts_all.csv")
# write.csv(count_data_GHo, "e12_counts_mid_hoxd13_sig.csv")
# write.csv(count_data_GG, "e12_counts_mid_hind_sig.csv")
# write.csv(count_data_GGHo, "e12_counts_venn_sig.csv")

# 
# write.csv(sig_GGHo_e12_unique, "e12_unique_sig_midhindhoxd13.csv")
# write.csv(count_data_GHo, "e14_counts_mid_hoxd13_sig.csv")
# write.csv(count_data_GG, "e14_counts_mid_hind_sig.csv")
# write.csv(count_data_GGHo, "e14_counts_venn_sig.csv")


# res_ids_GHo <- left_join(res_table_tb_GHo, annotations_ahb, by=c("gene"="gene_id"))
# res_ids_GG <- left_join(res_table_tb_GG, annotations_ahb, by=c("gene"="gene_id"))
# 
# # Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
# all_genes <- as.character(res_ids_GHo$gene)
# all_genes <- as.character(res_ids_GG$gene)
# 
# GGHo_GO <- res_ids_GG  %>%
#   
# # GGHo_GO <- res_ids_GHo  %>%
#   dplyr::filter(gene %in% sig_GGHo$gene)
# GGHo_GO_genes <- as.character(GGHo_GO$gene)
# ego_GGHo <- enrichGO(gene = GGHo_GO_genes, universe = all_genes, keyType = "ENSEMBL", OrgDb = org.Gg.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE, pool = TRUE)
# cluster_summary_GGHo <- data.frame(ego_GGHo)
# 
# dotplot(ego_GGHo)
# 
# 
# # ## Likelihood ratio test & gene clustering... note that this can find a different number of significant genes than the log2 fold analysis, so the gene vectors don't quite match.
# # LRT significant results tibble
# dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)
# res_LRT <- results(dds_lrt)
# res_LRT_tb <- res_LRT %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>%
#   as_tibble()
# sigLRT_genes <- res_LRT_tb %>%
#   filter(padj < padj.cutoff)
# # # # 
# # # 
# # # # Improve distances/clustering (also useful for PCA)
# rld <- rlog(dds, blind=TRUE)
# rld_mat <- assay(rld)
# # 
# # #PCA
# colnames(rld_mat) <- c('gfp_mg_1','gfp_mg_2','gfp_mg_3','gfp_hg_1','gfp_hg_2','gfp_hg_3','hox_mg_1','hox_mg_2','hox_mg_3')
# plotPCA(rld, intgroup="batchnum")
# rld_cor <- cor(rld_mat)
# heatmap(rld_cor)
# # # 
# # # # *** Subset results for faster cluster finding 
#  clustering_sig_genes <- sigLRT_genes %>%
#    arrange(padj)
#  cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
# # #
# # # # Gene clusters across samples
# clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col=NULL)
# # # 
# # # # Extract genes by group and annotate
# group1 <- clusters$df %>%
#   filter(cluster == 1)
# row.names(group1) <- 1:nrow(group1)
# group1_annotated <- genes(chick_ens, return.type = "data.frame")  %>%
#   dplyr::select(gene_id, gene_name, entrezid, gene_biotype) %>%
#   dplyr::filter(gene_id %in% group1$gene)
# #
# group2 <- clusters$df %>%
#   filter(cluster == 2)
# row.names(group2) <- 1:nrow(group2)
# group2_annotated <- genes(chick_ens, return.type = "data.frame")  %>%
#   dplyr::select(gene_id, gene_name, entrezid, gene_biotype) %>%
#   dplyr::filter(gene_id %in% group2$gene)
# 
# group3 <- clusters$df %>%
#   filter(cluster == 3)
# row.names(group3) <- 1:nrow(group3)
# group3_annotated <- genes(chick_ens, return.type = "data.frame")  %>%
#   dplyr::select(gene_id, gene_name, entrezid, gene_biotype) %>%
#   dplyr::filter(gene_id %in% group3$gene)
# 
# group4 <- clusters$df %>%
#   filter(cluster == 4)
# row.names(group4) <- 1:nrow(group4)
# group4_annotated <- genes(chick_ens, return.type = "data.frame")  %>%
#   dplyr::select(gene_id, gene_name, entrezid, gene_biotype) %>%
#   dplyr::filter(gene_id %in% group4$gene)


# 
# ###
# 
# # 
# # ## Gene ontology / hypergeometric testing... if it seems like it doesnt work its because there are no significant results
# # 
# Merge the AnnotationHub dataframe with the results

# # 
# # # Extract significant results
# # sigGO <- dplyr::filter(res_ids, padj < 0.05)
# # sigGO_500 <- sigGO %>%
# #   arrange(padj) %>%
# #   head(n=500)
# # sigGO_genes <- as.character(sigGO_500$gene)
# # 
# # # Run GO enrichment analysis, output a table (note that ont can be changed to reflect desired broad type of gene product function)
# # ego <- enrichGO(gene = sigGO_genes, universe = all_genes, keyType = "ENSEMBL", OrgDb = org.Gg.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE, pool = TRUE)
# # cluster_summary <- data.frame(ego)
# # 
# # dotplot(ego)
# # 
# # # GO enrichment but just with group 1 from previous clustering
# # group1_GO <- res_ids  %>%
# #   dplyr::filter(gene %in% group1$gene)
# # group1_GO_genes <- as.character(group1_GO$gene) 
# # ego1 <- enrichGO(gene = group1_GO_genes, universe = all_genes, keyType = "ENSEMBL", OrgDb = org.Gg.eg.db, ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE, pool = TRUE)
# # cluster_summary1 <- data.frame(ego1)
# # 
# # dotplot(ego1)
# # 
# # 
# # ###
# # 
# # 
# # ## GSEA
# Remove any NA values, entrez duplicates
# res_entrez <- dplyr::filter(res_ids_GG, entrezid != "NA")
# res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]
# res_entrez <- res_entrez  %>%
#   dplyr::filter(gene %in% sig_GGHo$gene)
# 
# # Extract the fold changes, sort in decreasing order
# foldchanges <- res_entrez$log2FoldChange
# names(foldchanges) <- res_entrez$entrezid
# foldchanges <- sort(foldchanges, decreasing = TRUE)
# 
# # GSEA using gene sets from KEGG pathways
# gseaKEGG <- gseKEGG(geneList = foldchanges,
#                     organism = "gga",
#                     pAdjustMethod = "none",
#                     minGSSize = 5, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
#                     pvalueCutoff = 0.05,
#                     verbose = FALSE)
# 
# # Extract the GSEA results
# gseaKEGG_results <- gseaKEGG@result
# 
# # Plot the GSEA plot for single enriched pathways
# gseaplot(gseaKEGG, geneSetID = 'pathway here')
# 
# # Use pathview
# Library(pathview)
# 
# # Output images for a single significant KEGG pathway
# pathview(gene.data = foldchanges,
#          pathway.id = "gga03010",
#          species = "gga",
#          limit = list(gene = 2, # value gives the max/min limit for foldchanges
#                       cpd = 1))
# 
# # GSEA using gene sets associated with BP Gene Ontology terms
# gseaGO <- gseGO(geneList = foldchanges,
#                 OrgDb = org.Gg.eg.db,
#                 ont = 'BP',
#                 minGSSize = 5,
#                 pvalueCutoff = 0.3,
#                 verbose = FALSE)
# 
# gseaGO_results <- gseaGO@result
# # # 
# # 
# # 
# # ###
# # 
# # 
# ## Volcano plot
# 
# sig_GGHo_e12_e14 <- read_excel("e12e14_venn_sig_midhindhoxd13.xlsx")
tgfbeta_ecm_genes <- read_excel("Tabin Lab/Molecular biology/rnaseq/tgfbeta_ecm_genes.xlsx")
tgfbeta_ecm_genes <- inner_join(tgfbeta_ecm_genes, sig_GGHo_e12_e14, by="gene_name.x.x") %>%
  dplyr::filter(log2FoldChange.x.x*log2FoldChange.y.x > 0) %>%
  arrange(-log2FoldChange.x.x)

tgfbeta_ecm_genes_2 <- read_excel("Tabin Lab/Molecular biology/rnaseq/tgfbeta_ecm_genes_2.xlsx")
tgfbeta_ecm_genes_2 <- inner_join(tgfbeta_ecm_genes_2, sig_GGHo_e14_unique, by="gene_name.x") %>%
  dplyr::filter(log2FoldChange.x*log2FoldChange.y > 0) %>%
  arrange(-log2FoldChange.x)

tgf_annot <- annotations %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::distinct()
tgfbeta_ecm_genes <- merge(tgfbeta_ecm_genes, tgf_annot, by.x="gene_name.x.x", by.y="gene_name")%>%
  arrange(-log2FoldChange.x.x)
# tgfbeta_ecm_genes <- merge(tgfbeta_ecm_genes, tgf_annot, by.x="gene_name.x", by.y="gene_name")%>%
#   arrange(-log2FoldChange.x)

res_table_tb_volc <- res_table_tb_GG %>%
  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.50)

# 
# ## Add all the gene symbols as a column from the galgal6 table using bind_cols()
galgal6annot <- annotations %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::distinct()
res_table_tb_volc <- bind_cols(res_table_tb_volc, gene_name=galgal6annot$gene_name[match(res_table_tb_volc$gene, galgal6annot$gene_id)])

# ## Create an empty column to indicate which genes to label
# res_table_tb_volc <- res_table_tb_volc %>% mutate(hoxhindsharedsymbol = "")
# res_table_tb_volc <- res_table_tb_volc %>% mutate(e12e14sharedcolor = "")
# res_table_tb_volc <- res_table_tb_volc %>% mutate(TGFgenelabels = "")


## Sort by padj values, filter out the top HOX genes for simplicity
res_table_tb_volc_e14_GHo <- res_table_tb_volc %>% 
  arrange(padj) %>%
  dplyr::filter(!gene %in% sig_GGHo$gene) %>%
  dplyr::filter(!gene %in% sig_GGHo_e12_e14$gene)
res_table_tb_volc_e14_GHo$padj <- -log10(res_table_tb_volc_e14_GHo$padj)
write.csv(res_table_tb_volc_e14_GHo, "volcano_e14_midhoxd13_only.csv")
# # filter(!gene_name %in% c("HOXB6", "HOXD13", "HOXB7", "HOXA10", "HOXA13", "HOXA3", "HOXD11", "HOXC6", "HOXA11", "HOXC9", "HOXD12", "HOXA9", "HOXD10", "HOXB3", "HOXD8", "HOXC8", "HOXB8", "HOXD9"))
# 
res_table_tb_volc_e14_GGHo <- res_table_tb_volc %>% 
  arrange(padj) %>%
  dplyr::filter(gene %in% sig_GGHo$gene) %>%
  dplyr::filter(!gene %in% sig_GGHo_e12_e14$gene) %>%
  dplyr::filter(!gene_name %in% tgfbeta_ecm_genes$gene_name.x.x)
res_table_tb_volc_e14_GGHo$padj <- -log10(res_table_tb_volc_e14_GGHo$padj)
write.csv(res_table_tb_volc_e14_GGHo, "volcano_e14_midhindhoxd13_noTGF_only.csv")

res_table_tb_volc_e12_GGHo_TGF <- res_table_tb_volc %>% 
  arrange(padj) %>%
  dplyr::filter(gene_name %in% tgfbeta_ecm_genes$gene_name.x.x)
res_table_tb_volc_e12_GGHo_TGF$padj <- -log10(res_table_tb_volc_e12_GGHo_TGF$padj)
# write.csv(res_table_tb_volc_e14_GGHo_TGF, "volcano_e12_midhindhoxd13_TGF.csv")

res_table_tb_volc_e12e14_GGHo <- res_table_tb_volc %>% 
  arrange(padj) %>%
  dplyr::filter(gene %in% sig_GGHo_e12_e14$gene) %>%
  dplyr::filter(!gene_name %in% tgfbeta_ecm_genes$gene_name.x.x)
res_table_tb_volc_e12e14_GGHo$padj <- -log10(res_table_tb_volc_e12e14_GGHo$padj)
write.csv(res_table_tb_volc_e12e14_GGHo, "volcano_e12e14_midhindhoxd13_noTGF.csv")

res_table_tb_volc_e12_all <- rbind(res_table_tb_volc_e14_GHo, res_table_tb_volc_e14_GGHo, res_table_tb_volc_e12e14_GGHo, res_table_tb_volc_e14_GGHo_TGF)
write.csv(res_table_tb_volc_e12_all, "volcano_e14_ordered.csv")

e14_TGF_unique <- tgfbeta_ecm_genes2 %>%
  dplyr::filter(!gene_name.x %in% tgfbeta_ecm_genes$gene_name.x.x)
# # Populate the genelabels column with contents of the gene symbols column for the first 10 rows, i.e. the top 10 most significantly expressed genes
# # res_table_tb_volc$genelabels[volcnum] <- as.character(res_table_tb_volc$gene_name[volcnum])
# 
# ## Label genes in common with other comparisons and/or timepoints
# res_table_tb_volc <- res_table_tb_volc %>% mutate(genelabels_e12_unique = "")
# res_table_tb_volc <- res_table_tb_volc %>% mutate(genelabels_e12_e14_venn = "")
res_table_tb_volc <- res_table_tb_volc %>% mutate(genelabels_venn = "")
res_table_tb_volc <- res_table_tb_volc %>% mutate(genelabels_venn_TF = "")
# 
# 
# # DEGs in common between mGFP midgut vs HoxD13 midgut & mGFP hindgut vs mGFP midgut
volc_num <- NULL
sig_list <- read_excel("Tabin Lab/Molecular biology/rnaseq/rnaseq output/sig_venn_midcontrolhindcontrol_midcontrolmidHOX_e14_-c(4,6,12).xlsx")
x <- c(1:nrow(sig_list))
for (val in x) {
  volc_num[val] <- rownames(res_table_tb_volc)[res_table_tb_volc$gene_name == sig_list$gene_name[val]]
}
volc_num <- as.numeric(volc_num)
res_table_tb_volc$genelabels_venn[volc_num] <- as.character(res_table_tb_volc$gene_name[volc_num])
zeros <- rep(0,17136)
ones <- rep(1,101)
res_table_tb_volc$genelabels_venn_TF <- zeros
res_table_tb_volc$genelabels_venn_TF[volc_num] <- ones
res_table_tb_volc$genelabels_venn_TF <- as.logical(res_table_tb_volc$genelabels_venn_TF)

# # DEGs unique to E14 in common between mGFP midgut vs HoxD13 midgut and mGFP hindgut vs mGFP midgut
# volc_num <- NULL
# sig_list <- read_excel("Tabin Lab/Molecular biology/rnaseq/rnaseq output/sig_unique_e12.xlsx")
# x <- c(1:nrow(sig_list))
# for (val in x) {
#   volc_num[val] <- rownames(res_table_tb_volc)[res_table_tb_volc$gene_name == sig_list$gene_name[val]]
# }
# volc_num <- as.numeric(volc_num)
# res_table_tb_volc$genelabels_e12_unique[volc_num] <- as.character(res_table_tb_volc$gene_name[volc_num])
# 
# # DEGs at both E12 and E14 in common between mGFP midgut vs HoxD13 midgut & mGFP hindgut vs mGFP midgut
# volc_num <- NULL
# sig_list <- read_excel("Tabin Lab/Molecular biology/rnaseq/rnaseq output/sig_venn_e12_e14.xlsx")
# x <- c(1:nrow(sig_list))
# for (val in x) {
#   volc_num[val] <- rownames(res_table_tb_volc)[res_table_tb_volc$gene_name == sig_list$gene_name[val]]
# }
# volc_num <- as.numeric(volc_num)
# res_table_tb_volc$genelabels_e12_e14_venn[volc_num] <- as.character(res_table_tb_volc$gene_name[volc_num])
# 
# 
# # Plot with different sets labeled
# pdf(file="volcano_e12_midcontrol_midhoxd13_e12_venn_labels.pdf")
ggplot(res_table_tb_volc, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = genelabels_venn_TF)) +
  geom_point(aes(colour = threshold_OE)) +
  geom_text_repel(aes(label = genelabels_venn), size = rel(5)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = rel(1.5), hjust = 0.5),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(2))) +
  ylim(-0.2,18) +
  xlim(-4.5,4)

# dev.off()
# #
# pdf(file="volcano_e12_midcontrol_midhoxd13_e12_unique_labels.pdf")
# ggplot(res_table_tb_volc, aes(x = log2FoldChange, y = -log10(padj))) +
#   geom_point(aes(colour = threshold_OE)) +
#   geom_text_repel(aes(label = genelabels_e12_unique), max.overlaps = 50) +
#   ggtitle(label = "E12 RCAS-HOXD13 midgut vs control midgut", subtitle = "Genes unique to E12 in common with control midgut vs hindgut") +
#   xlab("log2 fold change") +
#   ylab("-log10 adjusted p-value") +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.5), hjust = 0.5),
#         plot.subtitle = element_text(size = rel(1.25), hjust = 0.5),
#         axis.title = element_text(size = rel(1.25))) +
#   ylim(-0.2,18) +
#   xlim(-4,6)
# dev.off()
#
# pdf(file="volcano_e12_midcontrol_midhoxd13_e12_e14_venn_labels.pdf")
# ggplot(res_table_tb_volc, aes(x = log2FoldChange, y = -log10(padj))) +
#   geom_point(aes(colour = threshold_OE)) +
#   geom_text_repel(aes(label = genelabels_e12_e14_venn), max.overlaps = 100) +
#   ggtitle(label = "E12 RCAS-HOXD13 midgut vs control midgut", subtitle = "Genes in common at E14 with control midgut vs hindgut") +
#   xlab("log2 fold change") +
#   ylab("-log10 adjusted p-value") +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.5), hjust = 0.5),
#         plot.subtitle = element_text(size = rel(1.25), hjust = 0.5),
#         axis.title = element_text(size = rel(1.25))) +
#   ylim(-0.2,18) +
#   xlim(-4,6)
# dev.off()

### 


# ## Plot expression of only one gene
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

# 
# 
# ###
# 
# 
# ## Output
# 


# 
# #get top 10 most significant genes
# sig_padj <- sig %>%
#   arrange(padj)
# sig_padj <- head(sig_padj, 30)
# # write.csv(sig_padj, "./e12_output/midcontrol_midhoxd13/sig_padj_e12_midcontrol_midhoxd13.csv", row.names = FALSE)
# # 
# # #get top 10 most positive LFC
# sig_lfc_pos <- sig %>%
#   arrange(log2FoldChange)
# sig_lfc_pos <- tail(sig_lfc_pos, 30)
# # write.csv(sig_lfc_pos, "./e12_output/midcontrol_midhoxd13/sig_lfc_pos_e12_midcontrol_midhoxd13.csv", row.names = FALSE)
# # 
# # #get top ten most negative LFC
# sig_lfc_neg <- sig %>%
#   arrange(log2FoldChange)
# sig_lfc_neg <- head(sig_lfc_neg, 10)
# # write.csv(sig_lfc_neg, "./e12_output/midcontrol_midhoxd13/sig_lfc_neg_e12_midcontrol_midhoxd13.csv", row.names = FALSE)
# # 
# # #get clusters
# sig_group1 <- annotations  %>%
#   dplyr::select(gene_id, gene_name) %>%
#   dplyr::distinct() %>%
#   dplyr::filter(gene_id %in% group1$gene)
# sig_group1 <- merge(res_table_tb, group1_annotated, by.x="gene", by.y="gene_id")
# # sig_group1 <- sig_group1 %>% 
# #   select(-(8:9))
# # write.csv(sig_group1, "./e12_output/midcontrol_midhoxd13/sig_group1_e12_midcontrol_midhoxd13.csv", row.names = FALSE)
# # 
# sig_group2 <- annotations  %>%
#   dplyr::select(gene_id, gene_name) %>%
#   dplyr::distinct() %>%
#   dplyr::filter(gene_id %in% group2$gene)
# sig_group2 <- merge(res_table_tb, group2_annotated, by.x="gene", by.y="gene_id")
# # sig_group2 <- sig_group2 %>% 
# #   select(-(8:9))
# # write.csv(sig_group2, "./e12_output/midcontrol_midhoxd13/sig_group2_e12_midcontrol_midhoxd13.csv", row.names = FALSE)
# # 
# sig_group3 <- annotations  %>%
#   dplyr::select(gene_id, gene_name) %>%
#   dplyr::distinct() %>%
#   dplyr::filter(gene_id %in% group3$gene)
# sig_group3 <- merge(res_table_tb, group3_annotated, by.x="gene", by.y="gene_id")
# # sig_group3 <- sig_group3 %>% 
# #   select(-(8:9))
# # write.csv(sig_group3, "./e12_output/midcontrol_midhoxd13/sig_group3_e12_midcontrol_midhoxd13.csv", row.names = FALSE)
# # 
# sig_group4 <- annotations  %>%
#   dplyr::select(gene_id, gene_name) %>%
#   dplyr::distinct() %>%
#   dplyr::filter(gene_id %in% group4$gene)
# sig_group4 <- merge(res_table_tb, group4_annotated, by.x="gene", by.y="gene_id")
# sig_group4 <- sig_group4 %>%
# #   select(-(8:9))
# # write.csv(sig_group4, "./e12_output/midcontrol_midhoxd13/sig_group4_e12_midcontrol_midhoxd13.csv", row.names = FALSE)
# # 
# # #save images
# pdf(file="clusters_e12")
# degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col=NULL)
# dev.off()


# Gene expression heatmap - HG 11/7/21

# clusters_heatmap <- bind_rows(group1, group2, group3, group4)
# clusters_heatmap <- clusters_heatmap[clusters_heatmap$genes %in% sig_GG$gene,]

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
# test <- NULL
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

# grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))
# 
# sampleDendrogram_1 <- sampleDendrogram + scale_x_continuous(expand=c(.0085, .0085)) + scale_y_continuous(expand=c(0, 0))
# heatmap_1 <- heatmap + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0))
# 
# # Convert both grid based objects to grobs
# sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
# heatmapGrob <- ggplotGrob(heatmap)
# 
# # Check the widths of each grob
# sampleDendrogramGrob$widths
# heatmapGrob$widths
# 
# # Add in the missing columns
# sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
# sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)
# 
# # Make sure every width between the two grobs is the same
# maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths)
# sampleDendrogramGrob$widths <- as.list(maxWidth)
# heatmapGrob$widths <- as.list(maxWidth)
# 
# # Arrange the grobs into a plot
# finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))
# 
# # Draw the plot
# grid.draw(finalGrob)
# # # Convert the VST counts to long format for ggplot2
# # 
# # 
# # # First compare wide vs long version
# # deseq2VST_wide <- deseq2VST
# # deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))
# # 
# # head(deseq2VST_wide)
# # head(deseq2VST_long)
# 
# 
# # 
# # Convert the significant genes back to a matrix for clustering
# 
# 
# 
# # Compute a distance calculation on both dimensions of the matrix
# distanceGene <- dist(deseq2VSTMatrix)
# distanceSample <- dist(t(deseq2VSTMatrix))
# 
# # Cluster based on the distance calculations
# clusterGene <- hclust(distanceGene, method="average")
# clusterSample <- hclust(distanceSample, method="average")
# 
# # load the preliminary data into R
# # load(url("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/exercise_1/deseq2_exercise1_env.RData"))
# 
# ################################################################################
# ################# Step 1: create dendrogram for genes ##########################
# 
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
tgfe14unique <- tgfbeta_ecm_genes_2$gene_name.x

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements
# keyvals.shape <- ifelse(
#   res$gene_name %in% hindhox, 0,
#   ifelse(res$gene_name %in% tgfe12e14, 1,
#          2))
# keyvals.shape[is.na(keyvals.shape)] <- 3
# names(keyvals.shape)[keyvals.shape == 0] <- 'shared hind & hox'
# names(keyvals.shape)[keyvals.shape == 1] <- 'e12 e14 tgf'
# names(keyvals.shape)[keyvals.shape == 2] <- 'e12 e14'

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