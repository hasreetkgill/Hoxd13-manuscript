# Load libraries
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")
library("dplyr")
library("ggrepel")
library(dplyr)

ColorMap <- brewer.pal(11, "PRGn")

# Load Data
read.counts.input<-read.table("Hasreet_ReorderedLightSeq_MidHindHoxD13_05.07.22.csv", header=TRUE, sep=',', stringsAsFactors=FALSE)
row.names(read.counts.input)<-read.counts.input[,1]
read.counts.input<-read.counts.input[, -c(0:1)]

read.counts.innermes<-read.counts.input[, c(4,7,10,2,5,8,11,3,6,9,12)] #Inner mesenchyme only
read.counts.muscmuco<-read.counts.input[, c(16,19,22,14,17,20,23,15,18,21,24)] #Muscularis mucosa only
read.counts.innercirc<-read.counts.input[, c(28,31,34,26,29,32,35,27,30,33,36)] #Inner circular muscle only
read.counts.hind<-read.counts.input[, c(1,4,7,10,13,16,19,22,25,28,31,34)] #Hindgut only
read.counts.hox<-read.counts.input[, c(2,5,8,11,14,17,20,23,26,29,32,35)] #HoxD13 only
read.counts.mid<-read.counts.input[, c(3,6,9,12,15,18,21,24,27,30,33,36)] #Midgut only

sample.info.gen<-data.frame(chick=c("2","3","4","1","2","3","4","1","2","3","4"), condition=c(rep("Hindgut",3), rep("HoxD13+ Midgut",4),
                                                                                              rep("Midgut",4)), row.names=names(read.counts.muscmuco)) #Comparing genotypes
# sample.info.reg<-data.frame(chick=c("1","2","3","4","1","2","3","4","1","2","3","4"), condition=c(rep("Inner Mesenchyme",4), rep("Muscularis Mucosa",4),
#                                                                                               rep("Inner Circular",4)), row.names=names(read.counts.hind)) #Comparing regions
# Differential Expression analysis
# Create DESeq object
DESeq.ds<-DESeqDataSetFromMatrix(countData=read.counts.muscmuco, colData=sample.info.gen, design = ~ condition) #change which read.counts & colData
DESeq.ds <- DESeq(DESeq.ds)

## Comparing genotypes
results.mes.muco <- results(DESeq.ds, pAdjustMethod="BH", contrast = c("condition", "Inner Mesenchyme","Muscularis Mucosa"))
results.mes.circ <- results(DESeq.ds, pAdjustMethod="BH", contrast = c("condition", "Inner Mesenchyme","Inner Circular"))
results.muco.circ <- results(DESeq.ds, pAdjustMethod="BH", contrast = c("condition", "Muscularis Mucosa","Inner Circular"))
DGE.results <- c(results.mes.muco, results.mes.circ, results.muco.circ)

table(results.mes.muco$padj<0.05)
table(results.mes.circ$padj<0.05)
table(results.muco.circ$padj<0.05)

# Sort and obtain differentially expressed genes in a csv file
results.mes.muco.sorted <- results.mes.muco[order(results.mes.muco$log2FoldChange, decreasing=TRUE),]
results.mes.circ.sorted <- results.mes.circ[order(results.mes.circ$log2FoldChange, decreasing=TRUE),]
results.muco.circ.sorted <- results.muco.circ[order(results.muco.circ$log2FoldChange, decreasing=TRUE),]
# results.mes.muco.sorted <- results.mes.muco[order(results.mes.muco$padj),]
# results.mes.circ.sorted <- results.mes.circ[order(results.mes.circ$padj),]
# results.muco.circ.sorted <- results.muco.circ[order(results.muco.circ$padj),]

DGEgenes.mes.muco <- rownames(subset(results.mes.muco.sorted, padj<0.05))
DGEgenes.mes.circ <- rownames(subset(results.mes.circ.sorted, padj<0.05))
DGEgenes.muco.circ <- rownames(subset(results.muco.circ.sorted, padj<0.05))
All.DGEgenes <- c(DGEgenes.mes.muco, DGEgenes.mes.circ, DGEgenes.muco.circ)

DE_genes.mes.muco <- as.data.frame(results.mes.muco.sorted)
DE_genes.mes.circ <- as.data.frame(results.mes.circ.sorted)
DE_genes.muco.circ <- as.data.frame(results.muco.circ.sorted)
write.csv(DE_genes.mes.muco, "Midgut_InnerMesvsMuscMuco_byLFC.csv")
write.csv(DE_genes.mes.circ, "HMidgut_InnerMesvsInnerCirc_byLFC.csv")
write.csv(DE_genes.muco.circ, "Midgut_MuscMucovsInnerCirc_byLFC.csv")

#DE genes different between pairs of genotypes
mes.muco.only.Pos <- rownames(subset(results.mes.muco.sorted, log2FoldChange>0&padj<0.05))
mes.muco.only.Neg <- rownames(subset(results.mes.muco.sorted, log2FoldChange<0&padj<0.05))
mes.circ.only.Pos <- rownames(subset(results.mes.circ.sorted, log2FoldChange>0&padj<0.05))
mes.circ.only.Neg <- rownames(subset(results.mes.circ.sorted, log2FoldChange<0&padj<0.05))
muco.circ.only.Pos <- rownames(subset(results.muco.circ.sorted, log2FoldChange>0&padj<0.05))
muco.circ.only.Neg <- rownames(subset(results.muco.circ.sorted, log2FoldChange<0&padj<0.05))

#DE genes that are specific to each region
mes.Pos <- intersect(rownames(subset(results.mes.muco.sorted, log2FoldChange>0&padj<0.1)), rownames(subset(results.mes.muco.sorted, log2FoldChange>0&padj<0.1)))
mes.Neg <- intersect(rownames(subset(results.mes.muco.sorted, log2FoldChange>0&padj<0.1)), rownames(subset(results.mes.muco.sorted, log2FoldChange<0&padj<0.1)))

muco.Pos <- intersect(rownames(subset(results.mes.muco.sorted, log2FoldChange<0&padj<0.1)), rownames(subset(results.muco.circ.sorted, log2FoldChange>0&padj<0.1)))
muco.Neg <- intersect(rownames(subset(results.mes.muco.sorted, log2FoldChange>0&padj<0.1)), rownames(subset(results.muco.circ.sorted, log2FoldChange<0&padj<0.1)))

circ.Pos <- intersect(rownames(subset(results.mes.circ.sorted, log2FoldChange<0&padj<0.1)), rownames(subset(results.muco.circ.sorted, log2FoldChange<0&padj<0.1)))
circ.Neg <- intersect(rownames(subset(results.mes.circ.sorted, log2FoldChange>0&padj<0.1)), rownames(subset(results.muco.circ.sorted, log2FoldChange>0&padj<0.1)))

write.csv(as.data.frame(mes.Pos), "Midgut_InnerMesmarkers_byLFC.csv")
write.csv(as.data.frame(muco.Pos), "Midgut_MuscMucomarkers_byLFC.csv")
write.csv(as.data.frame(circ.Pos), "Midgut_InnerCircmarkers_byLFC.csv")

# #DE genes shared between pairs of genotypes
# hind.hox.Pos <- intersect(rownames(subset(results.hind.mid.sorted, log2FoldChange>0&padj<0.1)), rownames(subset(results.hox.mid.sorted, log2FoldChange>0&padj<0.1)))
# hind.hox.Neg <- intersect(rownames(subset(results.hind.mid.sorted, log2FoldChange<0&padj<0.1)), rownames(subset(results.hox.mid.sorted, log2FoldChange<0&padj<0.1)))
# hox.mid.Pos <- intersect(rownames(subset(results.hind.hox.sorted, log2FoldChange<0&padj<0.1)), rownames(subset(results.hind.mid.sorted, log2FoldChange<0&padj<0.1)))
# hox.mid.Neg <- intersect(rownames(subset(results.hind.hox.sorted, log2FoldChange>0&padj<0.1)), rownames(subset(results.hind.mid.sorted, log2FoldChange>0&padj<0.1)))
# hind.mid.Pos <- intersect(rownames(subset(results.hind.hox.sorted, log2FoldChange>0&padj<0.1)), rownames(subset(results.hox.mid.sorted, log2FoldChange<0&padj<0.1)))
# hind.mid.Neg <- intersect(rownames(subset(results.hind.hox.sorted, log2FoldChange<0&padj<0.1)), rownames(subset(results.hox.mid.sorted, log2FoldChange>0&padj<0.1)))
# 
# write.csv(as.data.frame(hind.hox.Pos), "InnerCirc_HindHoxShared_Up_byPadj.csv")
# write.csv(as.data.frame(hox.mid.Pos), "InnerCirc_HoxMidShared_Up_byPadj.csv")
# write.csv(as.data.frame(hind.mid.Pos), "InnerCirc_HindMidShared_Up_byPadj.csv")
# write.csv(as.data.frame(hind.hox.Neg), "InnerCirc_HindHoxShared_Down_byPadj.csv")
# write.csv(as.data.frame(hox.mid.Neg), "InnerCirc_HoxMidShared_Down_byPadj.csv")
# write.csv(as.data.frame(hind.mid.Neg), "InnerCirc_HindMidShared_Down_byPadj.csv")

rld <- rlog(DESeq.ds, blind = TRUE, fitType='local')
log.norm.counts <- assay(rld)

# Heatmap plot of all DGEs
DGE_Top<-unique(c(mes.muco.only.Pos, mes.muco.only.Neg, mes.circ.only.Pos, mes.circ.only.Neg,muco.circ.only.Pos, muco.circ.only.Neg))
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[-grep("ENSGAL", rownames(hm.mat_DGEgenes)), ]
colnames(hm.mat_DGEgenes)<-c("Mes_1","Mes_2", "Mes_3", "Mes_4","Muco_1","Muco_2", "Muco_3", "Muco_4","Circ_1","Circ_2", "Circ_3", "Circ_4")
pdf(file="Midgut_Top_DEGenesbyLFC_All_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap, fontsize=5)
dev.off()

# Heatmap plot of pairwise comparisons
DGE_Top<-c(mes.muco.only.Pos, mes.muco.only.Neg)
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,c(1:8)]
hm.mat_DGEgenes<-hm.mat_DGEgenes[-grep("ENSGAL", rownames(hm.mat_DGEgenes)), ]
colnames(hm.mat_DGEgenes)<-c("Mes_1", "Mes_2", "Mes_3", "Mes_4", "Muco_1", "Muco_2", "Muco_3", "Muco_4")
pdf(file="Midgut_Top_DEGenesbyLFC_InnerMesvMuscMuco_Pos&Neg_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

DGE_Top<-c(mes.circ.only.Pos, mes.circ.only.Neg)
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,c(1:4,9:12)]
hm.mat_DGEgenes<-hm.mat_DGEgenes[-grep("ENSGAL", rownames(hm.mat_DGEgenes)), ]
colnames(hm.mat_DGEgenes)<-c("Mes_1","Mes_2", "Mes_3", "Mes_4", "Circ_1", "Circ_2", "Circ_3", "Circ_4")
pdf(file="Midgut_Top_DEGenesbyLFC_InnerMesvsInnerCirc_Pos&Neg_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap, fontsize=4)
dev.off()

DGE_Top<-c(muco.circ.only.Pos, muco.circ.only.Neg)
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,c(5:12)]
hm.mat_DGEgenes<-hm.mat_DGEgenes[-grep("ENSGAL", rownames(hm.mat_DGEgenes)), ]
colnames(hm.mat_DGEgenes)<-c("Mucu_1","Muco_2", "Muco_3", "Muco_4", "Circ_1","Circ_2", "Circ_3", "Circ_4")
pdf(file="Midgut_Top_DEGenesbyPadj_MuscMucovsInnerCirc_Pos&Neg_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

# Heatmap plots of cell type markers
DGE_Top<-c(mes.Pos, muco.Pos, circ.Pos)
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[-grep("ENSGAL", rownames(hm.mat_DGEgenes)), ]
colnames(hm.mat_DGEgenes)<-c("Mes_1","Mes_2", "Mes_3", "Mes_4","Muco_1","Muco_2", "Muco_3", "Muco_4", "Circ_1","Circ_2", "Circ_3", "Circ_4")
pdf(file="Midgut_Top_DEGenesbyLFC_AllLayers_Pos_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()

DGE_Top<-c(mes.Neg, muco.Neg, circ.Neg)
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[-grep("ENSGAL", rownames(hm.mat_DGEgenes)), ]
colnames(hm.mat_DGEgenes)<-c("Mes_2", "Mes_3", "Mes_4","Muco_2", "Muco_3", "Muco_4", "Circ_2", "Circ_3", "Circ_4")
pdf(file="Hindgut_Top_DEGenesbyPadj_AllLayers_Neg_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
dev.off()


# ## Comparing genotypes
results.hind.hox <- results(DESeq.ds, pAdjustMethod="BH", contrast = c("condition", "Hindgut","HoxD13+ Midgut"))
results.hind.mid <- results(DESeq.ds, pAdjustMethod="BH", contrast = c("condition", "Hindgut","Midgut"))
results.hox.mid <- results(DESeq.ds, pAdjustMethod="BH", contrast = c("condition", "HoxD13+ Midgut","Midgut"))
DGE.results <- c(results.hind.hox, results.hind.mid, results.hox.mid)

# table(results.hind.hox$padj<0.1)
# table(results.hind.mid$padj<0.1)
# table(results.hox.mid$padj<0.1)
# 
# # Sort and obtain differentially expressed genes in a csv file
# # results.hind.hox.sorted <- results.hind.hox[order(results.hind.hox$log2FoldChange, decreasing=TRUE),]
# # results.hind.mid.sorted <- results.hind.mid[order(results.hind.mid$log2FoldChange, decreasing=TRUE),]
# # results.hox.mid.sorted <- results.hox.mid[order(results.hox.mid$log2FoldChange, decreasing=TRUE),]
results.hind.hox.sorted <- results.hind.hox[order(results.hind.hox$padj),]
results.hind.mid.sorted <- results.hind.mid[order(results.hind.mid$padj),]
results.hox.mid.sorted <- results.hox.mid[order(results.hox.mid$padj),]

DGEgenes.hind.hox <- rownames(subset(results.hind.hox.sorted, padj<0.01))
DGEgenes.hind.mid <- rownames(subset(results.hind.mid.sorted, padj<0.01))
DGEgenes.hox.mid <- rownames(subset(results.hox.mid.sorted, padj<0.01))
All.DGEgenes <- c(DGEgenes.hind.hox, DGEgenes.hind.mid, DGEgenes.hox.mid)
# 
# DE_genes.hind.hox <- as.data.frame(results.hind.hox.sorted)
# DE_genes.hind.mid <- as.data.frame(results.hind.mid.sorted)
# DE_genes.hox.mid <- as.data.frame(results.hox.mid.sorted)
# write.csv(DE_genes.hind.hox, "InnerCirc_HindvsHox_byPadj.csv")
# write.csv(DE_genes.hind.mid, "InnerCirc_HindvsMid_byPadj.csv")
# write.csv(DE_genes.hox.mid, "InnerCirc_HoxvsMid_byPadj.csv")
# 
# #DE genes different between pairs of genotypes
hind.hox.only.Pos <- rownames(subset(results.hind.hox.sorted, log2FoldChange>0&padj<0.01))
hind.hox.only.Neg <- rownames(subset(results.hind.hox.sorted, log2FoldChange<0&padj<0.01))
hind.mid.only.Pos <- rownames(subset(results.hind.mid.sorted, log2FoldChange>0&padj<0.01))
hind.mid.only.Neg <- rownames(subset(results.hind.mid.sorted, log2FoldChange<0&padj<0.01))
hox.mid.only.Pos <- rownames(subset(results.hox.mid.sorted, log2FoldChange>0&padj<0.01))
hox.mid.only.Neg <- rownames(subset(results.hox.mid.sorted, log2FoldChange<0&padj<0.01))

# #DE genes that are specific to each genotype
# # hind.Pos <- intersect(rownames(subset(results.hind.hox.sorted, log2FoldChange>0&padj<0.1)), rownames(subset(results.hind.mid.sorted, log2FoldChange>0&padj<0.1)))
# # hind.Neg <- intersect(rownames(subset(results.hind.hox.sorted, log2FoldChange>0&padj<0.1)), rownames(subset(results.hind.mid.sorted, log2FoldChange<0&padj<0.1)))
# # 
# # hox.Pos <- intersect(rownames(subset(results.hind.hox.sorted, log2FoldChange<0&padj<0.1)), rownames(subset(results.hox.mid.sorted, log2FoldChange>0&padj<0.1)))
# # hox.Neg <- intersect(rownames(subset(results.hind.hox.sorted, log2FoldChange>0&padj<0.1)), rownames(subset(results.hox.mid.sorted, log2FoldChange<0&padj<0.1)))
# # 
# # mid.Pos <- intersect(rownames(subset(results.hind.mid.sorted, log2FoldChange<0&padj<0.1)), rownames(subset(results.hox.mid.sorted, log2FoldChange<0&padj<0.1)))
# # mid.Neg <- intersect(rownames(subset(results.hind.mid.sorted, log2FoldChange>0&padj<0.1)), rownames(subset(results.hox.mid.sorted, log2FoldChange>0&padj<0.1)))
# # 
# # write.csv(as.data.frame(hind.Pos), "InnerMes_Hindmarkers.csv")
# # write.csv(as.data.frame(hox.Pos), "InnerMes_Hoxmarkers.csv")
# # write.csv(as.data.frame(mid.Pos), "InnerMes_Midmarkers.csv")
# 
# #DE genes shared between pairs of genotypes
hind.hox.Pos <- intersect(rownames(subset(results.hind.mid.sorted, log2FoldChange>0&padj<0.01)), rownames(subset(results.hox.mid.sorted, log2FoldChange>0&padj<0.01)))
hind.hox.Neg <- intersect(rownames(subset(results.hind.mid.sorted, log2FoldChange<0&padj<0.01)), rownames(subset(results.hox.mid.sorted, log2FoldChange<0&padj<0.01)))
hox.mid.Pos <- intersect(rownames(subset(results.hind.hox.sorted, log2FoldChange<0&padj<0.01)), rownames(subset(results.hind.mid.sorted, log2FoldChange<0&padj<0.01)))
hox.mid.Neg <- intersect(rownames(subset(results.hind.hox.sorted, log2FoldChange>0&padj<0.01)), rownames(subset(results.hind.mid.sorted, log2FoldChange>0&padj<0.01)))
hind.mid.Pos <- intersect(rownames(subset(results.hind.hox.sorted, log2FoldChange>0&padj<0.01)), rownames(subset(results.hox.mid.sorted, log2FoldChange<0&padj<0.01)))
hind.mid.Neg <- intersect(rownames(subset(results.hind.hox.sorted, log2FoldChange<0&padj<0.01)), rownames(subset(results.hox.mid.sorted, log2FoldChange>0&padj<0.01)))

write.csv(as.data.frame(hind.hox.Pos), "MuscMuco_HindHoxShared_Up_byPadj.csv")
write.csv(as.data.frame(hox.mid.Pos), "MuscMuco_HoxMidShared_Up_byPadj.csv")
write.csv(as.data.frame(hind.mid.Pos), "MuscMuco_HindMidShared_Up_byPadj.csv")
write.csv(as.data.frame(hind.hox.Neg), "MuscMuco_HindHoxShared_Down_byPadj.csv")
write.csv(as.data.frame(hox.mid.Neg), "MuscMuco_HoxMidShared_Down_byPadj.csv")
write.csv(as.data.frame(hind.mid.Neg), "MuscMuco_HindMidShared_Down_byPadj.csv")

rld <- rlog(DESeq.ds, blind = TRUE, fitType='local')
log.norm.counts <- assay(rld)
# 
# # Heatmap plot of pairwise comparisons
DGE_Top<-c(hind.mid.only.Pos[1:30], hind.mid.only.Neg[1:30])
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[-grep("ENSGAL", rownames(hm.mat_DGEgenes)), ]
colnames(hm.mat_DGEgenes)<-c("Hind_2", "Hind_3", "Hind_4", "Mid_1", "Mid_2", "Mid_3", "Mid_4")
pdf(file="InnerMes_Top_DEGenesbyPadj_HindvMid_Pos&Neg_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap, fontsize=7, cellwidth=12)
dev.off()
# 
DGE_Top<-c(hox.mid.only.Pos, hox.mid.only.Neg)
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[-grep("ENSGAL", rownames(hm.mat_DGEgenes)), ]
colnames(hm.mat_DGEgenes)<-c("Hox_1", "Hox_2", "Hox_3", "Hox_4", "Mid_1", "Mid_2", "Mid_3", "Mid_4")
pdf(file="MuscMuco_Top_DEGenesbyPadj_HoxvMid_Pos&Neg_Heatmap.pdf", onefile=FALSE)
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap, fontsize=12, cellwidth=15, cellheight=8)
dev.off()

# 
# # Heatmap plots of shared up & down genes between relevant genotype pairs
DGE_Top<-c(hind.hox.Pos[1:20], hind.hox.Neg[1:20])
hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
hm.mat_DGEgenes<-hm.mat_DGEgenes[-grep("ENSGAL", rownames(hm.mat_DGEgenes)), ]
colnames(hm.mat_DGEgenes)<-c("Hind_2", "Hind_3", "Hind_4", "Hox_1", "Hox_2", "Hox_3", "Hox_4", "Mid_1", "Mid_2", "Mid_3", "Mid_4")
# pdf(file="InnerMes_Top_SharedDEGenesbyPadj_Pos&Neg_HindHox_Heatmap.pdf", onefile=FALSE)
tiff('InnerMes_Top_SharedDEGenesbyPadj_Pos&Neg_HindHox_Heatmap.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap,fontsize=9, cellwidth=15)
dev.off()

# # Heatmap plot of top & bottom 20 differentially expressed genes for each population
# DGE_Top<-c(hind.Pos,hind.Neg)
# hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
# colnames(hm.mat_DGEgenes)<-c("Hind_1", "Hind_2", "Hind_3", "Hind_4", "Hox_1", "Hox_2", "Hox_3", "Hox_4", "Mid_1", "Mid_2", "Mid_3", "Mid_4")
# pdf(file="Top_DEGenes_Hind_Pos&Neg_Heatmap.pdf", onefile=FALSE)
# pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
# dev.off()
# 
# DGE_Top<-c(hox.Pos,hox.Neg)
# hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
# colnames(hm.mat_DGEgenes)<-c("Hind_1", "Hind_2", "Hind_3", "Hind_4", "Hox_1", "Hox_2", "Hox_3", "Hox_4", "Mid_1", "Mid_2", "Mid_3", "Mid_4")
# pdf(file="Top_DEGenes_Hox_Pos&Neg_Heatmap.pdf", onefile=FALSE)
# pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
# dev.off()
# 
# DGE_Top<-c(mid.Pos,mid.Neg)
# hm.mat_DGEgenes<-log.norm.counts[DGE_Top,]
# colnames(hm.mat_DGEgenes)<-c("Hind_1", "Hind_2", "Hind_3", "Hind_4", "Hox_1", "Hox_2", "Hox_3", "Hox_4", "Mid_1", "Mid_2", "Mid_3", "Mid_4")
# pdf(file="Top_DEGenes_Mid_Pos&Neg_Heatmap.pdf", onefile=FALSE)
# pheatmap(hm.mat_DGEgenes, cluster_rows=F, cluster_cols=F, scale="row", col = ColorMap)
# dev.off()
# 
# 

