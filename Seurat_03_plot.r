#Seurat3.1.1 Data visualization

setwd("H:/???ཿ ??ɽҽԺ/code")
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(plyr)
library(tidyr)
library(gridExtra)
library(monocle)
library(Seurat)
library(pheatmap)
library(Matrix)
library(MAST)
library(ggrepel)



########################### HeatMap for clusters
# read data
immune.combined <- readRDS(file = "H:/???ཿ ??ɽҽԺ/1/immune.combined.three.rds")
DefaultAssay(immune.combined) <- "RNA"
#cluster name
levels(immune.combined$seurat_clusters)
levels(immune.combined$seurat_clusters)[3] <- "RPLs"
levels(immune.combined$seurat_clusters)[1] <- "FGFBP2"
levels(immune.combined$seurat_clusters)[4] <- "TNFRSF4"
levels(immune.combined$seurat_clusters)[9] <- "NKG7"
levels(immune.combined$seurat_clusters)[12] <- "RPLs"
levels(immune.combined$seurat_clusters)[10] <- "NKG7"
levels(immune.combined$seurat_clusters)[7] <- "FGFBP2"
levels(immune.combined$seurat_clusters)[6] <- "TNFRSF4"
levels(immune.combined$seurat_clusters)[2] <- "DUSP2"
levels(immune.combined$seurat_clusters)[5] <- "CCR7"
levels(immune.combined$seurat_clusters)[6] <- "ISG15"
levels(immune.combined$seurat_clusters)[8] <- "MAP3K8"
levels(immune.combined$seurat_clusters)[9] <- "GZMK"
levels(immune.combined$seurat_clusters)[10] <- "PRF1"
levels(immune.combined$seurat_clusters)[11] <- "F5"
levels(immune.combined$seurat_clusters)[12] <- "CD79A"
levels(immune.combined$seurat_clusters)[13] <- "LYZ"
levels(immune.combined$seurat_clusters)[14] <- "PPBP"
levels(immune.combined$seurat_clusters)[15] <- "IFI44L"
levels(immune.combined$seurat_clusters)[16] <- "MALAT1"
cluster <- immune.combined$seurat_clusters
cl <- c("FGFBP2", "DUSP2", "RPLs", "TNFRSF4", "CCR7", "ISG15", "NKG7", "MAP3K8",  
        "GZMK", "PRF1", "F5", "CD79A", "LYZ", "PPBP", "IFI44L", "MALAT1")
IDX <- list()
for (i in 1:length(cl)){
  IDX[[i]] <- which(cluster==cl[i])
}
#gene select
geneList_selected <- c("FGFBP2", "GZMH", "GZMB", "DUSP2", 
                       "CD74", "S100A4", "HOPX", 
                      "TNFRSF4", "CORO1B", "CCR7", 
                      "SELL", "ISG15", "IFI6", "XAF1", "EIF2AK2", 
                      "CCL5", "KLRD1", "CST7", "CD160", 
                      "DUSP1", "GZMK", "JUN", "PMAIP1", 
                      "KLRB1", "CYB561", "S100A11", 
                      "LYAR", "GIMAP4", "NEAT1", "PRF1", 
                      "RTKN2", "IKZF2", "TNFRSF18", 
                      "DUSP4", "TIGIT", "MS4A1", "CD79A", "HLA-DRA", "CD79B", 
                      "FTL", "CTSS", "CST3", 
                      "PPBP", "TAGLN2", "FTH1", "PTGS1", "PF4", "CTSA", 
                      "GNG11", "MGLL", "RAB32", 
                      "EEF1B2", "IFITM3", 
                      "HBA2", "PARP9", "IFI44L")
idx <- matrix(1, 1, length(geneList_selected))
for (i in 1:length(geneList_selected)){
  idx[i] <- which(rownames(immune.combined)  %in%  geneList_selected[i])
}
#express matrix
expression <- matrix(0, length(geneList_selected), length(cl))
for (i in 1:length(cl)){
  buffer <- as.matrix(immune.combined@assays$RNA@data[,IDX[[i]]])
  buffer <- buffer[idx,]
  expression[,i] <- rowMeans(buffer, na.rm = TRUE)
}
rownames(expression) <- geneList_selected
colnames(expression) <- c("FGFBP2", "DUSP2", "RPLs", "TNFRSF4", "CCR7", "ISG15", "NKG7", "MAP3K8",  
                          "GZMK", "PRF1", "F5", "CD79A", "LYZ", "PPBP", "IFI44L", "MALAT1")
# HeatMap
pheatmap(expression, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, 
         show_colnames = TRUE, show_rownames = TRUE, gaps_row = c(4,16,20, 25, 33, 36), gaps_col = c(4, 6),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))






########################### ViolinPlot for samples
# sample
name <- c("Health", "Patient", "Treatment")
IDX <- list()
for (i in 1:length(name)){
  IDX[[i]] <- which(immune.combined$stim==name[i])
}
#gene select
geneList_selected2 <- c("CD8B", "CD4", "IL7R", "FOXP3", "IL2RA", "TNFRSF18", 
                        "CTLA4", "TIGIT", "CCR7", "TNFRSF4", "ICOS", "CD27", 
                        "BATF", "ISG15", "ISG20", "IFITM3", "IFI44L", "IFI6")
idx <- matrix(1, 1, length(geneList_selected2))
for (i in 1:length(geneList_selected2)){
  idx[i] <- which(rownames(immune.combined)  %in%  geneList_selected2[i])
}
#express array
sample <- c("health", "patient", "treatment")
gene <- c("CD8B", "CD4", "IL7R", "FOXP3", "IL2RA", "TNFRSF18", 
          "CTLA4", "TIGIT", "CCR7", "TNFRSF4", "ICOS", "CD27", 
          "BATF", "ISG15", "ISG20", "IFITM3", "IFI44L", "IFI6")
cell <- (1:12160)
expression <- array(data = NA, dim = c(12160, 3, length(geneList_selected2)), dimnames = list(cell, sample, gene))
for (i in 1:length(name)){
  buffer <- as.matrix(immune.combined@assays$RNA@data[,IDX[[i]]])
  buffer <- buffer[idx,]
  for(j in 1:ncol(buffer)){
    expression[j,i,] <- buffer[,j]
  }
}
# from array to data.frame
funforplot<-function(my) {
  mydata<-as.data.frame(my)
  mydata$cell<-rownames(my)
  mydata<-gather(mydata, key=sample, value = value, -cell)
  mydata
}
exp.frame <- adply(expression,3,.fun=funforplot)
exp.frame <- transform(exp.frame,time=as.character(gene))
exp.frame <- na.omit(exp.frame)
exp.frame <- exp.frame[,1:4]
colnames(exp.frame) = c("gene", "cell", "sample", "value")
# ViolinPlot
plot_list <- list()
data_list <- list()
title_list = c("CD8B", "CD4", "IL7R", "FOXP3", "IL2RA", "TNFRSF18", 
               "CTLA4", "TIGIT", "CCR7", "TNFRSF4", "ICOS", "CD27", 
               "BATF", "ISG15", "ISG20", "IFITM3", "IFI44L", "IFI6")
for(i in 1:18){
  data_list[[i]] <- exp.frame[(30890*i-30889):(30890*i),]
  data_list[[i]][1:8962,5] <- mean(data_list[[i]][1:8962,4])
  data_list[[i]][8963:21122,5] <- mean(data_list[[i]][8963:21122,4])
  data_list[[i]][21123:30890,5] <- mean(data_list[[i]][21123:30890,4])
  plot_list[[i]] <- ggplot(data_list[[i]]) + geom_violin(aes(x=sample, y=value, fill=V5)) + 
                            scale_fill_gradient(low="yellow", high="red")+
                            labs(title=title_list[i], x=NULL, y = NULL) +
                            theme(legend.position="none") +
                            theme(plot.title = element_text(hjust = 0.5))
}
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],plot_list[[5]],
             plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],plot_list[[10]],
             plot_list[[11]],plot_list[[12]],plot_list[[13]],plot_list[[14]],plot_list[[15]],
             plot_list[[16]],plot_list[[17]],plot_list[[18]], nrow = 6)






########################### UMAP by samples
# read data

immune.combined.markers <- FindAllMarkers(object = immune.combined, assay = "integrated", only.pos = TRUE, min.pct = 0.1, thresh.use = 0.25, test.use = "roc")
#write.table(immune.combined.markers , "Combined_ThreeMarkerGenesByCluster_immuneCombi_default.xls",sep="\t")
#top5 <- immune.combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_diff)
#DoHeatmap(immune.combined, features = top5$gene, group.by = 'ident') #+ NoLegend()
# Gene expression plot
FeaturePlot(immune.combined, features = c("PRF1", "GZMA", "GZMH", "FGFBP2", "ISG15", "ISG20", "IFITM3", "IFI44L", "IFI6"), split.by = "stim")

#FeaturePlot(immune.combined, features = c("PRF1"), split.by = "stim")






########################### VolcanoPlot for DEgenes
# read data
immune.combined <- readRDS(file = "H:/???ཿ ??ɽҽԺ/1/immune.combined.three.rds")
sample <- immune.combined$stim
cluster <- immune.combined$seurat_clusters
cluster <- matrix(cluster, ncol=1)
cluster <- as.numeric(cluster)
idx7 <- c(which(cluster==7))
immune.combined.7 <- immune.combined[,idx7]
DefaultAssay(immune.combined.7) <- "RNA"
idx7_h <- which(immune.combined.7$stim=="Health")
idx7_p <- which(immune.combined.7$stim=="Patient")
idx7_t <- which(immune.combined.7$stim=="Treatment")

# patient vs health
buffer <- cbind(as.matrix(immune.combined.7@assays$RNA@data[,idx7_h]), as.matrix(immune.combined.7@assays$RNA@data[,idx7_p]))
category <- c(rep("health", length(idx7_h)), rep("patient", length(idx7_p)))
cData <- category
cData <- data.frame(cData)
colnames(buffer) <- 1:dim(buffer)[2]
fData <- rownames(buffer)
fData <- data.frame(fData)
scaRaw <- as.matrix(buffer)
sca <- FromMatrix(scaRaw, cData, fData)
cond<-factor(colData(sca)$cData)
cond<-relevel(cond,"patient")
colData(sca)$cData<-cond
zlmCond <- zlm(~cData, sca)
summaryCond <- summary(zlmCond, doLRT='cDatahealth') 
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='cDatahealth' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='cDatahealth' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
p <- fcHurdle$`Pr(>Chisq)`
p.adj <- p.adjust(p, method = "fdr", n = length(p))
fcHurdle$`Pr(>Chisq)` <- p.adj
#write.table(fcHurdle, file = "Myloid patient vs health.csv", row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
idx <- which(!is.na(fcHurdle$coef))
fcHurdle <- fcHurdle[idx,]
idx <- which(!is.na(fcHurdle$`Pr(>Chisq)`))
fcHurdle <- fcHurdle[idx,]
idx_DE<- intersect(which(fcHurdle$`Pr(>Chisq)`<0.05), which(abs(fcHurdle$coef)>1))
logp <- -log10(fcHurdle$`Pr(>Chisq)`)
log2FC <- fcHurdle$coef
DE <- rep("No significant", length(log2FC))
DE[idx_DE] <- "Significant"
X <- data.frame(log10P=logp, log2FC= log2FC, DE = DE, geneID = fcHurdle$primerid)
rownames(X) <- fcHurdle$primerid
# VolcanoPlot
ggplot() + geom_point(data=X, aes(x=log2FC, y=log10P, color=DE))+ theme_classic() +
  scale_color_manual(values=c('grey', 'red'))+
  geom_text_repel(
    data = X[idx_DE,],
    aes(x=log2FC, y=log10P, label = geneID),size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))

#treatment vs health
buffer <- cbind(as.matrix(immune.combined.7@assays$RNA@data[,idx7_h]), as.matrix(immune.combined.7@assays$RNA@data[,idx7_t]))
category <- c(rep("health", length(idx7_h)), rep("treatment", length(idx7_t)))
cData <- category
cData <- data.frame(cData)
colnames(buffer) <- 1:dim(buffer)[2]
fData <- rownames(buffer)
fData <- data.frame(fData)
scaRaw <- as.matrix(buffer)
sca <- FromMatrix(scaRaw, cData, fData)
cond<-factor(colData(sca)$cData)
cond<-relevel(cond,"treatment")
colData(sca)$cData<-cond
zlmCond <- zlm(~cData, sca)
summaryCond <- summary(zlmCond, doLRT='cDatahealth') 
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='cDatahealth' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='cDatahealth' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
p <- fcHurdle$`Pr(>Chisq)`
p.adj <- p.adjust(p, method = "fdr", n = length(p))
fcHurdle$`Pr(>Chisq)` <- p.adj
#write.table(fcHurdle, file = "Myloid treatment vs health.csv", row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
idx <- which(!is.na(fcHurdle$coef))
fcHurdle <- fcHurdle[idx,]
idx <- which(!is.na(fcHurdle$`Pr(>Chisq)`))
fcHurdle <- fcHurdle[idx,]
idx_DE<- intersect(which(fcHurdle$`Pr(>Chisq)`<0.05), which(abs(fcHurdle$coef)>1))
logp <- -log10(fcHurdle$`Pr(>Chisq)`)
log2FC <- fcHurdle$coef
DE <- rep("No significant", length(log2FC))
DE[idx_DE] <- "Significant"
X <- data.frame(log10P=logp, log2FC= log2FC, DE = DE, geneID = fcHurdle$primerid)
rownames(X) <- fcHurdle$primerid
# VolcanoPlot
ggplot() + geom_point(data=X, aes(x=log2FC, y=log10P, color=DE))+
  scale_color_manual(values=c('grey', 'red'))+
  geom_text_repel(
    data = X[idx_DE,],
    aes(x=log2FC, y=log10P, label = geneID),size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))

#treatment vs patient
buffer <- cbind(as.matrix(immune.combined.7@assays$RNA@data[,idx7_p]), as.matrix(immune.combined.7@assays$RNA@data[,idx7_t]))
category <- c(rep("patient", length(idx7_p)), rep("treatment", length(idx7_t)))
cData <- category
cData <- data.frame(cData)
colnames(buffer) <- 1:dim(buffer)[2]
fData <- rownames(buffer)
fData <- data.frame(fData)
scaRaw <- as.matrix(buffer)
sca <- FromMatrix(scaRaw, cData, fData)
cond<-factor(colData(sca)$cData)
cond<-relevel(cond,"treatment")
colData(sca)$cData<-cond
zlmCond <- zlm(~cData, sca)
summaryCond <- summary(zlmCond, doLRT='cDatapatient') 
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='cDatapatient' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='cDatapatient' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
p <- fcHurdle$`Pr(>Chisq)`
p.adj <- p.adjust(p, method = "fdr", n = length(p))
fcHurdle$`Pr(>Chisq)` <- p.adj
#write.table(fcHurdle, file = "Myloid treatment vs patient.csv", row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
idx <- which(!is.na(fcHurdle$coef))
fcHurdle <- fcHurdle[idx,]
idx <- which(!is.na(fcHurdle$`Pr(>Chisq)`))
fcHurdle <- fcHurdle[idx,]
idx_DE<- intersect(which(fcHurdle$`Pr(>Chisq)`<0.05), which(abs(fcHurdle$coef)>1))
logp <- -log10(fcHurdle$`Pr(>Chisq)`)
log2FC <- fcHurdle$coef
DE <- rep("No significant", length(log2FC))
DE[idx_DE] <- "Significant"
X <- data.frame(log10P=logp, log2FC= log2FC, DE = DE, geneID = fcHurdle$primerid)
rownames(X) <- fcHurdle$primerid
# VolcanoPlot
ggplot() + geom_point(data=X, aes(x=log2FC, y=log10P, color=DE))+
  scale_color_manual(values=c('grey', 'red'))+
  geom_text_repel(
    data = X[idx_DE,],
    aes(x=log2FC, y=log10P, label = geneID),size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"))






########################### monocle
immune.combined <- readRDS(file = "H:/???ཿ ??ɽҽԺ/1/immune.combined.three.rds")


# CD8
cluster <- immune.combined$seurat_clusters
idx_CD8 <- c(which(cluster==0), which(cluster==4), which(cluster==6), which(cluster==7),
             which(cluster==8), which(cluster==9))
immune.combined.CD8 <- immune.combined[,idx_CD8]
DefaultAssay(immune.combined.CD8) <- "RNA"
sample <- immune.combined.CD8$stim
idx_sample <- c(which(sample=="Health"),which(sample=="Patient"),which(sample=="Treatment"))
# monocle workflow
buffer <- as.matrix(immune.combined.CD8@assays$RNA@data[,idx_sample])
buffer <- 2^buffer-1
SeuratP <- CreateSeuratObject(counts = buffer, min.cells = 3)
SeuratP <- subset(SeuratP, subset = nFeature_RNA > 700)
SeuratP <- NormalizeData(SeuratP, verbose = FALSE)
SeuratP <- FindVariableFeatures(SeuratP, selection.method = "vst", nfeatures = 2000)
geneList <- rownames(immune.combined.CD8@assays$RNA@data)
geneList_selected <- SeuratP@assays$RNA@var.features
idx_selected <- match(geneList_selected, geneList)
idx_selected <- na.omit(idx_selected)

# by sample
x <- 1:length(idx_sample)
idx <- sample(x, size = 3000, replace = FALSE)
idx_sample <- idx_sample[idx]
sample <- as.character(sample[idx_sample])
expression <- as.matrix(immune.combined.CD8@assays$RNA@data[idx_selected,idx_sample])
gene_feature <- rownames(expression)
gene_feature <- as.matrix(gene_feature)
rownames(gene_feature) <- gene_feature
colnames(gene_feature) <- "gene_short_name"
names(sample) <- colnames(expression)
names(cluster) <- colnames(expression)
gene_feature <- data.frame(gene_feature)
sample <- data.frame(sample)
cluster <- data.frame(cluster)
pd <- new("AnnotatedDataFrame", data = sample)
fd <- new("AnnotatedDataFrame", data = gene_feature)
cds <- newCellDataSet(expression, phenoData = pd, featureData = fd)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, ordering_genes$gene_id)
#cds_reduced <- reduceDimension(cds);
cds_reduced <- reduceDimension(cds, max_components = 2, num_dim = 6, reduction_method = 'DDRTree', verbose = T)
#cds_reduced <- clusterCells(cds_reduced);
cds_for_pseudo <- orderCells(cds_reduced, reverse = TRUE)
#cds_for_pseudo <- orderCells(cds_for_pseudo, root_state = which(cds_for_pseudo$cluster==1))
#write.table(pData(cds_for_pseudo),"NewPhenoData.txt",sep="\t",quote=F,row.names = T)
#saveRDS(cds_for_pseudo, file = "NKTcell_CD8_traj_randIdx6K_cds_for_pseudo.rds")
#saveRDS(idx_CD8, file = "NKTcell_CD8_randIdx6K_idx_CD8.rds")
#cds_for_pseudo <- readRDS(file = "NKTcell_CD8_traj.rds") # 10000 random picked cells
#idx_CD8 <- readRDS(file = "NKTcell_CD8_randIdx10K.rds") # 10000 random picked cells
plot_cell_trajectory(cds_for_pseudo, markers_linear = TRUE, color_by = "sample") +  ggtitle("CD8+")+ scale_color_brewer(palette="Set1")
plot_cell_trajectory(cds_for_pseudo, color_by = "sample", show_branch_points = T, cell_size = 2, cell_link_size = 0.3) + 
  facet_wrap(~sample, nrow = 1) +theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.position="right", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  theme (legend.position="top", legend.title=element_blank())

# by cluster
cds_for_pseudo$cluster <- immune.combined.CD8$seurat_clusters[match(colnames(cds_for_pseudo), colnames(immune.combined.CD8))]
plot_cell_trajectory(cds_for_pseudo, markers_linear = TRUE, color_by = "cluster") +  ggtitle("CD8+")+ scale_color_brewer(palette="Set1")
plot_cell_trajectory(cds_for_pseudo, color_by = "cluster", show_branch_points = T, cell_size = 2, cell_link_size = 0.3) + 
  facet_wrap(~cluster, nrow = 1) +theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.position="right", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  theme (legend.position="top", legend.title=element_blank())






# CD4
cluster <- immune.combined$seurat_clusters
idx_CD4 <- c(which(cluster==2), which(cluster==3), which(cluster==5), which(cluster==11),
             which(cluster==12))
immune.combined.CD4 <- immune.combined[,idx_CD4]
DefaultAssay(immune.combined.CD4) <- "RNA"
sample <- immune.combined.CD4$stim
idx_sample <- c(which(sample=="Health"),which(sample=="Patient"),which(sample=="Treatment"))

# monocle workflow
buffer <- as.matrix(immune.combined.CD4@assays$RNA@data[,idx_sample])
buffer <- 2^buffer-1
SeuratP <- CreateSeuratObject(counts = buffer, min.cells = 3)
SeuratP <- subset(SeuratP, subset = nFeature_RNA > 700)
SeuratP <- NormalizeData(SeuratP, verbose = FALSE)
SeuratP <- FindVariableFeatures(SeuratP, selection.method = "vst", nfeatures = 2000)
geneList <- rownames(immune.combined.CD4@assays$RNA@data)
geneList_selected <- SeuratP@assays$RNA@var.features
idx_selected <- match(geneList_selected, geneList)
idx_selected <- na.omit(idx_selected)

# by sample
x <- 1:length(idx_sample)
idx <- sample(x, size = 3000, replace = FALSE)
idx_sample <- idx_sample[idx]
sample <- as.character(sample[idx_sample])
expression <- as.matrix(immune.combined.CD4@assays$RNA@data[idx_selected,idx_sample])
gene_feature <- rownames(expression)
gene_feature <- as.matrix(gene_feature)
rownames(gene_feature) <- gene_feature
colnames(gene_feature) <- "gene_short_name"
names(sample) <- colnames(expression)
gene_feature <- data.frame(gene_feature)
sample <- data.frame(sample)
pd <- new("AnnotatedDataFrame", data = sample)
fd <- new("AnnotatedDataFrame", data = gene_feature)
cds <- newCellDataSet(expression, phenoData = pd, featureData = fd)
cds<- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds)
ordering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, ordering_genes$gene_id)
#cds_reduced <- reduceDimension(cds);
cds_reduced <- reduceDimension(cds, max_components = 2, num_dim = 6, reduction_method = 'DDRTree', verbose = T)
#cds_reduced <- clusterCells(cds_reduced);
cds_for_pseudo <- orderCells(cds_reduced, reverse = TRUE)
#cds_for_pseudo <- orderCells(cds_for_pseudo, root_state = which(cds_for_pseudo$cluster==1))
#write.table(pData(cds_for_pseudo),"NewPhenoData.txt",sep="\t",quote=F,row.names = T)
#saveRDS(cds_for_pseudo, file = "NKTcell_CD8_traj_randIdx6K_cds_for_pseudo.rds")
#saveRDS(idx_CD8, file = "NKTcell_CD8_randIdx6K_idx_CD8.rds")
#cds_for_pseudo <- readRDS(file = "NKTcell_CD8_traj.rds") # 10000 random picked cells
#idx_CD8 <- readRDS(file = "NKTcell_CD8_randIdx10K.rds") # 10000 random picked cells
plot_cell_trajectory(cds_for_pseudo, markers_linear = TRUE, color_by = "sample") +  ggtitle("CD4+")+ scale_color_brewer(palette="Set1")
#+scale_color_gradient2(midpoint=mid, low="blue", mid="white", high="red", space ="Lab" )
#plot_cell_trajectory(cds_for_pseudo , markers_linear = TRUE, color_by = "Pseudotime")  +  ggtitle("CD8+") + scale_color_gradient(low="blue", high="red")
plot_cell_trajectory(cds_for_pseudo, color_by = "sample", show_branch_points = T, cell_size = 2, cell_link_size = 0.3) + 
  facet_wrap(~sample, nrow = 1) +theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.position="right", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  theme (legend.position="top", legend.title=element_blank())

# by cluster
cds_for_pseudo$cluster <- immune.combined.CD4$seurat_clusters[match(colnames(cds_for_pseudo), colnames(immune.combined.CD4))]
plot_cell_trajectory(cds_for_pseudo, markers_linear = TRUE, color_by = "cluster") +  ggtitle("CD4+")+ scale_color_brewer(palette="Set1")
plot_cell_trajectory(cds_for_pseudo, color_by = "cluster", show_branch_points = T, cell_size = 2, cell_link_size = 0.3) + 
  facet_wrap(~cluster, nrow = 1) +theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(legend.position="right", legend.title=element_blank()) +
  stat_density2d(color='black', h = 8, alpha=I(0.25), size=I(0.25)) +
  theme (legend.position="top", legend.title=element_blank())
