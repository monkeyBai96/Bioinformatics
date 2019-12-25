#Seurat3.1.1 多样本整合与统计

########################### 1.数据准备与标准流程
#setwd("/Users/mkbai/Documents/Seurat/")
#rm(list=ls())
#graphics.off()
#library(Seurat)
#library(cowplot)
#library(dplyr)
#library(R.matlab)
#测序assay不同 其路径名也有差别
health.data <-  Read10X(data.dir = "/Users/mkbai/fsdownload/sample1/filtered_gene_bc_matrices/hg19/")
patient.data <- Read10X(data.dir = "/Users/mkbai/fsdownload/sample2/filtered_feature_bc_matrix/")
treatment.data <- Read10X(data.dir = "/Users/mkbai/fsdownload/sample3/filtered_feature_bc_matrix/")
# health样本标准化
health <- CreateSeuratObject(counts = health.data, project = "Health", min.cells = 3, min.features = 200)
health$stim <- "Health"
health[["percent.mt"]] <- PercentageFeatureSet(object = health, pattern = "^MT-")
health <- subset(health, subset = nFeature_RNA > 500 & nFeature_RNA < 2100 & percent.mt < 5)
health <- NormalizeData(health, verbose = FALSE)
health <- FindVariableFeatures(health, selection.method = "vst", nfeatures = 2000)
# patient样本标准化
patient <- CreateSeuratObject(counts = patient.data, project = "Patient", min.cells = 3, min.features = 200)
patient$stim <- "Patient"
patient[["percent.mt"]] <- PercentageFeatureSet(object = patient, pattern = "^MT-")
patient <- subset(patient, subset = nFeature_RNA > 800 & nFeature_RNA < 2500 & percent.mt < 10)
patient <- NormalizeData(patient, verbose = FALSE)
patient <- FindVariableFeatures(patient, selection.method = "vst", nfeatures = 2000)
# treatment样本标准化
treatment <- CreateSeuratObject(counts = treatment.data, project = "Treatment", min.cells = 3, min.features = 200)
treatment$stim <- "Treatment"
treatment[["percent.mt"]] <- PercentageFeatureSet(object = treatment, pattern = "^MT-")
treatment <- subset(treatment, subset = nFeature_RNA > 800 & nFeature_RNA < 2100 & percent.mt < 10)
treatment <- NormalizeData(treatment, verbose = FALSE)
treatment <- FindVariableFeatures(treatment, selection.method = "vst", nfeatures = 2000)

########################### combined samples
immune.anchors <- FindIntegrationAnchors(object.list = list(health, patient, treatment), dims = 1:20, anchor.features = 2000)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
#temp <- intersect(geneList_sOGFSC, immune.anchors@anchor.features)
#immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20, features.to.integrate = temp)
DefaultAssay(immune.combined) <- "integrated"

########################### standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
# UMAP/t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 1.0)
#保存
saveRDS(immune.combined, file = "immune.combined.three.rds")
#读取
immune.combined <- readRDS(file = "immune.combined.three.rds")




########################### 2.整合数据的基本统计
cluster <- immune.combined$seurat_clusters
ident <- matrix("0", dim(immune.combined)[2],1)
cellTypeLabel <- ident
# cluster identify
levels(immune.combined@active.ident)[3] <- "RPLs"
levels(immune.combined@active.ident)[1] <- "FGFBP2"
levels(immune.combined@active.ident)[4] <- "TNFRSF4"
levels(immune.combined@active.ident)[9] <- "NKG7"
levels(immune.combined@active.ident)[12] <- "RPLs"
levels(immune.combined@active.ident)[10] <- "NKG7"
levels(immune.combined@active.ident)[7] <- "FGFBP2"
levels(immune.combined@active.ident)[6] <- "TNFRSF4"
levels(immune.combined@active.ident)[2] <- "DUSP2"
levels(immune.combined@active.ident)[5] <- "CCR7"
levels(immune.combined@active.ident)[6] <- "ISG15"
levels(immune.combined@active.ident)[8] <- "MAP3K8"
levels(immune.combined@active.ident)[9] <- "GZMK"
levels(immune.combined@active.ident)[10] <- "PRF1"
levels(immune.combined@active.ident)[11] <- "F5"
levels(immune.combined@active.ident)[12] <- "CD79A"
levels(immune.combined@active.ident)[13] <- "LYZ"
levels(immune.combined@active.ident)[14] <- "PPBP"
levels(immune.combined@active.ident)[15] <- "IFI44L"
levels(immune.combined@active.ident)[16] <- "MALAT1"
for(i in 1:length(ident)) {
  immune.combined@active.ident[i] <- ident[i]
}
# Visualization
# 按照样本名着色
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
# 按照cluster着色
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
write.table(p1$data, "CombinedSeveralCluster.xls", sep="\t", row.names = TRUE, col.names = TRUE)
# find markers
immune.combined.markers <- FindAllMarkers(object = immune.combined, assay = "integrated", only.pos = TRUE, min.pct = 0.1, thresh.use = 0.25, test.use = "roc")
#immune.combined.markers <- FindAllMarkers(object = immune.combined, only.pos = TRUE, min.pct = 0.1, thresh.use = 0.25, test.use = "roc")
#write.table(immune.combined.markers , "MarkerGenesByCluster_immuneCombi_sOGFSC.xls",sep="\t")
write.table(immune.combined.markers , "Combined_ThreeMarkerGenesByCluster_immuneCombi_default.xls",sep="\t")
top5 <- immune.combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_diff)
DoHeatmap(immune.combined, features = top5$gene, group.by = 'ident') #+ NoLegend()

########################### 比较各个样本中关键基因
# health要经过NormalizeData，FindVariableFeatures，ScaleData，
# RunPCA/RunUMAP，FindNeighbors，FindClusters等步骤
FeaturePlot(health, features = c("PRF1", "GZMA", "GZMH", "FGFBP2", "ISG15", "ISG20", "IFITM3", "IFI44L", "IFI6"))
FeaturePlot(patient, features = c("PRF1", "GZMA", "GZMH", "FGFBP2", "ISG15", "ISG20", "IFITM3", "IFI44L", "IFI6"))
FeaturePlot(treatment, features = c("PRF1", "GZMA", "GZMH", "FGFBP2", "ISG15", "ISG20", "IFITM3", "IFI44L", "IFI6"))


########################### 统计cluster中各个sample的比例
x <- matrix(0, 3, 19)
for (i in 1:19){
  idx <- which(immune.combined$seurat_clusters == i-1)
  x[1,i] <- length(which(immune.combined$stim[idx]== "Health"))
  x[2,i] <- length(which(immune.combined$stim[idx]== "Patient"))
  x[3,i] <- length(which(immune.combined$stim[idx]== "Treatment"))
}
write.table(x, file = "/Users/mkbai/Documents/Seurat/statistics.csv")
#install.packages("RColorBrewer")
library(ggplot2)
library(RColorBrewer)
cluster <- read.csv("/Users/mkbai/Documents/Seurat/cluster_frame.csv", header = TRUE)
data <- data.frame(cluster)

ggplot(data=data, aes(x=cluster, y=cell.per, fill=sample)) +  geom_bar(stat="identity") + 
  scale_colour_manual(values=brewer.pal(12, "Paired")) + labs(x="cluster")

