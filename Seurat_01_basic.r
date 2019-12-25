#Seurat3.1.1 基本分析

########################### 1.加载数据集(pbmc外周血)
#install.packages("Seurat")
#install.packages("dplyr")
#install.packages("ggsci")
library(Seurat)
#packageVersion("Seurat")
library(dplyr)
library(ggsci)
list.files("/Users/mkbai/fsdownload/filtered_gene_bc_matrices/hg19")
??Reads10X
pbmc.data <- Read10X(data.dir = "/Users/mkbai/fsdownload/filtered_gene_bc_matrices/hg19")
pbmc.data

########################### 2.用raw data初始化Seurat变量
#The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc#pbmc
#An object of class Seurat 
#15651 features across 9383 samples within 1 assay 
#Active assay: RNA (15651 features)

########################### 3.质控QC和过滤细胞-clear data
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
head(pbmc@meta.data, 10)
#过滤金标准
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by
# the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
#筛选
pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 5)

########################### 4.预处理
#标准化
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#确定变异系数大的基因
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
#Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = pbmc), 10)
#plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
#数据比例缩放
#all.genes <- rownames(pbmc)
#pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc)

########################### 5.PCA降维
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
#确定聚类需要的PCA维度
#5.1 JackStrawPlot函数
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
#5.2 Elbow Plot法
ElbowPlot(pbmc)

########################### 6.细胞聚类
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

########################### 7.非线性降维
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

########################### 8.寻找差异表达基因
#寻找clusterX中的差异表达基因
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, ident.2 = c(1, 3), min.pct = 0.25)
head(cluster2.markers, n = 10)
#寻找每个cluster与其他cluster相比后的差异基因marker
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(pbmc, features = c("LYZ", "CST3", "CD74"))
##SCTransform()函数可代替2.0版本的NormalizeData, ScaleData, FindVariableFeature
##pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
