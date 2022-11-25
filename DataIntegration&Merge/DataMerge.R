######***********************######
##First time use,install packages##
######***********************######
###################################
# source("/home/scripts/R_scripts/newest scripts/functions/source.R") ##调用函数与包
# .libPaths("/home/R/x86_64-redhat-linux-gnu-library/3.6")
# .libPaths("/home/R/x86_64-pc-linux-gnu-library/4.0")

######***********************######
######**Install R packages **######
######***********************######
install.packages("hdf5r")

######***********************######
######**library R packages **######
######***********************######
library(ggplot2)
library(stringr)
library(Seurat)
library(magrittr)
library(dplyr)
library(hdf5r)

###################################
setwd("/home/singlecell/rawdata")
DefaultAssay(a.object)
DefaultAssay(b.object)

SPT_A_20210627$samples=rep("S1",ncol(SPT_A_20210627))
SPTB$samples=rep("S2",ncol(SPTB))
SPTC$samples=rep("S3",ncol(SPTC))
SPTD$samples=rep("S4",ncol(SPTD))

SPTABCD= merge(SPT_A_20210627,SPTB)
SPTABCD= merge(SPTABCD,SPTC)
SPTABCD= merge(SPTABCD,SPTD)
?merge
##作图，初步看下数据分布
VlnPlot(SPTABCD, features = "nCount_Spatial",pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(SPTABCD, features = "nCount_Spatial") + theme(legend.position = "right")

########################################################################
##########################引用单细胞测序降维方式########################
#########################数据预处理
###推荐sctransform
SPTABCD <- SCTransform(SPTABCD, assay = "Spatial", verbose = FALSE)

###UMAP降维
SPTABCD <- RunPCA(SPTABCD, assay = "SCT", verbose = FALSE) %>% #运行PCA运算
  FindNeighbors(reduction = "pca", dims = 1:30) %>% #选择dims
  FindClusters(verbose = FALSE) %>% #选择resolution
  RunUMAP(reduction = "pca", dims = 1:30) #tsne运算
#看降维效果
DimPlot(SPTABCD, reduction = "umap", label = TRUE, cols = colors_aa)
SpatialDimPlot(SPTABCD, label = FALSE, label.size = 1, cols = colors_aa)
#寻找高变基因
SPTD <- FindSpatiallyVariableFeatures(SPTD, assay = "SCT", features = VariableFeatures(SPTD)[1:1000], 
                                      selection.method = "markvariogram")
SPTD.top.features <- head(SpatiallyVariableFeatures(SPTD, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(SPTD, features = SPTD.top.features, ncol = 3, alpha = c(0.1, 1),crop = FALSE)