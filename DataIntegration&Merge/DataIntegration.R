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

umi.list <- list(a.object, b.object)
umi <- FindIntegrationAnchors(object.list = umi.list, 
                              dims = 1:40,reduction = "cca")
# umi <- FindIntegrationAnchors(object.list = umi.list, 
#                               dims = 1:40,reduction= "rpca")
umi.integrated <- IntegrateData(anchorset = umi, dims = 1:30)
SeuratObject::DefaultAssay(umi.integrated) <- "integrated"

umi.integrated <- ScaleData(umi.integrated, verbose = T)
umi.integrated <- RunPCA(umi.integrated, verbose = T)
DimPlot(object = umi.integrated, reduction = "pca")
ElbowPlot(object = umi.integrated, ndims = 50)
umi.integrated <- FindNeighbors(object = umi.integrated, 
                                dims = 1:5)
umi.integrated <- FindClusters(object = umi.integrated, 
                               resolution = 1)
umi.integrated <- RunTSNE(object = umi.integrated, dims = 1:5, 
                          check_duplicates = FALSE)
umi.integrated <- RunUMAP(umi.integrated, reduction = "pca", 
                          dims = 1:5)
integrated.object = umi.integrated
integrated.object <- FindClusters(object = integrated.object,  ###调整分辨率，重新作图
                             resolution = 0.7)
integrated.object <- RunUMAP(integrated.object, reduction = "pca", 
                        dims = 1:5)
integrated.object <- RunTSNE(integrated.object, reduction = "pca", 
                        dims = 1:5)

###绘图
levels(integrated.object@active.ident)
DimPlot(object = integrated.object, reduction = "umap", split.by = "samples",cols= colors)
DimPlot(object = integrated.object, reduction = "tsne", group.by = "samples",cols= colors)
SpatialDimPlot(integrated.object,label = FALSE,label.size = 1,crop = FALSE,pt.size.factor = 1,cols = colors_figure)

#导出pdf
FeaturePlot(integrated.object, features = 'nFeature_SCT')
FeaturePlot(integrated.object, features = 'nCount_SCT')
pdf(paste0("20221109_integrated.object_5PCA_DimPlot2.pdf"),width = 9,height = 8)
p=DimPlot(integrated.object,reduction = 'umap',cols = color2)
print(p)
dev.off()

levels(integrated.object)
color2 <- colors_figure
names(color2)=levels(integrated.object)
pdf(paste0("20221115_integrated.object_5PCA_DimPlot2.pdf"),width = 30,height = 8)
p=SpatialDimPlot(integrated.object,label = FALSE,label.size = 1,crop = FALSE,pt.size.factor = 1.35,cols = color2)
print(p)
dev.off()

#计算比例
Get_percentplot_split(integrated.object,filename = '20221109_integrated.object_5PCA')
#导出pdf
levels(integrated.object)
names(color2)=levels(integrated.object)
pdf(paste0("20221110_integrated.object_5PCA_SpatialDimPlot.pdf"),width = 11.69,height = 8.27)
p=SpatialDimPlot(integrated.object,label = TRUE,label.size = 1,crop = FALSE,pt.size.factor = 1,cols = color2)
print(p)
dev.off()
###由于SpatialDimPlot本身不适合换颜色，将颜色加标签，后改进
levels(integrated.object)
names(color2)=levels(integrated.object)
SpatialDimPlot(integrated.object,label = TRUE,label.size = 1,crop = FALSE,pt.size.factor = 1.4,cols = colors_figure)
DimPlot(integrated.object,reduction = 'umap',cols = colors_aa)
