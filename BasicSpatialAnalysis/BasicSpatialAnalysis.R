######***********************######
##First time use,install packages##
######***********************######
###################################
# source("/home/longjie/scripts/R_scripts/newest scripts/functions/source.R") ##调用函数与包
# .libPaths("/home/yankai/R/x86_64-redhat-linux-gnu-library/3.6")
# .libPaths("/home/R/x86_64-pc-linux-gnu-library/4.0")
### 调色方案
colors_aa=c("#F2D778","#B785BF","#BD8D8F","#B8DDFD","#FF6A2A","#FFC8BC","#FFB700","#FF6B8F","#24504F","#6DA1F2",
            "#398C95","#68C589","#F8BEEF","#FFF7A4","#F262BA","#AF19D2","#7B4A42","#D886B3","#F7BBD7","#F88ED0",
            "#73C3DA","#DFDBF0","#FFE3E0","#6F56D4","#F6B659","#3B95C8","#E6B7D7","#5E75B7",
            "#FF0090","#918DC0","#00C4FF","#FED293","#D7EBF8","#C4C4B4","#FF8985","#F2998A",
            "#FFB9C9","#67D9FB","#D8B385","#F2DFE7","#6A88D9","#C0DAE0","#977195","#D566AB","#BD86BF","#FFE3AF",
            "#EDD5ED","#DAB490","#E8696B")
colors_spt =c ( "#D7EBF8",  "#73C3DA",  "#67D9FB",  "#FFB9C9",  "#39B5F9",  "#F7F5AF",  "#6A88D9",  "#00C4FF",  "#24504F",  "#FFB700", 
                "#FF5EA9",  "#FFE3E0",  "#6F56D4",  "#FFE3AF",  "#D8B385",  "#C4C4B4",  "#ACA7C9",  "#BD8D8F",  "#398C95",  "#B8DDFD", 
                "#7B4A42",  "#FF0090",  "3A43245",  "#F8BEEF",  "#BD86BF",  "#F5B680",  "#FF6A2A",  "#AF19D2",  "#D886B3",  "#EDD5ED", 
                "#6DA1F2",  "#9AC6A2",  "#F2D778",  "#DAB490",  "#F2DFE7",  "#F7BBD7",  "#FED293",  "#DA705A",  "#F6B659",  "#FBE2DA",
                "#FF8985",  "#FFC8BC",  "#D566AB",  "#FFF7A4",  "#5E75B7",  "#F2998A",  "#E8696B",  "#C0DAE0",  "#F262BA",  "#B785BF", 
                "#DFDBF0",  "#918DC0",  "#FF6B8F",  "#E6B7D7",  "#B597B7",  "#3B95C8",  "#977195",  "#F88ED0",  "#F8BEEF", 
                "#68C589",  "#E0EEBB" ,  "#00BFFF", "#AB82FF", "#CDAA7D", "#FF7256", "#66CD00", "#6495ED", "#FFB90F", "#BF3EFF", "#EE1289", "#00BFFF", "#EE2C2C", "#228B22",
                "#ADFF2F", "#EE9A49", "#EE6AA7", "#FF6A6A", "#EE5C42", "#FFF68F", "#EE7942", "#FF0000", "#FFA500", "#EED2EE", "#FF4500", "#FFFF00",
                "#32CD32", "#4876FF" ,"#2E8B57", "#87CEEB", "#87CEFA", "#A020F0", "#6A5ACD" ,"#DDA0DD" ,"#FFC0CB", "#DB7093", "#AFEEEE", "#FF83FA")



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

##调整工作目录为函数脚本所放目录
setwd('/home/SPT')
## load('object.RData') ###针对Rdata这种，read针对Rds等

### 数据读取与保存
# saveRDS(object, 'object.RDS')
#读取
object <- readRDS('object.RDS')


##读取数据
object = Load10X_Spatial(
  data.dir = '/home/singlecell/rawdata/A/object/outs',
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper = FALSE
)

##作图，初步看下数据分布
plot1 = VlnPlot(object, features = "nCount_Spatial", group.by = 'TILS',pt.size = 0.1) + NoLegend()
plot2 = SpatialFeaturePlot(object, features = "nCount_Spatial") + theme(legend.position = "right")
plot_grid(plot1, plot2)

###引入TILS的变量
group_list <- read.csv('object-TLS EPI.csv',header = T,row.names = 1)
head(group_list)
head(object@meta.data)
object@meta.data = object@meta.data[,c(1:12)] ###删除多余的变量
object <- AddMetaData(object, group_list)
#查看TILS变量
table(object$TILS)
#降维前画门方式
plot3 = VlnPlot(object, features = "nCount_Spatial", group.by = 'TILS',pt.size = 0.1) + NoLegend()

###引入TISSUE的变量
group_list2 <- read.csv('object-TISSUE.csv',header = T,row.names = 1)
head(group_list2)
head(object@meta.data) 
object <- AddMetaData(object, group_list2)
#查看TISSUE变量
table(object$TISSUE)
#降维前画门方式
plot4 = VlnPlot(object, features = "nCount_Spatial", group.by = 'TISSUE',pt.size = 0.1) + NoLegend()
plot4
plot_grid(plot3, plot4)


########################################################################
##########################引用单细胞测序降维方式########################
#########################数据预处理
###推荐sctransform
object <- SCTransform(object, assay = "Spatial", verbose = FALSE)

###UMAP降维
object <- RunPCA(object, assay = "SCT", verbose = FALSE) %>% #运行PCA运算
  FindNeighbors(reduction = "pca", dims = 1:30) %>% #选择dims
  FindClusters(verbose = FALSE) %>% #选择resolution
  RunUMAP(reduction = "pca", dims = 1:30) #tsne运算
#看降维效果
p1 <- DimPlot(object, reduction = "umap", label = TRUE, cols = colors_aa)
p2 <- SpatialDimPlot(object, label = TRUE, label.size = 3, cols = colors_zzm)
p1
p2
p3 <- DimPlot(object, reduction = "umap", label = TRUE, split.by = 'TISSUE', cols = colors_aa)
p3
p4 <- SpatialDimPlot(object, label = TRUE, label.size = 3, 
                     group.by = "TISSUE", cols.highlight = c("#DE2D26", "grey50"),
                     pt.size.factor = 1.0)
p4 
p5 <- SpatialDimPlot(object, label = TRUE, label.size = 3, 
                     group.by = "TISSUE", cols.highlight = c("#DE2D26", "grey50"),
                     pt.size.factor = 0, cols = colors_zzm)
p5

#######看亚群的基因表达
#计算确定
DimPlot(object,reduction = 'tsne',cols = colorful[[1]])
object_markers = FindAllMarkers(object)

dim(object_markers)
table(object_markers$cluster)


Get_percentplot2(object,filename = '20210612object')
TLS.object2 <- FindVariableFeatures(TLS.object2)
TLS.object2@assays$RNA@var.features


DefaultAssay(TLS.object2) <- 'RNA'
TLS.object2@assays$SCT
TLS.object2=GSVA_score_insert(TLS.object2, gspath="/home/yankai/0007_SPT/20210624reanalysis/select_cd8t.gmt",
                            species="human",min.sz=2,parallel.sz=10)

library(Seurat)
TLS.object2=FindVariableFeatures(TLS.object2)
TLS.object2 = TLS.object2[which(rownames(TLS.object2)%in%VariableFeatures(TLS.object2)),]
rownames(TLS.object2)

total_matrix=as.data.frame(TLS.object2@assays$RNA@data)
#写表
write.csv(object_markers,'20210612object.csv')
#各分群的比例
library(LIANLAB)
COLOR <- data.frame(celltype = levels(object),cols = colors_zzm[1:11])
COLOR$cols <- as.character(COLOR$cols)
Get_percentplot(object,filename = '20210615object',colorm = COLOR)

###根据病理信息分组
object$tissue_cluster <- paste0(object$TISSUE,'_',object$seurat_clusters)
SpatialDimPlot(object, label = TRUE, label.size = 3,pt.size.factor = 1.8, 
               group.by = 'tissue_cluster', cols = colors_aa)

###根据病理信息分组后看亚群的基因表达
object@active.ident <- as.factor(object$seurat_clusters) ###替换活动表
object_markers <- FindAllMarkers(object)
#### object_markers = myfindmarkers(object,filename = '20210615object_markers_tissue_cluster')
#### FindMarkers(object,group.by = 'tissue_cluster')
dim(object_markers)
#写表
write.csv(object_markers,'20210805_object_cluster.csv')
#各分群的比例
library(LIANLAB)
COLOR <- data.frame(celltype = levels(object),cols = colors_aa[1:11])
COLOR$cols <- as.character(COLOR$cols)
## Get_percentplot_split(object,filename = '20210805_object_cluster',split_by = "tissue_cluster",colorm = COLOR)
Get_percentplot(object,filename = '20210805_object_cluster',colorm = COLOR)

###计算头部基因
library(LIANLAB,lib.loc = "/home/R/x86_64-pc-linux-gnu-library/4.0")
top_object_markers <- top_m(object_markers,20)

write.csv(top_object_markers,'20210805 object_cluster_top_object_markers.csv')