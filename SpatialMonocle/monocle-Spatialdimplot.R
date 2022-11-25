library(LIANLAB)
library(Seurat)
library(GSVA)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(pheatmap)
library(limma)
library(biomaRt)
library(monocle3)
library(SPATA2)

source('E:/011 BioInformation/singlecell_pipeline/singlecell_pepiline/singlecell_basic_analyze.R')
setwd("E:/0000 空间转录组/20220328 SPT Analysis-PTC")
save.image("20220405_object.RData")

#20220405
#################################################################
########################筛选合适的assay标签#############################
#################################################################
SpatialDimPlot(object)
table(object@meta.data$seurat_clusters) 

#对此数据进行作图
table(object@meta.data$seurat_clusters)

#
#重命名,将FC Area当一个整体分析
levels(object@active.ident)
object@active.ident <- as.factor(object$seurat_clusters)
levels(object@active.ident)
table(object@active.ident)
new.cluster.ids <- c(
  "FC Area",  "FC Area","FC Area","FC Area","FC Area","FC Area","Immune","NA","Tumor Area"
)
names(x = new.cluster.ids) <- levels(x = object)
object <- RenameIdents(object = object, new.cluster.ids)
levels(object@active.ident)
SpatialDimPlot(object)

SpatialDimPlot(object,cols = colors)
SpatialDimPlot(object,cols = colors_spt)
SpatialDimPlot(object,cols = colors1)

#################################################################
########################Dimplot/SPATA2/monocle3#############################
#################################################################
object_TFCAFC_spata <- transformSeuratToSpata(object_TFCAFC,sample_name = "SPTA",assay_name = "SCT")
# compile a cell-data-set
cortex_cds <- transformSpataToCDS(object = object_TFCAFC_spata)

#要跟踪不同平台的相应对象，您可以将保存它们的目录存储在 spata-object 中。要获取当前存储的目录，请使用getDirectoryInstructions(). 通过 将目录添加到相应的 cell_data_set 中adjustDirectoryInstruction()。
object_TFCAFC_spata2 <- adjustDirectoryInstructions(object = object_TFCAFC_spata,
                                                              to = "cell_data_set",
                                                              directory_new = "E:/0000 空间转录组/20220328 SPT Analysis-PTC/object_TFCAFC_spata_cell_data_set")

# obtain the directories that are currently stored in the spata-object
getDirectoryInstructions(object = object_TFCAFC_spata2)

#要在存储在 spata-object 中的默认目录下保存或加载相应的 cell_data_set 使用saveCorrespondingCDS()和loadCorrespondingCDS()
# use the directory information of the spata-object to conveniently save the cell_data_set
saveCorrespondingCDS(cds = cortex_cds, object = object_TFCAFC_spata2)

# use the directory information of the spata-object to conveniently load the cell_data_set
cortex_cds <- loadCorrespondingCDS(object = object_TFCAFC_spata2) 

plotSurface(object_TFCAFC_spata2, color_by = "seurat_clusters", pt_clrp = "npg", pt_size = 1.4) + 
  legendBottom()

plot_cells(cortex_cds,  reduction_method = "UMAP",
           color_cells_by = "seurat_clusters" # visualize a transferred variable
) +  legendNone() 

#作图
cell_size = 2
pdf(paste0("20220108_object_TFCAFC_mococle——plotcells.pdf"),width=10,height=8)
plot_cells(
  cortex_cds,  x = 1,  y = 2,
  reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
  color_cells_by = "seurat_clusters",
  group_cells_by = c("cluster"),
  genes = NULL,
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "grey28",
  trajectory_graph_segment_size = 0.75,
  norm_method = c( "size_only"),
  label_cell_groups = F,  label_groups_by_cluster = TRUE,
  group_label_size = 5,  labels_per_group = 21,
  label_branch_points = TRUE,  label_roots = TRUE,
  label_leaves = TRUE,  graph_label_size = 2,
  alpha = 1,  min_expr = 0.1,  cell_size = 1,
  cell_stroke = I(cell_size/2),
  rasterize = FALSE,  scale_to_range = FALSE,
  label_principal_points = FALSE
)
dev.off()

#################################################################
########################Dimplot/SPATA2/monocle3#############################
#################################################################
a.object_spata

# compile a cell-data-set
cortex_cds <- transformSpataToCDS(object = a.object_spata)

#要跟踪不同平台的相应对象，您可以将保存它们的目录存储在 spata-object 中。要获取当前存储的目录，请使用getDirectoryInstructions(). 通过 将目录添加到相应的 cell_data_set 中adjustDirectoryInstruction()。
a.object_spata2 <- adjustDirectoryInstructions(object = a.object_spata,
                                               to = "cell_data_set",
                                               directory_new = "E:/0000 空间转录组/20220328 SPT Analysis-PTC/object_TFCAFC_spata_cell_data_set")

# obtain the directories that are currently stored in the spata-object
getDirectoryInstructions(object = a.object_spata2)

#要在存储在 spata-object 中的默认目录下保存或加载相应的 cell_data_set 使用saveCorrespondingCDS()和loadCorrespondingCDS()
# use the directory information of the spata-object to conveniently save the cell_data_set
saveCorrespondingCDS(cds = cortex_cds, object = a.object_spata2)

# use the directory information of the spata-object to conveniently load the cell_data_set
# cortex_cds <- loadCorrespondingCDS(object = a.object_spata2) 

plotSurface(a.object_spata2, color_by = "seurat_clusters", pt_clrp = "npg", pt_size = 1.4) + 
  legendBottom()

plot_cells(cortex_cds,  reduction_method = "UMAP",
           color_cells_by = "seurat_clusters" # visualize a transferred variable
) +  legendNone() 

#################################################################
########################筛选合适的assay标签#############################
#################################################################
SpatialDimPlot(object)
table(object@meta.data$seurat_clusters) 
table(object@meta.data$Tumorcluster20211031)
table(object@meta.data$Seurat_clusters3)
table(object@meta.data$TnT_cluster20211031)
table(object@meta.data$TnT_Tcluster20211101)
table(object@meta.data$seurat_clusters)
table(object@meta.data$group9)


#作图
cell_size = 2
pdf(paste0("20220108_object_TFCAFC_mococle——plotcells.pdf"),width=10,height=8)
plot_cells(
  cortex_cds,  x = 1,  y = 2,
  reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
  color_cells_by = "seurat_clusters",
  group_cells_by = c("cluster"),
  genes = NULL,
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "grey28",
  trajectory_graph_segment_size = 0.75,
  norm_method = c( "size_only"),
  label_cell_groups = F,  label_groups_by_cluster = TRUE,
  group_label_size = 5,  labels_per_group = 21,
  label_branch_points = TRUE,  label_roots = TRUE,
  label_leaves = TRUE,  graph_label_size = 2,
  alpha = 1,  min_expr = 0.1,  cell_size = 1,
  cell_stroke = I(cell_size/2),
  rasterize = FALSE,  scale_to_range = FALSE,
  label_principal_points = FALSE
) + scale_color_manual( values=colors1)
dev.off()

colors1

#重命名,将FC Area当一个整体分析
levels(object@active.ident)
object@active.ident <- as.factor(object$seurat_clusters)
levels(object@active.ident)
table(object@active.ident)
new.cluster.ids <- c(
  "FC Area",  "FC Area","FC Area","FC Area","FC Area","FC Area","Immune","NA","Tumor Area"
)
names(x = new.cluster.ids) <- levels(x = object)
object <- RenameIdents(object = object, new.cluster.ids)
levels(object@active.ident)
SpatialDimPlot(object,cols = colors1,group.by = "seurat_clusters")

table(object@meta.data$seurat_clusters)

#作图
#FC area作为整体
cell_size = 2
pdf(paste0("20220405_object_mococle——plotcells.pdf"),width=10,height=8)
plot_cells(
  cortex_cds,  x = 1,  y = 2,
  reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
  color_cells_by = "seurat_clusters",
  group_cells_by = c("cluster"),
  genes = NULL,
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "grey28",
  trajectory_graph_segment_size = 0.75,
  norm_method = c( "size_only"),
  label_cell_groups = F,  label_groups_by_cluster = TRUE,
  group_label_size = 5,  labels_per_group = 21,
  label_branch_points = TRUE,  label_roots = TRUE,
  label_leaves = TRUE,  graph_label_size = 2,
  alpha = 1,  min_expr = 0.1,  cell_size = 1,
  cell_stroke = I(cell_size/2),
  rasterize = FALSE,  scale_to_range = FALSE,
  label_principal_points = FALSE
) + scale_color_manual( values=c("#6DA1F2","#6DA1F2","#6DA1F2","#6DA1F2","#6DA1F2","#6DA1F2","#D7EBF8","#FBE2DA","#E8696B"))
dev.off()

#作图
#FC area拆分
cell_size = 2
pdf(paste0("20220405_object_mococle—seurat_clusters—plotcells-4.pdf"),width=10,height=8)
plot_cells(
  cortex_cds,  x = 1,  y = 2,
  reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"),
  color_cells_by = "seurat_clusters",
  group_cells_by = c("cluster"),
  genes = NULL,
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "grey28",
  trajectory_graph_segment_size = 0.75,
  norm_method = c( "size_only"),
  label_cell_groups = F,  label_groups_by_cluster = TRUE,
  group_label_size = 5,  labels_per_group = 21,
  label_branch_points = TRUE,  label_roots = TRUE,
  label_leaves = TRUE,  graph_label_size = 2,
  alpha = 1,  min_expr = 0.1,  cell_size = 1,
  cell_stroke = I(cell_size/2),
  rasterize = FALSE,  scale_to_range = FALSE,
  label_principal_points = FALSE
) + scale_color_manual( values=c("#AB82FF","#6F56D4","#EE2C2C","#68C589", "#FED293", "#00B8E2","#D7EBF8","#FBE2DA","#E8696B")) #此色阶很重要，贯穿Figure2
dev.off()
# values=c("#AB82FF","#6F56D4","#977195","#68C589", "#FED293", "#00B8E2","#D7EBF8","#FBE2DA","#E8696B")
# values=c("#AB82FF","#6F56D4","#977195","#68C589", "#FED293", "#C4C4B4","#D7EBF8","#FBE2DA","#E8696B")
# values=c("#AB82FF","#6F56D4","#977195","#68C589", "#FED293", "#EE2C2C","#D7EBF8","#FBE2DA","#E8696B")
# values=c("#AB82FF","#6F56D4","#977195","#68C589", "#FED293", "#00BFFF","#D7EBF8","#FBE2DA","#E8696B")
# values=c("#AB82FF","#6F56D4","#977195","#68C589", "#FED293", "#67D9FB","#D7EBF8","#FBE2DA","#E8696B")
# values=c("#AB82FF","#6F56D4","#977195","#68C589", "#FED293", "#0B8E2","#D7EBF8","#FBE2DA","#E8696B")

#20220405 FC area拆分作图
table(object@meta.data$seurat_clusters)
table(object@active.ident)
object$TFCI_area <- as.factor(object@active.ident) #将活动变量存入TFCI_area
table(object$TFCI_area)
object@active.ident <- as.factor(object$seurat_clusters) #将seurat_clusters作为活动变量
table(object@active.ident)

#
SpatialDimPlot(object)
object_FC_area <- subset(object,idents =c("AFC", "AFCI",  "FC1","FC1I",  "FC2",  "FC2I"))
SpatialDimPlot(object_FC_area,cols = colors_aa)
pdf(paste0(i,"20220405_object_FC_area_SpatialDimPlot2.pdf"), width=9, height=8)
p=SpatialDimPlot(object_FC_area,cols = c("#AB82FF","#6F56D4","#977195","#68C589", "#FED293", "#00B8E2"),
                 crop = F,pt.size.factor = 1.35)
p=SpatialDimPlot(object_FC_area,cols = c("#AB82FF","#6F56D4","#EE2C2C","#68C589", "#FED293", "#00B8E2"),
                 crop = F,pt.size.factor = 1.35)
print(p)
dev.off()


