.libPaths("/home/complxmat/R/x86_64-pc-linux-gnu-library/4.1")
setwd("/media/sftp_root/shared/yankai")
#sftp://yankai@192.168.31.250/shared

# load packages
library(SPATA2)
library(tidyverse)
library(Seurat)

# object <- readRDS("file.RDS")

object <-transformSeuratToSpata(object,"SPT",method = "spatial",coords_from = "umap", 
                                 assay_name = "SCT", assay_slot = NULL, image_name = NULL, 
                                 gene_set_path = NULL, verbose = TRUE)
str(object,max.level = 2)
head(object@fdata)
# if you want to save a spata object from your session use
# saveSpataObject(object = object,file_name = "object") 
getGeneSets(object)
length(unique(object@used_genesets$ont))
head(unique(object@used_genesets$ont))
getGeneSetsInteractive(object)

##
plotSurface(object = object,
            of_sample = "SPT",
            color_by = "predicted.id_T",
            pt_clrp = "npg",
            pt_size = 1.5) + labs(color = "Clusters")
plotSurface(object = object, color_by = "seurat_clusters")
plotSurface(object = object, color_by = "group1")
# plot gene-set expression 
plotSurface(object = object,
            of_sample = "SPT",
            color_by  = "TFF3",
            pt_size = 1.5,
            pt_clrsp = "magma")
# plot gene-set expression (spatially smoothed)
plotSurface(object = object,
            color_by =  "TFF3",
            pt_size = 1.6,
            pt_clrsp = "magma",
            smooth = TRUE,
            smooth_span = 0.01)

####
# ?plotSurfaceComparison
# Character value. The method according to which gene sets will be handled specified as a character of length one. This can be either 'mean' or one of 'gsva', 'ssgsea', 'zscore', or 'plage'. The latter four will be given to gsva::GSVA().
object <- createSegmentation(object = object)
plotSegmentation(object)
str(object)
object@fdata$SPT$segmentation
object@fdata$SPT$group1
object[['tumor_segment']]=paste(object@fdata$SPT$group1,object@fdata$SPT$segmentation,sep='__')
saveSpataObject(object = object,file = "20211008_object.RDS") 
saveSpataObject(object)
object.obj <- transformSpataToSeurat(object)

