# load packages
library(Seurat)

#######***********************######
######**Cell Label Transfer**######
######***********************######

########以A样本为准cell label transfer
table(object$group) #此数据正确
SpatialDimPlot(object2)
to.b.object<-FindTransferAnchors(reference=object,query=object2, 
                                 normalization.method = "SCT",dims=1:10)
b.predictions<-TransferData(anchorset=to.b.object,refdata = object$T_Tsub_FC_I_20220426,dims=1:10)
b.srtobj<-AddMetaData(object2,metadata = b.predictions)
table(b.srtobj$predicted.id)
SpatialDimPlot(b.srtobj,group.by = "predicted.id",cols = colorful[[1]],pt.size.factor = 1.6)
head(b.srtobj@meta.data)

#学术作图
pdf(paste0("filename.pdf"))
SpatialDimPlot(b.srtobj,group.by = "predicted.id",cols = colors_figure5EF,crop = F,pt.size.factor = 1.35)
dev.off()
object2 <- b.srtobj