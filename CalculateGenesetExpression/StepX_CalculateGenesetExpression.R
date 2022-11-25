######Install qusage package to use read.gmt###
#BiocManager::install("qusage")
###############################################

######Load library######
library(qusage)
########################

######Set work dir######
setwd('E:/A--Lianlabwork/A018--YkSpt20211206/Step08--GenesetContour')
########################

######Read geneset
genesets<-read.gmt('E:/0000 空间转录组/20211122 SPT Analysis-PTC/GenesetContourMaps/c3.tft.v7.4.symbols.gmt')
length(genesets)
names(genesets)
genesets[1]
genesets[[(names(genesets)[2])]]
#################

######Load Seurat data######
load('E:/A--Lianlabwork/A018--YkSpt20211206/Step02--ContourAnalysis/Step00_SpData.RData')
#Copy the gene expression matrix to global variable gene.exp
gene.exp<-a.object_rename2@assays$SCT@scale.data
#Copy the genes in expression matrix to global variable genes.with.exp
genes.with.exp<-rownames(gene.exp)
#A function use global variable genes.with.exp to filter genesets
#Only genesets with at least gsN number of genes will be used
gsN=5
find_available_genes<-function(genes){
  ags<-intersect(genes,genes.with.exp)
  if(length(ags)<gsN){
    return(NA)
  }else{
    return(ags)
  }
}
available.genesets<-sapply(genesets,find_available_genes)
available.genesets<-available.genesets[!is.na(available.genesets)]
#The function calculating one geneset will use gene.exp
calculate_geneset_exp<-function(genes){
  sum.exp<-colSums(gene.exp[genes,])
  return(sum.exp)
}
geneset.expression.matrix<-sapply(available.genesets,calculate_geneset_exp)
geneset.expression.matrix<-t(geneset.expression.matrix)
write.csv(geneset.expression.matrix,"20220107_c3.tftpathway_geneset.expression.matrix.csv")
colnames(geneset.expression.matrix)

a.object_rename2[["c3tftpathway"]] <- CreateAssayObject(data = geneset.expression.matrix)
a.object_rename2@assays$c3tftpathway@data
head(a.object_rename2@assays$hallpathway@data)
a.object_rename2@assays$hallpathway@data['KEGG-PRIMARY-IMMUNODEFICIENCY',]
SpatialFeaturePlot(a.object_rename2,features = "WP-OXIDATIVE-PHOSPHORYLATION")
str(a.object_rename2)

###################清除多余Assay######################
a.object_rename2[['c2.cp_pathway']] <- NULL
saveRDS(a.object_rename2,"20210107_a.object_rename2_with_pathway.RDS")
