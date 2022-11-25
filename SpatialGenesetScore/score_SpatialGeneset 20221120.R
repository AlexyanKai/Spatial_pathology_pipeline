######***********************######
######**根据基因列表计算打分 **######
######***********************######
# genelist <- read.table('/home/yankai/0007_SPT/20210918_SPTA_rename2_subset_markers_top20.txt',header = T,fill = T)
genelist <- read.table('E:/0000_SPT/20210919 spt and sRNA gene TOP20(1) with FC2.csv',header = T,fill = T,sep = ",")
genelist
merged_data <- as.data.frame(SPTE@assays$SCT@scale.data)
socre_data <- data.frame()
for (i in colnames(genelist)) {
  a <- merged_data[which(rownames(merged_data)%in%genelist[,i]),]
  dim(a)    
  socre <- data.frame(row.names = colnames(a),socres = apply(a, 2, 'mean'))
  colnames(socre) <- i
  if (i==colnames(genelist)[1]) {
    socre_data <- socre
  }else{
    socre_data <- cbind(socre_data,socre)
  }
}
SPTE[['TYPE']] <- CreateAssayObject(data = as.matrix(t(socre_data)))
rownames(SPTE@assays$TYPE@data)
#密度图
rownames(SPTE@assays$TYPE@data)
b <- as.data.frame(t(as.data.frame(SPTE@assays$TYPE@data)))
b$seurat_clusters <- SPTE$seurat_clusters
colnames(b)

ggplot(b,aes(x = spt.tumor)) + geom_density(aes(color = seurat_clusters,fill = seurat_clusters,alpha = 0 ))  + theme_bw()
