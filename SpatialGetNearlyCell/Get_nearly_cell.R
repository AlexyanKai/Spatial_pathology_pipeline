#----最短距离-------
Get_nearly_cell <- function(SPT.obj,types = 'Coordinates',
                            cluster.1 = NULL,cluster.2 = NULL,
                            cell.1 = NULL,cell.2 = NULL,
                            point.size = 3, cols=c('red','blue'), 
                            draw_line = T,
                            line.size = 0.02,line.cols = 'black'){
  
  if (types == 'Coordinates') {
    xy_info <- GetTissueCoordinates(SPT.obj)  
  }else{
    if (types == 'spotcenters') {
      xy_info <- SPT.obj@images$slice1@coordinates
      colnames(xy_info)
      xy_info <- xy_info[,c("row","col" )]
    }else{
      if (types =='imagerow') {
        xy_info <- SPT.obj@images$slice1@coordinates
        colnames(xy_info)
        xy_info <- xy_info[,c("imagerow","imagecol")]
      }else{
        print('Please select a mode.')
        stop()
      }
    }
  }
  
  
  dim(xy_info)
  if (!is.null(cluster.1)) {
    cell.1 <- rownames(subset(SPT.obj@meta.data,seurat_clusters==cluster.1 ))
    cell.2 <- rownames(subset(SPT.obj@meta.data,seurat_clusters==cluster.2 ))
  }
  type1 <- xy_info[which(rownames(xy_info)%in%cell.1),]
  type2 <- xy_info[which(rownames(xy_info)%in%cell.2),]
  
  
  plan(strategy = 'multisession',workers = 20)
  i = 1
  system.time(  res <- future_lapply(1:nrow(type1),function(i){
    j = 1
    while (j < nrow(type2) + 1) {
      if (j==1) {
        distancess <- sqrt((type1[i,1] - type2[j,1])^2 + (type1[i,2] - type2[j,2])^2)
        barcodes <- rownames(type2[j,])
      }else{
        if (distancess > sqrt((type1[i,1] - type2[j,1])^2 + (type1[i,2] - type2[j,2])^2) ) {
          distancess = sqrt((type1[i,1] - type2[j,1])^2 + (type1[i,2] - type2[j,2])^2)
          barcodes <- rownames(type2[j,])
        }
      }
      j = j + 1
    }
    
    x <- data.frame(row.names = c('near_cell_barcodes','distances'), c(barcodes,distancess))
    colnames(x) <- rownames(type1[i,])
    return(x)
  }
  )
  )
  
  plan(strategy = 'multisession',workers = 1)
  merged_distancess <- data.frame()
  for (i in 1:length(res)) {
    if (i == 1) {
      merged_distancess <- res[[i]]
    }else {
      merged_distancess <- cbind(merged_distancess,res[[i]])
    }
    
  }
  merged_distancess <- as.data.frame(t(merged_distancess))
  head(merged_distancess)
  
  type1$type <- rep('cluster1',nrow(type1))
  type2$type <- rep('cluster2',nrow(type2))
  draw_basic_info <- rbind(type1,type2)
  colnames(draw_basic_info)[1:2] <- c('imagerow','imagecol') 
  
  merged_distancess$barcode <- rownames(merged_distancess)
  type1$barcode <- rownames(type1)
  
  draw_info <- merge(type1,merged_distancess,by = 'barcode')
  
  colnames(draw_info)[2:3] = c('imagerow.1','imagecol.1')
  
  draw_info <- draw_info[,-4]
  select_type2 <- type2[which(rownames(type2)%in%draw_info$near_cell_barcodes),]
  colnames(select_type2)[1:2] = c('imagerow.2','imagecol.2')
  select_type2 <- select_type2[,-3]
  select_type2$near_cell_barcodes <- rownames(select_type2)
  
  draw_info <- merge(draw_info,select_type2, by = 'near_cell_barcodes',all.x = T)
  dim(draw_info)
  
  head(draw_basic_info)
  head(draw_info)
  
  p <- ggplot(draw_basic_info,aes(x = imagerow, y = imagecol,color = type, fill = type)) + geom_point(size = point.size) + 
    scale_color_manual(values = cols) + theme_bw() 
  print(p)
  
  if (draw_line == T){
    for (i in 1:nrow(draw_info)) {
      p <- p + geom_segment(x=draw_info$imagerow.1[i],y=draw_info$imagecol.1[i],
                            xend=draw_info$imagerow.2[i],yend = draw_info$imagecol.2[i], size = line.size, color = line.cols)
    }
  }
  print(p)
  
  pdf('cell_near_cell.pdf',width = 10,height = 10)
  print(p)
  dev.off()
  
  
  return(merged_distancess)
  
}