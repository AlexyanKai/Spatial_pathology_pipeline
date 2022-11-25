######Load libraries###########
library(Seurat)
###############################

######Set work dir#############
setwd('E:/A--Lianlabwork/A020--YkSpt20211228')
###############################

######Load 4 sample ST data###
load('E:/A--Lianlabwork/A018--YkSpt20211206/Step02--ContourAnalysis/Step00_SpData.RData')
ls()
a.srtobj<-a.object_rename2
rm(a.object_rename2,b.object,c.object,d.object)
load('20220105_bcd.prediction.RData')
ls()
#a\b\c\d.srtobj
##############################

###############################################
#Try method introduced in
#http://www.360doc.com/content/21/0416/18/19913717_972654728.shtml
library(tidyverse)
library(ggdendro)
library(cowplot)
library(ggtree)
library(patchwork)
######Try example######################
# gene_cluster<-read.table('scRNA_dotplot_data.tsv.gz',sep='\t',header=TRUE)
# gene_cluster %>% sample_n(5)
# gene_cluster %>% mutate('exppec'=(cell_exp_ct/cell_ct)*100) %>%
#   ggplot(aes(x=cluster,y=Gene,color=count,size=exppec))+geom_point()
###############################################
#Make a table like this format
#Tx Sample Pct AvgPbt
#T1 A 0.14  1
#T1 B 0.36  0.7
#Pct means percentage
#Pbt means probability, a=100
#It will has Tx=13*Samples4=52 rows
prediction.table<-data.frame(matrix(nrow=0,ncol=4))
Txs<-paste0('T',seq(1,13))
a.txs.spots<-colnames(a.srtobj)[which(a.srtobj@meta.data$TnT_cluster20211031 %in% Txs)]
b.txs.spots<-colnames(b.srtobj)[which(b.srtobj@meta.data$predicted.id %in% Txs)]
c.txs.spots<-colnames(c.srtobj)[which(c.srtobj@meta.data$predicted.id %in% Txs)]
d.txs.spots<-colnames(d.srtobj)[which(d.srtobj@meta.data$predicted.id %in% Txs)]

for (Tx in Txs){
  a.tx.spots<-colnames(a.srtobj)[which(a.srtobj@meta.data$TnT_cluster20211031==Tx)]
  a.pct<-100.0*length(a.tx.spots)/length(a.txs.spots)
  a.avgpbt<-100
  prediction.table<-rbind(prediction.table,c(Tx,'A',a.pct,a.avgpbt))
  
  prediction.colname<-paste0('prediction.score.',Tx)
  
  if (Tx %in% b.srtobj@meta.data$predicted.id){
    b.tx.spots<-colnames(b.srtobj)[which(b.srtobj@meta.data$predicted.id==Tx)]
    b.pct<-100.0*length(b.tx.spots)/length(b.txs.spots)
    b.avgpbt<-mean(b.srtobj@meta.data[b.tx.spots,prediction.colname])*100
    prediction.table<-rbind(prediction.table,c(Tx,'B',b.pct,b.avgpbt))
  }else{
    prediction.table<-rbind(prediction.table,c(Tx,'B',0,0))
  }
  
  if (Tx %in% c.srtobj@meta.data$predicted.id){
    c.tx.spots<-colnames(c.srtobj)[which(c.srtobj@meta.data$predicted.id==Tx)]
    c.pct<-100.0*length(c.tx.spots)/length(c.txs.spots)
    c.avgpbt<-mean(c.srtobj@meta.data[c.tx.spots,prediction.colname])*100
    prediction.table<-rbind(prediction.table,c(Tx,'C',c.pct,c.avgpbt))
  }else{
    prediction.table<-rbind(prediction.table,c(Tx,'C',0,0))
  }
  
  if (Tx %in% d.srtobj@meta.data$predicted.id){
    d.tx.spots<-colnames(d.srtobj)[which(d.srtobj@meta.data$predicted.id==Tx)]
    d.pct<-100.0*length(d.tx.spots)/length(d.txs.spots)
    d.avgpbt<-mean(d.srtobj@meta.data[d.tx.spots,prediction.colname])*100
    prediction.table<-rbind(prediction.table,c(Tx,'D',d.pct,d.avgpbt))
  }else{
    prediction.table<-rbind(prediction.table,c(Tx,'D',0,0))
  }
}
colnames(prediction.table)<-c('Tx','Sample','Pct','AvgPbt')
prediction.table$Pct<-as.numeric(prediction.table$Pct)
prediction.table$AvgPbt<-as.numeric(prediction.table$AvgPbt)
prediction.table%>%ggplot(aes(x=Sample,y=Tx,color=AvgPbt,size=Pct))+geom_point()

pp<-ggplot(prediction.table,aes(x=Sample,y=Tx,color=AvgPbt,size=Pct),)+geom_point()+scale_color_gradient(low = "cyan",high = "red")
pp + theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(),axis.line = element_line(colour = "black"))
