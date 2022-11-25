score_RNA_Clinical <- function(x,cutoff_immune,cutoff_tumor,cutoff_Fc){
  spt_tumor = x[1]
  spt_fc = x[2]
  spt_fc2 = x[3]
  spt_afc = x[4]
  spt_immune = x[5]
  G2=7.863990+2.136296*spt_tumor+14.159850*spt_fc+3.706542*spt_fc2-3.244212*spt_afc+7.876117*spt_immune #FC
  G3=-3.691816+2.931509*spt_tumor-1.499818*spt_fc-0.261779*spt_fc2+1.336733*spt_afc-1.243596*spt_immune #AFC
  G4=4.700847+3.716903*spt_tumor+9.680457*spt_fc+5.703404*spt_fc2-4.806712*spt_afc+2.869975*spt_immune #FC2
  G5=0 #other
  
  P_FC=exp(G2)/(exp(G2)+exp(G3)+exp(G4)+exp(G5))
  P_AFC=exp(G3)/(exp(G2)+exp(G3)+exp(G4)+exp(G5))
  P_FC2=exp(G4)/(exp(G2)+exp(G3)+exp(G4)+exp(G5))
  
  score <- c(P_FC,P_AFC,P_FC2)
  names(score) <- c('FC','AFC','FC2')  #不把immune放进去 immune的分在第一轮用过了
  RNA_Clinical_Type2 = names(score[which(score == max(score))])
  
  XY_tumor=5.485692-3.067855*spt_tumor+1.743714*spt_fc+1.624768*spt_fc2+0.817176*spt_immune #Tumor
  # XY_tumor=3.698178-2.022181*spt_tumor+2.307638*spt_fc+0.901418*spt_fc2+1.554894*spt_afc+1.329410*spt_immune #Tumor
  XY_immune=33.454-4.906*spt_immune+5.690*spt_tumor+15.310*spt_fc+7.938*spt_afc #immune
  
  if (XY_immune < cutoff_immune ) {
    RNA_Clinical_Type <- 'I' 
  }else{
    RNA_Clinical_Type <- 'nI' #not Immune的意思
  } #这是标签2
  
  if (XY_tumor < cutoff_tumor ) {
    RNA_Clinical_Type3 <- 'T' #not Immune的意思
  }else{
    RNA_Clinical_Type3 <- 'nT'
  } #这是标签3
  
  if (max(score) < cutoff_Fc) {
    RNA_Clinical_Type2 <- 'NA'
  } #这是标签1
  
  return(paste0(RNA_Clinical_Type3,'_',RNA_Clinical_Type2,'_',RNA_Clinical_Type))
}

rownames(SPTA2@assays$TYPE3@data)
SPTA2_spt_score <- as.data.frame(t(SPTA2@assays$TYPE3@data[1:5,])) #1-5对应 "spt.tumor" "spt.FC" "spt.FC2" "spt.AFC" "spt.immune"
head(SPTA2_spt_score)

# RNA_Clinical_Type_calculation <- apply(SPTA2_spt_score, 1, function(x) score_RNA_Clinical(x,cutoff_immune = 0.872,cutoff_other = 0.5) ) # 算分  阈值
RNA_Clinical_Type_calculation <- apply(SPTA2_spt_score, 1, function(x) score_RNA_Clinical(x,cutoff_immune = 0.5,cutoff_Fc = 0.5,cutoff_tumor = 0.5) ) # 算分  阈值

merged_tissue <- as.data.frame(RNA_Clinical_Type_calculation)

head(merged_tissue)
table(merged_tissue$RNA_Clinical_Type_calculation)  # 结果总览

SPTA2 <- AddMetaData(SPTA2,merged_tissue) # 结果回输
SpatialDimPlot(SPTA2,group.by = "RNA_Clinical_Type_calculation",cols = colors_figure)