######Set word dir and load data######
setwd('E:/A--Lianlabwork/A018--YkSpt20211206/Step02--ContourAnalysis')
#setwd('/media/complxmat/LabWorkMain/A--Lianlabwork/A018--YkSpt20211206/Step02--ContourAnalysis')
load('Step00_SpData.RData')
ls()
#"a.object_rename2" "b.object"         "c.object"         "d.object"
slotNames(a.object_rename2)
# [1] "assays"       "meta.data"    "active.assay" "active.ident" "graphs"       "neighbors"   
# [7] "reductions"   "images"       "project.name" "misc"         "version"      "commands"    
# [13] "tools"
names(a.object_rename2@meta.data)
######################################


#############Load libraries###########
library(Seurat)
library(CellChat)
######################################

#######Define function for generating ligand2receptor score######
#This function uses a global variable gene_exp
l2r_score_calculation<-function(l2r) {
  genes<-unlist(strsplit(l2r,'_'))
  if (genes[1] %in% rownames(gene_exp) & genes[2] %in% rownames(gene_exp))
    l2r_score<-gene_exp[genes[1],]+gene_exp[genes[2],]
  else
    l2r_score<-rep(0,dim(gene_exp)[2])
  names(l2r_score)<-colnames(gene_exp)
  return(l2r_score)
}
#################################################################

############Load CellChatDB for Ligand2Receptor score############
options(stringsAsFactors = FALSE)
CellChatDB <- CellChatDB.human #if working on the human dataset
interaction_input <- CellChatDB$interaction
ligand2receptors<-interaction_input$interaction_name#Ligand_Receptor1_Receptor2_...
#length(ligand2receptors)
#Convert to ligand2receptor format
ligand2receptor_list<-c()
for (l2rs in ligand2receptors){
  genes<-unlist(strsplit(l2rs,'_'))
  l2r_list<-paste(genes[1],genes[-1],sep='_')#Ligand_Receptor1,Ligand_Receptor2,...
  for (l2r in l2r_list){
    if (l2r %in% ligand2receptor_list)
      next
    else
      ligand2receptor_list<-append(ligand2receptor_list,l2r)
  }
}
#length(ligand2receptor_list)
################################################################

######Extract gene exp and spot coordinates from each seurat object
######Start with a.object_rename2
#Read Subjective defined region type data
subj_type<-read.csv('Step00_SubjectiveDefinedTypes/20210918 SPTA.csv',row.names = 1)
a.object_rename2[['SubjType']]<-subj_type
a.object_rename2[['SubjType']][is.na(a.object_rename2[['SubjType']])]<-'NotDefined'
table(a.object_rename2@meta.data$SubjType)
# #
# ATPICAL CELL FOLLICULAR CELL          immune      NotDefined           tumor 
# 128             584              70            1563              11 
# #
Idents(a.object_rename2)<-'SubjType'
#Extract Image coordinates
coordi_matrix<-a.object_rename2@images$slice1@coordinates
coordi_matrix['SubjType']<-a.object_rename2[['SubjType']]
write.table(coordi_matrix,'Step01_CoordinateAndCelltype/A_CoordinateAndSubjType.txt',sep='\t',quote=FALSE)
#Extract Image coordinates with altered identities TnT
table(a.object_rename2[['TnT_20211031']])
# AFC    AFCI     FC1    FC1I     FC2    FC2I       I nTnFCnI       T 
# 501     156     707      76     267     132     269      85     163
coordi_matrix<-a.object_rename2@images$slice1@coordinates
coordi_matrix['TnT20211031']<-a.object_rename2[['TnT_20211031']]
write.table(coordi_matrix,'Step01_CoordinateAndCelltype/A_CoordinateAndTnT.txt',sep='\t',quote=FALSE)
#Extract Gene Exp to file
write.table(a.object_rename2@assays$RNA@scale.data,'Step01_GeneScaleData/A_GeneExpScaled.txt',sep='\t',quote=FALSE)
#Calculate Ligand2Receptor score and write to file
gene_exp<-a.object_rename2@assays$RNA@scale.data
l2r_score_matrix<-t(as.data.frame(sapply(ligand2receptor_list,l2r_score_calculation)))
l2r_score_matrix<-l2r_score_matrix[which(rowSums(l2r_score_matrix) > 0),]
write.table(l2r_score_matrix,'Step01_L2rScore/A_L2rSumScore.txt',sep='\t',quote=FALSE)

######b.object
write.table(b.object@assays$SCT@scale.data,'Step01_GeneScaleData/B_GeneExpScaled.txt',sep='\t',quote=FALSE)
subj_type<-read.csv('Step00_SubjectiveDefinedTypes/20210918 SPTB.csv',row.names = 1)
b.object[['SubjType']]<-subj_type
b.object[['SubjType']][is.na(b.object[['SubjType']])]<-'NotDefined'
table(b.object@meta.data$SubjType)
# blood           FIBRO FOLLICULAR CELL      NotDefined           tumor 
# 6             419             932            1902             490 
Idents(b.object)<-'SubjType'
coordi_matrix<-b.object@images$slice1@coordinates
coordi_matrix['SubjType']<-b.object[['SubjType']]
write.table(coordi_matrix,'Step01_CoordinateAndCelltype/B_CoordinateAndSubjType.txt',sep='\t',quote=FALSE)
gene_exp<-b.object@assays$SCT@scale.data
l2r_score_matrix<-t(as.data.frame(sapply(ligand2receptor_list,l2r_score_calculation)))
l2r_score_matrix<-l2r_score_matrix[which(rowSums(l2r_score_matrix) > 0),]
write.table(l2r_score_matrix,'Step01_L2rScore/B_L2rSumScore.txt',sep='\t',quote=FALSE)

######c.object
write.table(c.object@assays$SCT@scale.data,'Step01_GeneScaleData/C_GeneExpScaled.txt',sep='\t',quote=FALSE)
subj_type<-read.csv('Step00_SubjectiveDefinedTypes/20210918 SPTC.csv',row.names = 1)
c.object[['SubjType']]<-subj_type
c.object[['SubjType']][is.na(c.object[['SubjType']])]<-'NotDefined'
table(c.object@meta.data$SubjType)
# FIBRO FOLLICULAR CELL      NotDefined 
# 14             195            1109
Idents(c.object)<-'SubjType'
coordi_matrix<-c.object@images$slice1@coordinates
coordi_matrix['SubjType']<-c.object[['SubjType']]
write.table(coordi_matrix,'Step01_CoordinateAndCelltype/C_CoordinateAndSubjType.txt',sep='\t',quote=FALSE)
gene_exp<-c.object@assays$SCT@scale.data
l2r_score_matrix<-t(as.data.frame(sapply(ligand2receptor_list,l2r_score_calculation)))
l2r_score_matrix<-l2r_score_matrix[which(rowSums(l2r_score_matrix) > 0),]
write.table(l2r_score_matrix,'Step01_L2rScore/C_L2rSumScore.txt',sep='\t',quote=FALSE)

######d.object
write.table(d.object@assays$SCT@scale.data,'Step01_GeneScaleData/D_GeneExpScaled.txt',sep='\t',quote=FALSE)
subj_type<-read.csv('Step00_SubjectiveDefinedTypes/20210918 SPTD.csv',row.names = 1)
d.object[['SubjType']]<-subj_type
d.object[['SubjType']][is.na(d.object[['SubjType']])]<-'NotDefined'
table(d.object@meta.data$SubjType)
# blood           FIBRO FOLLICULAR CELL          immune      NotDefined           tumor 
# 52              27             115              27            1791              19 
Idents(d.object)<-'SubjType'
coordi_matrix<-d.object@images$slice1@coordinates
coordi_matrix['SubjType']<-d.object[['SubjType']]
write.table(coordi_matrix,'Step01_CoordinateAndCelltype/D_CoordinateAndSubjType.txt',sep='\t',quote=FALSE)
gene_exp<-d.object@assays$SCT@scale.data
l2r_score_matrix<-t(as.data.frame(sapply(ligand2receptor_list,l2r_score_calculation)))
l2r_score_matrix<-l2r_score_matrix[which(rowSums(l2r_score_matrix) > 0),]
write.table(l2r_score_matrix,'Step01_L2rScore/D_L2rSumScore.txt',sep='\t',quote=FALSE)

######Manual check featureplot from A
SpatialFeaturePlot(a.object_rename2,'TMSB4X')
