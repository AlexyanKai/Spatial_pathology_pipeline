#加载包，清空
library(Seurat)
library(tidyverse)
library(ggplot2)
library(infercnv)
library(ComplexHeatmap)
library(ggpubr)
setwd('./infercnv')

# rm(list=ls())
#加载数据
# scRNA<- readRDS("object.RDS")
scRNA<- object

scRNA_harmony <- scRNA
scRNA_harmony <- subset(scRNA_harmony,idents=c("T","AFC","FC1"))
#########################################################################################################################################################
#思考一下 如何推断cnv的  gene的表达量（count表达矩阵），gene的位置（基因染色体信息），比较信息（分组信息）
#inferCNV需要三个文件1.count表达矩阵，2.分组信息，3.基因染色体信息

#制作基因染色体位置信息 和提取表达矩阵
scRNA_harmony@assays$Spatial@counts['ISG15',1:3]
scRNA_harmony@assays$SCT@counts['ISG15',1:3]

dat <- GetAssayData(scRNA_harmony,assay = "Spatial",slot = "counts")
dat <- as.data.frame(dat)
dim(dat)
# expFile=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv")  #这个是inferCNV自带参比数据集  用来做infer
# geneFile=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv")   #这个是inferCNV自带参比数据集  用来做infer
# groupFiles=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv")    #这个是inferCNV自带参比数据集  用来做infer


# install.packages('AnnoProbe')

library(AnnoProbe)
geneInfor=annoGene(rownames(dat),"SYMBOL",'human')

colnames(geneInfor)

geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]      
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)

dat=dat[match(geneInfor[,1], rownames(dat)),]    #将表达矩阵的基因排列顺序与geneInfor的基因排列顺序弄成一致
rownames(geneInfor) <- geneInfor$SYMBOL   
geneInfor <- geneInfor[,-1]     #这样我们就制作好了染色体位置信息和排列顺序好的count表达矩阵

#制作mate信息
head(scRNA_harmony@meta.data)

scRNA_harmony$celltype <- as.character(scRNA_harmony@active.ident)

DefaultAssay(scRNA_harmony)
table(scRNA_harmony@active.ident)

DimPlot(scRNA_harmony)
head(scRNA_harmony@meta.data)

iso_cell <- sample(rownames(subset(scRNA_harmony@meta.data,celltype =='FC1')),200) #随机挑选 iso cell 
scRNA_harmony$celltype[which(colnames(scRNA_harmony)%in%iso_cell)] <- 'iso'
table(scRNA_harmony$celltype)

meta <- subset(scRNA_harmony@meta.data,select = c("celltype"))   #假如你后面是想分析每一个群的CNV的话  这里要改成seruat_cluster
table(meta$celltype)

#验证1   表达矩阵的列名要与meta的行名一致
identical(colnames(dat),rownames(meta))  
#验证2   表达矩阵的行名要与geneInfor的行名一致
identical(rownames(dat),rownames(geneInfor))

#因此三个输入数据准备好了   dat-表达矩阵  meta-分组信息  geneInfor-基因染色体信息


#二步法构建对象
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=dat,
                                    annotations_file=meta,
                                    delim="\t",
                                    gene_order_file=geneInfor,
                                    ref_group_names=c("iso") )   #选基础的细胞  或者样本 看meta的输入类型   也可以不选择 算法根据平均值来自己算



infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff= 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="cnv-results_TFCAFC/",  #每次重算记得改位置
                             cluster_by_groups=TRUE,  # 选择TRUE是按样本分组 改为FALSE会进行按另一个参数k_obs_groups给出的分组数（默认为1）进行分组
                             denoise=F,     #是否去噪
                             num_threads = 16, #使用核core
                             HMM=T)   # 是否基于HMM预测CNV,True的话时间很久
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff= 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="cnv-results_TFCAFC2/",  #每次重算记得改位置
                             cluster_by_groups=TRUE,  # 选择TRUE是按样本分组 改为FALSE会进行按另一个参数k_obs_groups给出的分组数（默认为1）进行分组
                             denoise=T,     #是否去噪
                             num_threads = 16, #使用核core
                             HMM=T)   # 是否基于HMM预测CNV,True的话时间很久

#结果输出在当前工作路径下的out_dir文件夹下（最终会输出很多文件在out_dir的目录下）   可以直接用里面的热图
summary(scRNA_harmony$nFeature_Spatial)



#第一步：提取inferCNV的结果  
meta <- rownames_to_column(meta)
obs <- read.table("cnv-results_TFCAFC/infercnv.observations.txt",header = T,check.names = F)     #首先这个obs是每一个基因在每一个细胞的拷贝信息,相当于该基因的拷贝量
dim(obs)
nrow(obs)
obs[1:3,1:4]

genelist <- c('BRAF','KRAS','NRAS','HRAS','RAF1','AKT1','AKT2','AKT3','PTEN','PIK3CA','PAX8' ,'PPARG','CTNNB1','RET','TERT') 
genelist = c('BRAF','KRAS','NRAS','HRAS',
             'AKT1','AKT2','AKT3','PTEN',
             'PIK3CA','PAX8' ,'PPARG')
as.numeric(obs['PTEN',])
score = data.frame(row.names = colnames(obs),
                   score = as.numeric(obs['PTEN',])) #针对单个基因的

score=as.data.frame(colSums(obs[which(rownames(obs)%in%genelist),]))   #把obs的每一个基因拷贝量加起来，就是这个细胞的总拷贝数obs
colnames(score) <- 'score'

#提取meta信息  celltype,seurat_clusters,orig.ident
head(scRNA_harmony@meta.data)
meta <- subset(scRNA_harmony2@meta.data,select = c("celltype"))  #提取该细胞的其他的meta信息

#将meta信息添加给score
head(meta)
head(score)
meta <- rownames_to_column(meta)
score <- rownames_to_column(score)
score <- merge(score,meta,by.x = "rowname",by.y = "rowname")    #这里会可能损失一些细胞   为什么呢，因为在前面infer时，有一些细胞无法推断，会被删掉，但是总体上问题不大

head(score)
table(score$celltype)

ggboxplot(score,x="celltype",y="score",fill="celltype")


#投影到seurat中
dim(obs)
dim(scRNA_harmony)
scRNA_harmony1 <- subset(scRNA_harmony , cells = colnames(obs))
scRNA_harmony1[['CNV']] <- CreateAssayObject(data = as.matrix(obs))

DefaultAssay(scRNA_harmony1) <- 'CNV'
scRNA_harmony1 <- ScaleData(scRNA_harmony1,features = rownames(scRNA_harmony1))
levels(scRNA_harmony1) <- c("FC1","AFC","T")
VlnPlot(scRNA_harmony1,features = c('BRAF','KRAS','NRAS','HRAS','RAF1','AKT1','AKT2','AKT3','PTEN','PIK3CA','PAX8' ,'PPARG','CTNNB1','RET','TERT')  )

# #取重点关注的亚群绘图
# scRNA_harmony2 <- subset(scRNA_harmony1,idents=c("T","AFC","FC1"))
# levels(scRNA_harmony2) <- c("FC1","AFC","T")
# VlnPlot(scRNA_harmony2,features = c('BRAF','KRAS','NRAS','HRAS',
#                                     'AKT1','AKT2','AKT3','PTEN',
#                                     'PIK3CA','PAX8' ,'PPARG'),
#         pt.size = 0,stack = T)
# VlnPlot(scRNA_harmony2,features = c('BRAF','KRAS','NRAS','HRAS',
#                                     'RAF1','AKT1','AKT2','AKT3','PTEN',
#                                     'PIK3CA','PAX8' ,'PPARG','CTNNB1','RET'),
#         pt.size = 0,cols = c("#6DA1F2","#D7EBF8", "#FBE2DA"))
# 
# stacked_violin_plot(
#   gene= c('BRAF','KRAS','NRAS','HRAS', 'RAF1','AKT1','AKT2','AKT3','PTEN',
#          'PIK3CA','PAX8' ,'PPARG','CTNNB1','RET','TERT'),
#   scRNA_harmony2,  cluster = NULL,  limits.max = 2,  width = 13,  height = 10.3,
#   flip = T,  filename = "",  text.size = 10,  Mean = T,  col = NULL)


# heatmap 跟常规的差太远
if(F){
  #第一步
  #我们先解决细胞注释的表格
  #由于我们要定义CNV的程度，所以首先要将数值变成注释
  #加change列,标记inferCNV的变化
  sel_gene_cnv <- obs[which(rownames(obs)%in%genelist),]
  score=as.data.frame(colSums(sel_gene_cnv))   #把obs的每一个基因拷贝量加起来，就是这个细胞的总拷贝数obs
  colnames(score) <- 'scoress'
  
  max(score$scoress)          #根据最大最小值来定义
  min(score$scoress)          #一般变化倍数大于2倍的可以定义为拷贝数变异
  
  score$mutate <-  ifelse(score$scoress < 14,"normal",ifelse(score$scoress > 14.5,"High","Low"))  #要看懂这段代码  ifelse的条件句  假如符合K1标准为normal，不符合就按照ifelse(k2,"High","Low")这个条件再分，假如符合K2标准就High，不符合就Low
  table(score$mutate)   #这样就做好了一个anno
  
  score <- score[order(score$scoress,decreasing = T),]   #对数据框进行升序排列  如果加decreasing = T就实现降序排列,这样就会对后面热图的细胞顺序进行排序，为什么想对细胞排序呢，因为我想让CNV属于normal，High，Low的细胞排一起，
  
  score <- score[,-c(1)]  #删除多余的
  
  #第二步，解决细胞每一个基因的表达量
  #由于前面score跟meta merge了，所以这里要对dat进行匹配
  dat=dat[,match(rownames(score), colnames(dat))]
  identical(rownames(score), colnames(dat))#这样每一个细胞的anno和它的基因表达量就对应了，记住一定要对应，要想明白哈，有些东西虽然不报错，但是你对应错了，结果也是错的。
  dat <- log2(dat+1)#这一步的目的是给表达矩阵取log2，增加基因之间的差异
  
  #第三步 画图啦
  #想自定义annotation的颜色       
  ann_colors = list(
    change = c(normal="black", Low="white",High="red"))   #这里也可以改变多个anno的颜色  
  
  #画热图
  library(pheatmap)
  pheatmap(dat,
           show_colnames =F,
           show_rownames = F,
           scale = "row",
           cluster_cols = F,                                        #这些参数建议去Google
           cluster_rows = F, 
           annotation_col=score,
           annotation_colors = ann_colors,
           colorRampPalette(c("#00FF00", "white", "#DC143C"))(75),   #改热图的颜色
           breaks = seq(-3,3,length.out = 100))
  
}