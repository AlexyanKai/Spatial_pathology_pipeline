######Load libraries###########
library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(nichenetr)
library(tidyr)
library(tibble)
###############################

######Set work dir#############
setwd('F:/A--Lianlabwork/A020--YkSpt20211228/Step01_RegionNichenet')
###############################

######Read data################
# a.object<-readRDS('F:/A--Lianlabwork/A020--YkSpt20211228/Step00_YkAnalyzed/Step00_SampleA.rds')
DimPlot(a.object_rename2)
table(a.object_rename2$TnT_cluster20211031)
a.object_rename2@active.ident <- as.factor(a.object_rename2$TnT_cluster20211031)
SpatialDimPlot(a.object_rename2,label=TRUE)
a.object_rename2_Tsubset <- subset(a.object_rename2,idents=c("AFC","AFCI","FC1", "FC1I","FC2", "FC2I", "I", "nTnFCnI","T1","T2", "T3","T4","T5","T6","T7","T8"))
table(a.object_rename2_Tsubset@active.ident)
SpatialDimPlot(a.object_rename2_Tsubset,label=TRUE)
ls(a.object_rename2_Tsubset@assays)
#Only spatial
# a.object_rename2_Tsubset <- SCTransform(a.object_rename2_Tsubset,assay="Spatial",variable.features.n = 30000) #本来就有，此步骤在这个里面不需要做
#dim(a.object_rename2_Tsubset@assays$SCT@scale.data)

###############################

######Load Nichenet Database#########
ligand_target_matrix=readRDS('NichenetResult20211229/NicheNetDbs/ligand_target_matrix.rds')
lr_network=readRDS('NichenetResult20211229/NicheNetDbs/lr_network.rds')
weighted_networks=readRDS('NichenetResult20211229/NicheNetDbs/weighted_networks.rds')
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
#####################################

###!!! We are using Nichenet for spot/region analysis. !!!###


######Analysis-1#####################
#Set receiver(being regulated by ligands from others) identity
receiver='T2'
#Find expressed genes and background genes in receiver
expressed_genes_receiver = get_expressed_genes(receiver, a.object, pct = 0.5)
# The expressed genes are defined based on the SCT slot, but this can be changed via the assay_oi parameter.
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
#Find expressed genes in sender
sender_celltypes=c('AFCI')
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, a.object, 0.5) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
#Find differential expressed genes between T2 and T4
condition_oi='T2'
condition_reference='T4'
DE_table_receiver = FindMarkers(object = a.object,
                                ident.1 = condition_oi,
                                ident.2 = condition_reference,
                                min.pct = 0.50) %>% rownames_to_column("gene")
#Choose significant genes and they should be in database
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
#Define  potential ligands, ligands expressed in sender cell population
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
#Run prediction
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
#Choose best
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
#Make dot plot
#DotPlot(a.object, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
#Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
#p_ligand_target_network
#Receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
#p_ligand_receptor_network
#Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
#p_ligand_receptor_network_strict

######Try integrated function######
source('NichenetResult20211229/nichenet4spacial.R')
nichenet4spatial(a.object_rename2_Tsubset,'T2','T4',c('Imm1'),'Imm1_T2-T4',0.5,0.5,0.5,0.25,20)
nichenet4spatial(a.object_rename2_Tsubset,'T2','T4',c('Imm2'),'Imm2_T2-T4',0.5,0.5,0.5,0.25,20)
nichenet4spatial(a.object_rename2_Tsubset,'T2','T4',c('Imm4'),'Imm4_T2-T4',0.5,0.5,0.5,0.25,20)
nichenet4spatial(a.object_rename2_Tsubset,'T2','T4',c('Imm5'),'Imm5_T2-T4',0.5,0.5,0.5,0.25,20)
nichenet4spatial(a.object_rename2_Tsubset,'T2','T4',c('Imm6'),'Imm6_T2-T4',0.5,0.5,0.5,0.25,20)
#20110111
nichenet4spatial(a.object_rename2_Tsubset,'T2'#目的群
                 ,'T4' #参照群
                 ,c('AFC') #受到谁的影响，配体在的群
                 ,'AFC_T2-T4',0.5,0.5,0.5,0.25,20)
nichenet4spatial(a.object_rename2_Tsubset,'T2','T4',c('AFCI'),'AFCI_T2-T4',0.5,0.5,0.5,0.25,20)
nichenet4spatial(a.object_rename2_Tsubset,'T2','T4',c('AFC','AFCI'),'AFCandAFCI_T2-T4',0.5,0.5,0.5,0.25,20)
nichenet4spatial(a.object_rename2_Tsubset,'T2','T4',c('FC1'),'FC1_T2-T4',0.5,0.5,0.5,0.25,20)
nichenet4spatial(a.object_rename2_Tsubset,'T2','T4',c('FC1I'),'FC1I_T2-T4',0.5,0.5,0.5,0.25,20)
nichenet4spatial(a.object_rename2_Tsubset,'T2','T4',c('FC1','FC1I'),'FC1andFC1I_T2-T4',0.5,0.5,0.5,0.25,20)

