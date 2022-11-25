######Load libraries###########
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(nichenetr)
library(tidyr)
library(tibble)
###############################
nichenet4spatial<-function(srtobj,
                           receiver,
                           reference,
                           sender_celltypes,
                           figure_dir,
                           receiver_expression_threshold,
                           sender_expression_threshold,
                           diff_exp_exp_threshold,
                           diff_exp_logfc_threshold,
                           ligand_topn){
  dir.create(figure_dir)
  #Set receiver(being regulated by ligands from others) identity
  #receiver='T2'
  #Find expressed genes and background genes in receiver
  expressed_genes_receiver = get_expressed_genes(receiver, srtobj, pct = receiver_expression_threshold)
  # The expressed genes are defined based on the SCT slot, but this can be changed via the assay_oi parameter.
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  #Find expressed genes in sender
  #sender_celltypes=c('Imm1')
  list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, srtobj, sender_expression_threshold) # lapply to get the expressed genes of every sender cell type separately here
  expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
  #Find differential expressed genes between T2 and T4
  #condition_oi='T2'
  #condition_reference='T4'
  DE_table_receiver = FindMarkers(object = srtobj, 
                                  ident.1 = receiver, 
                                  ident.2 = reference, 
                                  min.pct = diff_exp_exp_threshold) %>% rownames_to_column("gene")
  #Choose significant genes and they should be in database
  geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= diff_exp_logfc_threshold) %>% pull(gene)
  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  #Define  potential ligands, ligands expressed in sender cell population
  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  #Run prediction
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))
  #Choose best
  best_upstream_ligands = ligand_activities %>% top_n(ligand_topn, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
  #Make dot plot
  png(paste0(figure_dir,'/Figure1_LigandActivityAmongRegions.png'),
      width=200+30*ligand_topn,
      height=100+30*length(levels(srtobj@active.ident)))
  print(DotPlot(srtobj, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis())
  dev.off()
  #Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
  png(paste0(figure_dir,'/Figure2_Ligand2Target.png'),
      width=100+20*length(levels(p_ligand_target_network$data$x)),
      height=200+20*length(levels(p_ligand_target_network$data$y)))
  print(p_ligand_target_network)
  dev.off()
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
  png(paste0(figure_dir,'/Figure3_Ligand2Receptor.png'),
      width=100+20*length(levels(p_ligand_receptor_network$data$x)),
      height=200+20*length(levels(p_ligand_receptor_network$data$y)))
  print(p_ligand_receptor_network)
  dev.off()
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
  png(paste0(figure_dir,'/Figure4_Ligand2ReceptorPPI.png'),
      width=100+20*length(levels(p_ligand_receptor_network_strict$data$x)),
      height=200+20*length(levels(p_ligand_receptor_network_strict$data$y)))
  print(p_ligand_receptor_network_strict)
  dev.off()
}