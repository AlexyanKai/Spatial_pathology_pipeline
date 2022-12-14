Get_percentplot_split3 <-
  function (seurat_object, filename = NULL, order = "increasing", 
          split_by = NULL, colorm = NULL) {
  variable <- count <- value <- Sample <- Percent <- Cluster <- NULL
  num <- table(seurat_object@active.ident)
  if (is.null(split_by)) {
    split_by <- "samples"
    sample_name = unique(seurat_object@meta.data$samples)
  }
  else {
    sample_name = unique(seurat_object@meta.data[, split_by])
  }
  frequency_matrix <- data.frame(num)
  colnames(frequency_matrix) = c("type", "num")
  for (i in sample_name) {
    cell_name = row.names(seurat_object@meta.data[seurat_object@meta.data[, 
                                                                          split_by] == i, ])
    subset_data = subset(seurat_object, cells = cell_name)
    frequency = table(subset_data@active.ident)/ncol(subset_data)
    frequency1 = data.frame(frequency)
    colnames(frequency1) = c("type", i)
    frequency_matrix = merge(frequency_matrix, frequency1, 
                             by = "type", all = T)
  }
  frequency_matrix[is.na(frequency_matrix)] <- 0
  o = order(frequency_matrix$type, decreasing = F)
  frequency_matrix = frequency_matrix[o, ]
  rownames(frequency_matrix) = frequency_matrix$type
  colnames(frequency_matrix)[1] = "cluster"
  frequency_matrix = frequency_matrix[, -2]
  frequency_matrix[] <- lapply(frequency_matrix, as.character)
  frequency_matrix[is.na(frequency_matrix)] <- 0
  frequency_matrix[, 2:ncol(frequency_matrix)] <- as.data.frame(lapply(frequency_matrix[, 
                                                                                        2:ncol(frequency_matrix)], as.numeric))
  frequency_matrix$cluster = as.character(frequency_matrix$cluster)
  if (order == "increasing") {
    o = order(colnames(frequency_matrix), decreasing = F)
  }
  else {
    o = order(colnames(frequency_matrix), decreasing = T)
  }
  frequency_matrix = frequency_matrix[, o]
  write.csv(frequency_matrix, paste0(filename, "_frequency_matrix.csv"), 
            row.names = T)
  data_rownames <- rownames(frequency_matrix)
  frequency_matrix_m <- melt(frequency_matrix, id.vars = c("cluster"))
  frequency_matrix_m <- frequency_matrix_m %>% group_by(variable) %>% 
    mutate(count = sum(value)) %>% mutate(freq = round(100 * 
                                                         value/count, 2))
  colnames(frequency_matrix_m) = c("Cluster", "Sample", 
                                   "Frequency", "count", "Percent")
  if (is.null(colorm)) {
    colorful <- LIANLAB::colorful
    colors 
    colorm <- data.frame(frequency_matrix$cluster, colors[1:length(rownames(frequency_matrix))])
    row.names(colorm) <- colorm[, 1]
  }
  p <- ggplot(frequency_matrix_m, aes(x = Sample, y = Percent,
                                      group = factor(Cluster, levels=c("tumor",  "ant",    "fibro", "FC2", "immune", "other")))) + geom_bar(stat = "identity",
                                      position = "fill", aes(fill = Cluster)) + scale_fill_manual(values = colorm[,
                                      2], breaks = colorm[, 1]) + theme_bw() + theme(panel.border = element_blank(),
                                      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      axis.line = element_line(colour = "black"), axis.text.x = element_text(vjust = 0.5,
                                      hjust = 0.5, angle = 45))
     pdf(paste0("percent_plot_", filename, ".pdf"), 
      width = ncol(frequency_matrix) * 0.8 + 2, height = nrow(frequency_matrix) * 
        0.25 + 1)
  print(p)
  dev.off()
  }

