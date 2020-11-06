enriched_keys <- function(file = F, path, enrichment_keys = NULL, cluster_info){
  output <- list()
  clusters <- unique(cluster_info$cluster)
  
  # if file == F, the enrichment_keys must be provided and are checked against GO enrichment terms
  if(!file){
    if(is.null(enrichment_keys)){
      print("file is set to FALSE but no enrichment keys were provided")
      return()
    }
    enriched_list <- list()
    cell_types_per_cluster <- NULL
    for(c in clusters){
      
      genes <- cluster_info$gene
      
      enrich <- clusterProfiler::enrichGO(genes,
                                          OrgDb = "org.Hs.eg.db",
                                          keyType = "SYMBOL",
                                          ont = "BP",
                                          pvalueCutoff = 0.05)
      enriched_list[[c]] <- enrich
      
      cell_enrich <- list(counts = list(), genes = list())
      
      hits <- 0
      top_10 <- enrich@result[order(enrich@result$Count, decreasing = T),][1:20,]
      for(type in enrichment_keys){
        
        tmp <- dplyr::filter(top_10, grepl(type, top_10$Description, fixed = TRUE) == T) %>%
          dplyr::pull(., "geneID") %>%
          paste0(., collapse = "/") %>% 
          base::strsplit(., split = "/") %>% 
          unlist(.) %>% 
          unique(.)
        hits <- hits + length(tmp)
        cell_enrich$counts[[type]] <- length(tmp)
        cell_enrich$genes[[type]] <- tmp
        
      }
      tmp <- data.frame(matrix(unlist(cell_enrich$counts), ncol= length(cell_enrich$counts), byrow=T) %>% t(),stringsAsFactors=FALSE)%>%
        cbind(., names(cell_enrich$counts))%>%
        cbind(., rep(c, length(cell_enrich$counts)))
      
      sink(file = paste0(working_directory, "/enrichedInKeys_",c,".txt"))
      for (n in names(cell_enrich$genes)){
        cat(NULL, sep = "\n")
        cat(n, sep = "\n\n")
        cat(cell_enrich$genes[[n]], sep = "\n")
      }
      
      sink()
      colnames(tmp) <- c("count", "cell_type", "cluster")
      tmp$hits <- rep(hits, nrow(tmp))
      
      
      cell_types_per_cluster <- rbind(cell_types_per_cluster, tmp)
    }
    output[["enrichlist"]] <- enriched_list
  }else{
    f <- read.csv(path)
    enrichment_keys <- colnames(f)
    cell_types_per_cluster <- NULL
    
    for(c in clusters){
      print(c)
      genes <- cluster_info[cluster_info$cluster == c,]%>%
        dplyr::pull(., "gene")
      
      cell_enrich <- list(counts = list(), genes = list())
      hits <- 0
      
      for (type in enrichment_keys){
        tmp <- genes[genes %in% f[,c(type)]]
        print(length(tmp))
        hits <- hits + length(tmp)
        cell_enrich$counts[[type]] <- length(tmp)
        cell_enrich$genes[[type]] <- tmp
      }
      
      tmp <- data.frame(matrix(unlist(cell_enrich$counts), ncol= length(cell_enrich$counts), byrow=T) %>% t(),stringsAsFactors=FALSE)%>%
        cbind(., names(cell_enrich$counts))%>%
        cbind(., rep(c, length(cell_enrich$counts)))
      
      sink(file = paste0(working_directory, "/enrichedInKeys_",c,".txt"))
      for (n in names(cell_enrich$genes)){
        cat(NULL, sep = "\n")
        cat(n, sep = "\n\n")
        cat(cell_enrich$genes[[n]], sep = "\n")
      }
      
      sink()
      colnames(tmp) <- c("count", "cell_type", "cluster")
      if(hits == 0){
        tmp$count <- 0
      }else{
        tmp$count <- (tmp$count/hits)*100
      }
      
      tmp$hits <- rep(hits, nrow(tmp))
      
      
      cell_types_per_cluster <- rbind(cell_types_per_cluster, tmp)
      
    }
  }
  output[["cell_types_per_cluster"]] <-cell_types_per_cluster
  return(output)
}




signature_comparison <- function(counts, anno, cluster_info, IDcol, voi, ctrl, 
                                 col_order = NULL, row_order = NULL){
  
  output <- list()
  
  #extract signature genes:
  genes <- cluster_info$gene
  
  
  # filter data for signature genes:
  counts_filt <- counts[rownames(counts) %in% genes,]
  

  # calculate GFCs across conditions:
  counts_filt[is.na(counts_filt)] <- 0
  
  
  tmp_GFCS <- compare_external_signature(sample_file = counts_filt, 
                                         anno_file =  anno,
                                         grpvar = voi,
                                         cluster_info = cluster_info)
  
  
  tmp_GFCS <- tmp_GFCS[rownames(tmp_GFCS) %in% rownames(counts_filt),]
  colnames(tmp_GFCS) <- c(as.character(unique(dplyr::pull(anno, voi))), "Gene")
  
  
  # without controls:
  
  output[["GFCs_no_ctrls"]] <- tmp_GFCS[, colnames(tmp_GFCS)[grepl(ctrl, colnames(tmp_GFCS))==F]]
  
  
  return(output)
}





make_rownames_unique <- function(counts_new){
  
  counts_new <- counts_new[!duplicated(counts_new$SYMBOL),] %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(., "SYMBOL")
  
  # remove all non-numeric columns (description, gene id, etc.):
  for (x in colnames(counts_new)){
    if(!is.numeric(counts_new[[x]])){
      counts_new[[x]] <- NULL
    }
  }
  if(ncol(counts_new) == 0){
    print("All columns were deleted when removing non-numeric columns. 
          Please check the data type of your expression values.")
  }
  
  return(counts_new)
}





calculate_GFCs <- function(expressions, anno, genes, grpvar){
  
# filter expression data for genes of which the GFC shall be calculated
  expressions <- tibble::rownames_to_column(expressions, var = "SYMBOL")%>%
    dplyr::filter(., SYMBOL %in% genes)%>%
    tibble::column_to_rownames(., "SYMBOL")

# calculate GFCs for groups
  if(!grpvar %in% colnames(anno)){
    print("The column you specified to contain the grouping variables does not exists in the annotation table.")
    return()
  }
  
  GFCs <- NULL
  for (g in genes){
# extract expression values
    g_df <- expressions[rownames(expressions) == g,]
    colnames(g_df) <- colnames(expressions)
    
# vector of group means:
    mean_vec <- NULL
    for (var in unique(anno[[grpvar]])){
      var_anno <- anno[anno[[grpvar]] == var,]
      var_exp <- g_df[, colnames(g_df) %in% var_anno$ID]
      if(nrow(var_anno) == 1){
        mean_vec <- c(mean_vec,mean(var_exp %>% as.numeric()))
      }else{
        mean_vec <- c(mean_vec,mean(var_exp[1,] %>% as.numeric()))
      }
      
      
    }
    
# overall expression mean of current gene:
    g_mean <- mean(mean_vec)
    GFC_vec <- NULL
    for (x in 1:length(mean_vec)){
      GFC_vec <- c(GFC_vec, gtools::foldchange(mean_vec[x], g_mean) %>% 
                     ifelse(.>2, 2,.) %>%
                     ifelse(.< (-2), -2,.))
    }
    
# add this genes GFC vector to GFCs
    GFCs <- rbind( GFCs, GFC_vec)
  }
  
  GFCs <- as.data.frame(GFCs) %>%
    round(., digits = 3)
  
# add column with gene names
  GFCs$Gene <- genes
  colnames(GFCs) <- c(unique(anno[[grpvar]]), "Gene") 
  rownames(GFCs) <- GFCs$"Gene"
  return(GFCs)
}



cluster_partition_and_plot <- function(cluster_info, GFCs, col_order, cell_types_per_cluster1, 
                                       cell_types_per_cluster2, row_order){
  
# filter for included clusters (non-white)
  mat_heatmap <- NULL
  
  
  if(!is.null(row_order)){
    for (c in row_order){
      
#get genes from the original cluster
      genes <- cluster_info[cluster_info$cluster == c,] %>%
        dplyr::pull(., "gene")
      
# GFCs of new data set, where genes are found in original cluster
      c_GFCs <- dplyr::filter(GFCs, Gene %in% genes)
      c_GFC_means <- apply(c_GFCs[, c(1:(ncol(c_GFCs)-1))], 2, mean)
      
      mat_heatmap <- rbind(mat_heatmap, c_GFC_means)
      
    }
    rownames(mat_heatmap) <- row_order
    
  }else{
    for (c in unique(cluster_info$cluster)){
      
#get genes from the original cluster
      genes <- cluster_info[cluster_info$cluster == c,] %>%
        dplyr::pull(., "gene")
      
      
      
# GFCs of new data set, where genes are found in original cluster
      c_GFCs <- dplyr::filter(GFCs, Gene %in% genes)
      c_GFC_means <- apply(c_GFCs[, c(1:(ncol(c_GFCs)-1))], 2, mean)
      
      mat_heatmap <- rbind(mat_heatmap, c_GFC_means)
      
    }
    rownames(mat_heatmap) <- unique(cluster_info$cluster)
  }
  
  
  colnames(mat_heatmap) <- colnames(GFCs)[1:(ncol(GFCs)-1)]
  
  if(!is.null(col_order)){
    mat_heatmap <- mat_heatmap %>% as.data.frame()
    mat_heatmap <- mat_heatmap[, col_order] %>% as.matrix()
  }
  
  
# creating the data for the barplots in the heatmap annotation:
  enrich_mat1 <- list()
  enrich_count1 <- list()
  enrich_mat2 <- list()
  enrich_count2 <- list()
  
  if(!is.null(row_order)){
    if(!is.null(cell_types_per_cluster1)){
      for(x in row_order){
        enrich_mat1[[x]] <- dplyr::filter(cell_types_per_cluster1, cluster == x)%>%
          dplyr::pull(., count)
        enrich_count1[[x]] <- dplyr::filter(cell_types_per_cluster1, cluster == x)%>%
          dplyr::pull(., hits)%>%
          dplyr::first(.)
      }
    }
    if(!is.null(cell_types_per_cluster2)){
      for(x in row_order){
        enrich_mat2[[x]] <- dplyr::filter(cell_types_per_cluster2, cluster == x)%>%
          dplyr::pull(., count)
        enrich_count2[[x]] <- dplyr::filter(cell_types_per_cluster2, cluster == x)%>%
          dplyr::pull(., hits)%>%
          dplyr::first(.)
      }
    }
    
  }else{
    for(x in unique(cell_types_per_cluster1$cluster)){
      enrich_mat1[[x]] <- dplyr::filter(cell_types_per_cluster1, cluster == x)%>%
        dplyr::pull(., count)
      enrich_count1[[x]] <- dplyr::filter(cell_types_per_cluster1, cluster == x)%>%
        dplyr::pull(., hits)%>%
        dplyr::first(.)
    }
    for(x in unique(cell_types_per_cluster2$cluster)){
      enrich_mat2[[x]] <- dplyr::filter(cell_types_per_cluster2, cluster == x)%>%
        dplyr::pull(., count)
      enrich_count2[[x]] <- dplyr::filter(cell_types_per_cluster2, cluster == x)%>%
        dplyr::pull(., hits)%>%
        dplyr::first(.)
    }
  }
  
  if(!length(enrich_mat1) == 0){
    enrich_mat1 <- matrix(unlist(enrich_mat1), nrow =length(enrich_mat1), byrow = T)
    enrich_count1 <- unlist(enrich_count1)
    
  }
  if(all(enrich_mat1 == 0) == T){
    enrich_mat1 <- list()
  }
  
  if(!length(enrich_mat2) == 0){
    
    enrich_mat2 <- matrix(unlist(enrich_mat2), nrow =length(enrich_mat2), byrow = T)
    enrich_count2 <- unlist(enrich_count2)
  }
  if(all(enrich_mat2 == 0) == T){
    enrich_mat2 <- list()
  }
  
  c_df <- NULL
  for(c in unique(cluster_info$cluster)){
    c_df <- rbind(c_df, data.frame(color = c, gene_no = length(cluster_info$cluster[cluster_info$cluster == c])))
  }
# setting the order of rows and creating the respective colour vector
  if(!is.null(row_order)){
    cluster_colors <- factor(row_order)
    names(cluster_colors) <- row_order
    c_df <- c_df[match(rowvec, c_df$color),]
  }else{
    cluster_colors <- factor(c_df$color)
    names(cluster_colors) <- c_df$color
    row_order <- unique(c_df$color)
  }
  
  
  if(length(enrich_mat1) == 0 & length(enrich_mat2) == 0){
    ha <- ComplexHeatmap::HeatmapAnnotation(clusters = anno_simple(row_order, col = cluster_colors, 
                                                                   simple_anno_size = unit(0.5, "cm")),
                                            genes = anno_barplot(c_df$gene_no, width = unit(2.5, "cm")),
                                            
                                            
                                            which = "row", 
                                            width = unit(4.5, "cm"),
                                            annotation_name_side = "top",
                                            gap = unit(2, "mm"), 
                                            annotation_name_rot = 0,
                                            annotation_name_gp = gpar(fontsize = 8))
    
    lgd_list <- list(
      
    )
  }
  else if(length(enrich_mat1) > 0 & length(enrich_mat2) == 0){
    ha <- ComplexHeatmap::HeatmapAnnotation(clusters = anno_simple(row_order, col = cluster_colors, 
                                                                   simple_anno_size = unit(0.5, "cm")),
                                            genes = anno_barplot(c_df$gene_no, width = unit(2.5, "cm")),
                                            
                                            enriched_count = anno_text(paste0(enrich_count1, "/", c_df$gene_no), width = unit(1.5, "cm")),
                                            enriched = anno_barplot(enrich_mat1,
                                                                    width = unit(3, "cm"),
                                                                    gp = gpar(fill = ggsci::pal_d3(palette = "category20")(20),
                                                                              col = ggsci::pal_d3(palette = "category20")(20))),
                                            which = "row", 
                                            width = unit(9, "cm"),
                                            annotation_name_side = "top",
                                            gap = unit(2, "mm"), 
                                            annotation_name_rot = 0,
                                            annotation_name_gp = gpar(fontsize = 8))
    
    lgd_list <- list(
      
      ComplexHeatmap::Legend(labels = unique(cell_types_per_cluster1$cell_type), title = "enriched",
                             legend_gp = gpar(col = ggsci::pal_d3(palette = "category20")(20)),
                             type = "points", pch = 15)
    )
  }
  else if(length(enrich_mat1) == 0 & length(enrich_mat2) > 0){
    ha <- ComplexHeatmap::HeatmapAnnotation(clusters = anno_simple(row_order, col = cluster_colors, 
                                                                   simple_anno_size = unit(0.5, "cm")),
                                            genes = anno_barplot(c_df$gene_no, width = unit(1.5, "cm")),
                                            
                                            enriched_count = anno_text(paste0(enrich_count2, "/", c_df$gene_no), width = unit(1.5, "cm"),
                                                                       gp = gpar(fontsize = 8)),
                                            enriched = anno_barplot(enrich_mat2,
                                                                    width = unit(5, "cm"),
                                                                    gp = gpar(fill = ggsci::pal_d3(palette = "category20")(20),
                                                                              col =ggsci::pal_d3(palette = "category20")(20))),
                                            
                                            which = "row", 
                                            width = unit(9, "cm"),
                                            annotation_name_side = "top",
                                            gap = unit(2, "mm"), 
                                            annotation_name_rot = 0,
                                            annotation_name_gp = gpar(fontsize = 8))
    
    lgd_list <- list(
      
      
      ComplexHeatmap::Legend(labels = unique(cell_types_per_cluster2$cell_type), title = "enriched",
                             legend_gp = gpar(col = ggsci::pal_d3(palette = "category20")(20)),
                             type = "points", pch = 15)
    )
  }else{
    ha <- ComplexHeatmap::HeatmapAnnotation(clusters = anno_simple(row_order, col = cluster_colors, 
                                                                   simple_anno_size = unit(0.25, "cm")),
                                            genes = anno_barplot(c_df$gene_no, width = unit(0.75, "cm")),
                                            
                                            enriched_count_1 = anno_text(paste0(enrich_count1, "/", c_df$gene_no),
                                                                         width = unit(0.75, "cm"),
                                                                         gp = gpar(fontsize = 8)),
                                            enriched_1 = anno_barplot(enrich_mat1,
                                                                      width = unit(2, "cm"),
                                                                      gp = gpar(fill = ggsci::pal_d3(palette = "category20")(20),
                                                                                col = ggsci::pal_d3(palette = "category20")(20)),
                                                                      baseline = 0),
                                            enriched_count_2 = anno_text(paste0(enrich_count2, "/", c_df$gene_no), 
                                                                         width = unit(0.75, "cm"),
                                                                         gp = gpar(fontsize = 8)),
                                            enriched_2 = anno_barplot(enrich_mat2,
                                                                      width = unit(2, "cm"),
                                                                      gp = gpar(fill = ggsci::pal_d3(palette = "category20")(20),
                                                                                col = ggsci::pal_d3(palette = "category20")(20)),
                                                                      baseline = 0),
                                            which = "row", 
                                            width = unit(12, "cm"),
                                            annotation_name_side = "top",
                                            gap = unit(2, "mm"), 
                                            annotation_name_rot = 0,
                                            annotation_name_gp = gpar(fontsize = 8))
    
    lgd_list <- list(
      
      
      ComplexHeatmap::Legend(labels = unique(cell_types_per_cluster1$cell_type), title = "enriched_1",
                             legend_gp = gpar(col = ggsci::pal_d3(palette = "category20")(20)),
                             type = "points", pch = 15),
      ComplexHeatmap::Legend(labels = unique(cell_types_per_cluster2$cell_type), title = "enriched_2",
                             legend_gp = gpar(col = ggsci::pal_d3(palette = "category20")(20)),
                             type = "points", pch = 15)
    )
  }
  
  Cairo(file = paste0(working_directory, "/Heatmap_external_sig.pdf"),
        width = 50, 
        height = 30,
        pointsize=11,
        dpi=300,
        type = "pdf",
        units = "in")
  
  
  hm <- ComplexHeatmap::Heatmap(mat_heatmap, 
                                right_annotation = ha,
                                col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(length(seq(-2, 2, by = .1))),
                                clustering_distance_rows = "euclidean",
                                clustering_distance_columns = "euclidean",
                                clustering_method_rows = "complete",
                                clustering_method_columns = "complete",
                                cluster_columns = F,
                                cluster_rows = F,
                                column_names_rot = 90,
                                column_names_centered = F,
                                row_names_gp = gpar(fontsize = 8),
                                column_names_gp = gpar(fontsize = 8),
                                rect_gp = grid::gpar(col = "grey"),
                                heatmap_legend_param = list(title = "", legend_height = unit(3, "cm")))
  
  hm_w_lgd <- ComplexHeatmap::draw(hm, annotation_legend_list = lgd_list, merge_legends = T, 
                                   padding = unit(c(2, 2, 2, 30), "mm"))
  
  dev.off()
  
  print(hm_w_lgd)
  return(hm_w_lgd)
  
}

compare_external_signature <- function(sample_file,
                                       anno_file,
                                       grpvar,
                                       cluster_info){
  
  
  tmp_GFCS <- calculate_GFCs(sample_file, anno_file, genes =  cluster_info$gene, grpvar = grpvar)
  return(tmp_GFCS)
  
}

# creating column vector for heatmap:

colvec <- c("Covid-19 Cluster 2", "Covid-19 Cluster 4", "Covid-19 Cluster 6", "Covid-19 Cluster 1", "Covid-19 Cluster 5", "Covid-19 Cluster 3",
            "Chikungunya acute", "Zika Early Acute", "Zika Late Acute", "Ebola vaccine",
            "HIV", "Tuberculosis HIV", "Tuberculosis HIV ART", 
            "Tuberculosis Lung", "Tuberculosis MTP", "Tuberculosis", "Latent tuberculosis", "uncomplicated sepsis", "severe sepsis", "septic shock", 
            "sepsis death", "SIRS",  "NLRC4-MAS anakinra treatment", "NLRC4-MAS No treatment", "NOMID anakinra treatment", 
            "NOMID No treatment", "Rheumatoid arthritis Unknown", "Autoimmune disease")

# creating row vector for heatmap:

rowvec <- c("indianred", "maroon", "darkorange", "steelblue", "gold", "lightgreen", "darkgreen", "orchid",
            "pink", "darkgrey")

working_directory <- "SET WORKING DIRECTORY"
cluster_information <- read.csv("INSERT PATH TO: cluster_to_gene.csv")

cell_types_per_cluster <- enriched_keys(file =  T, path ="INSERT PATH TO: immune_sig_m.csv",
                                        cluster_info = cluster_information)

cell_types_per_cluster2 <- enriched_keys(file =  T, path = "INSERT PATH TO: selected_neutro_markers_pos.csv",
                                         cluster_info = cluster_information)

ext_data <- read.csv("INSERT PATH TO: norm_counts_covgenes_noRP_noMT_noMEN_noFAT.txt", sep = "\t", row.names = 1, 
                     check.names = F)

ext_anno <- read_delim("INSERT PATH TO: norm_anno_covgenes_noRP_noMT_noMEN_noFAT.txt", 
                       "\t", escape_double = FALSE, trim_ws = TRUE)

compared_signatures <- signature_comparison(counts = ext_data, 
                                        anno = ext_anno, 
                                        cluster_info = cluster_information, 
                                        IDcol = "ID",
                                        voi = "merged",
                                        ctrl = "control",
                                        col_order = colvec,
                                        row_order = rowvec
)

hm <- cluster_partition_and_plot(cluster_info = cluster_information,
                                 
                                 GFCs = compared_signatures$GFCs_no_ctrls,
                                 col_order = colvec,
                                 cell_types_per_cluster1 = cell_types_per_cluster$cell_types_per_cluster,
                                 cell_types_per_cluster2 = cell_types_per_cluster2$cell_types_per_cluster,
                                 row_order = rowvec)

