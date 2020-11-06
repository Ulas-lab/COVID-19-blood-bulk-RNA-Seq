#create the lncRNA interaction data table
CreateTargetTable <- function(){
  

  all_genes <- unique(c(filt_cutoff_data$V1, filt_cutoff_data$V2))

  
  if(exists("lncRNA_DNA") == FALSE){
  lncRNA_DNA <- read.csv(paste0(working_directory, "data/", "RNAInter_interaction_RD.txt") ,
                         header = FALSE ,
                         check.names = FALSE,
                         sep = "\t",
                         stringsAsFactors = F)
  }
  colnames(lncRNA_DNA) <- c("RNAInterID", "Interactor1", "ID1", "Category1", "Species1",
                            "Interactor2", "ID2", "Category2", "Species2", "Score")
  lncRNA_DNA <- lncRNA_DNA[lncRNA_DNA$Species1 == "Homo sapiens" & 
                             lncRNA_DNA$Species2 == "Homo sapiens" & lncRNA_DNA$Category1 == "lncRNA", ]
  lncRNA_DNA <- lncRNA_DNA[,2:10]
  
  if(exists("lncRNA_RNA") == FALSE){
  lncRNA_RNA <- read.csv(paste0(working_directory, "data/", "RNA-RNA.txt") ,
                         header = FALSE ,
                         check.names = FALSE,
                         sep = "\t",
                         stringsAsFactors = F)
  }
  colnames(lncRNA_RNA) <- c("RNAInterID", "Interactor1", "ID1", "Category1", "Species1",
                            "Interactor2", "ID2", "Category2", "Species2", "Score")
  lncRNA_RNA <- lncRNA_RNA[lncRNA_RNA$Species1 == "Homo sapiens" & 
                             lncRNA_RNA$Species2 == "Homo sapiens" & lncRNA_RNA$Category1 == "lncRNA" &
                             lncRNA_RNA$Category2 == "mRNA", ]
  lncRNA_RNA <- lncRNA_RNA[,2:10]
  
  output <- merge(lncRNA_DNA, lncRNA_RNA, all = T)
  output <- output[output$Interactor2 %in% all_genes, ]
  tmp <- !duplicated(t(apply(output[c("Interactor1", "Interactor2")], 1, sort)))
  output <- output[tmp,]
  return(output)
}

lincr <- function(target_type = character(0), score){
  
  clusters_included <- cluster_information[cluster_information$cluster_included == "yes",]
  all_genes <- c(0)
  k <- 0
  for (i in 1:length(clusters_included$color)) {
    tmp <- unlist(strsplit(clusters_included$gene_n[i], ","))
    for (m in 1:length(tmp)) {
      k <- k + 1
      all_genes[k] <- tmp[m]
    }
  }

top_enrich <- list()
lnc_anno <- norm_anno[norm_anno$GENETYPE == "lncRNA",]
lnc_anno <- lnc_anno$SYMBOL
lncTarget_table <- lncTarget_table[lncTarget_table$Interactor1 %in% lnc_anno,]
lncTarget_table <- lncTarget_table[lncTarget_table$Interactor1 %in% c("HOTAIR", "HOTAIRM1","NEAT1","HOXA-AS2","PVT1","EGOT","CYTOR"),]
lncTarget_table <- lncTarget_table[lncTarget_table$Score >= score,]

lncTargets <- lncTarget_table[,c("Interactor1","Interactor2")]
colnames(lncTargets) <- c("ont", "gene")


for (i in 1:nrow(clusters_included)){
  cluster_oi <- unlist(strsplit(clusters_included$gene_n[i], ","))

top_enrich[[i]] <- clusterProfiler::enricher(gene=cluster_oi, 
                                        TERM2GENE = lncTargets,
                                        minGSSize = 1,
                                        pvalueCutoff = 0.05,
                                        universe = all_genes,
                                        pAdjustMethod = "none",
                                        qvalueCutoff = 1.0)
}
names(top_enrich) <- clusters_included$color
return(top_enrich)
}

lincr_stats_enrichLnc <- function(cluster, enrichedLnc, corr_cut_off = 0.25){

  clusters_included <- cluster_information[cluster_information$cluster_included == "yes",]
  cluster_number <- match(cluster, clusters_included$color)
  top_pathway <- c()
  top_enrich <- list()
  top_lnc <- c()
  enrichedLnc[[cluster_number]]@result <- enrichedLnc[[cluster_number]]@result[order(enrichedLnc[[cluster_number]]@result$p.adjust, decreasing = F),]
  cluster_genes <- unlist(strsplit(clusters_included$gene_n[cluster_number], ","))
  cluster_genes <- clusterProfiler:: bitr(cluster_genes, 
                                      fromType="SYMBOL", 
                                      toType="ENTREZID", 
                                      OrgDb="org.Hs.eg.db", 
                                      drop = T)
  k <- 0
  lncRNA <- c()
  cor_vec_all <- list()
  all_genes <- unique(c(filt_cutoff_data$V1, filt_cutoff_data$V2))
  all_genes <- clusterProfiler:: bitr(all_genes, 
                                      fromType="SYMBOL", 
                                      toType="ENTREZID", 
                                      OrgDb="org.Hs.eg.db", 
                                      drop = T)
  enrichedLnc[[cluster_number]]@result <- enrichedLnc[[cluster_number]]@result[order(enrichedLnc[[cluster_number]]@result$Count, decreasing = T),]
  for (i in 1:nrow(enrichedLnc[[cluster_number]]@result)){
    if(enrichedLnc[[cluster_number]]@result$p.adjust[i] <= 0.1){

    genes_oi <- unlist(strsplit(enrichedLnc[[cluster_number]]@result$geneID[i],"/"))
    lnc_vector <- c(genes_oi, enrichedLnc[[cluster_number]]@result$ID[i] )
    
    count_file_name <- count_table[rownames(count_table) %in% lnc_vector,]
    
    ds = count_file_name
    dd2 <- head(ds,topvar_genes)
    dd2 = t(dd2)
    voi_id <- "new_cluster"
    corresp_info = info_dataset[rownames(dd2)%in%rownames(info_dataset),]
    corresp_info$grpvar =purrr::pmap(corresp_info[intersect(voi_id,colnames(info_dataset))],
                                       paste, sep="-") %>% unlist()
    
    
    norm_data_anno = merge(dd2, corresp_info["grpvar"], by="row.names",all.x=T)
    norm_data_anno = norm_data_anno[,-1]
    norm_data_anno <- norm_data_anno[ , c(ncol(norm_data_anno) , 1:(ncol(norm_data_anno)-1))]
    
    cor_norm <- norm_data_anno
    cor_norm <- cor_norm[,2:ncol(cor_norm)]
    cor_vec_spearman <- list()
    cor_vec_pearson <- list()
    cor_vec <- c()
    cor_p <- c()
    cor_vec_all_pre <- list()
    for (l in 1:(ncol(cor_norm)-1)){
      cor_vec_spearman[[l]] <- cor.test(x = cor_norm[,ncol(cor_norm)], y= cor_norm[,l],alternative = "two.sided", method = "spearman")
      cor_vec[l] <- cor_vec_spearman[[l]][["estimate"]][["rho"]]
      cor_p[l] <- cor_vec_spearman[[l]][["p.value"]]
      cor_vec_all_pre[[l]]  <- cor_vec_spearman[[l]]
    
    }
    cor_norm <- cor_norm[,1:(ncol(cor_norm)-1)]
    cor <- data.frame(r = cor_vec, p.value = cor_p)
    rownames(cor) <- colnames(cor_norm)
    names(cor_vec_all_pre) <- colnames(cor_norm)
    
    cor <- cor[which(abs(cor$r) >= corr_cut_off),]
    cor <- cor[which(cor$p.value <= 0.1),]
    
    trans_norm <- setNames(data.frame(t(norm_data_anno[ , -1])) , norm_data_anno[,1])
    if(data_in_log==TRUE){
      antilog <- function(lx , base) {
        
        lbx <- lx/log(exp(1) , base = base)
        result <- exp(lbx)
        result
      }
      trans_norm <- antilog(trans_norm , 2)
    }
    
    trans_norm <- t(apply(trans_norm , 1 , function(x) tapply(x , colnames(trans_norm) , mean)))
    trans_norm <- cbind(trans_norm , rowMeans(trans_norm))
    colnames(trans_norm)[ncol(trans_norm)] <- "group_mean"
    grplist=colnames(trans_norm)[-(ncol(trans_norm))]
    
    
    group_means= trans_norm[,"group_mean"]
    
    gfc_calc=function(grp){
      df1= trans_norm[,grp]
      df2=gtools::foldchange(df1, group_means) %>% ifelse(.>2.0, 2.0,.) %>% ifelse(.< (-2.0), -2.0,.) %>% as.data.frame()
      colnames(df2) = paste0("GFC_",grp)
      return(df2)
    }
    
    GFC_all_genes=do.call("cbind", lapply(grplist, gfc_calc))
    GFC_all_genes= round(GFC_all_genes,3)
    GFC_all_genes$Gene = rownames(GFC_all_genes)
    
    if (length(cor_vec) > 0){
    genes_oi <- genes_oi[genes_oi %in% rownames(cor)]
    GFC_pathway <- GFC_all_genes[GFC_all_genes$Gene %in% c(rownames(cor_vec),enrichedLnc[[cluster_number]]@result$ID[i]),]
    
      k <- k + 1
      lncRNA[k] <- enrichedLnc[[cluster_number]]@result$ID[i]
      cor_vec_all[[k]] <- cor

    }

  }
  }
  names(cor_vec_all) <- lncRNA
  return(cor_vec_all)
}

figure_2F <- function(){
  lnc_vector <- c("CYTOR", "VIM","PIK3CB")
  count_file_name <- norm_anno[norm_anno$SYMBOL %in% lnc_vector,]
  count_file_name <- count_file_name[,c(2:(ncol(count_file_name)-3))]
  count_file_name_symbol <- count_file_name[!duplicated(count_file_name[,ncol(count_file_name)]), ]
  count_file_name_symbol <- count_file_name_symbol[complete.cases(count_file_name_symbol), ]
  row.names(count_file_name_symbol) <- count_file_name_symbol[,ncol(count_file_name)]
  count_file_name <- count_file_name_symbol[,c(-(ncol(count_file_name)))]
  
  ds = count_file_name
  dd2 <- head(ds,topvar_genes)
  dd2 = (t(dd2))
  
  corresp_info = info_dataset[rownames(dd2)%in%rownames(info_dataset),]
  voi_id = "new_cluster"
  
  corresp_info$grpvar =purrr::pmap(corresp_info[intersect(voi_id,colnames(info_dataset))],
                                   paste, sep="-") %>% unlist()
  norm_data_anno = merge(dd2, corresp_info["grpvar"], by="row.names",all.x=T)
  norm_data_anno = norm_data_anno[,-1]
  norm_data_anno <- norm_data_anno[ , c(ncol(norm_data_anno) , 1:(ncol(norm_data_anno)-1))]
  for (i in 1:nrow(norm_data_anno)){
    if (norm_data_anno$grpvar[i] ==2){
      norm_data_anno$grpvar[i] = "G1"
    } else if(norm_data_anno$grpvar[i] == 4){
      norm_data_anno$grpvar[i] = "G2"
    } else if(norm_data_anno$grpvar[i] == 6){
      norm_data_anno$grpvar[i] = "G3"
    }else if(norm_data_anno$grpvar[i] == 1){
      norm_data_anno$grpvar[i] = "G4"
    }else if(norm_data_anno$grpvar[i] == 5){
      norm_data_anno$grpvar[i] = "G5"
    }else if(norm_data_anno$grpvar[i] == 3){
      norm_data_anno$grpvar[i] = "G6"
    }
  }
  norm_data_anno$grpvar <- factor(norm_data_anno$grpvar, levels = c("G1","G2","G3","G4","G5","G6"))
  rownames(norm_data_anno) <- info_dataset$ID
  
  norm_data_anno <- norm_data_anno[order(norm_data_anno$grpvar, decreasing = F),]
  
  norm_anno_matrix <- t(scale(as.matrix((norm_data_anno[,2:ncol(norm_data_anno)])), center = T))
  
  test <- plot_annotation[,2,drop = F]
  test$new_cluster <- as.character(levels(test$new_cluster)[test$new_cluster])
  
  for (i in 1:nrow(test)){
    if (test$new_cluster[i] ==2){
      test$new_cluster[i] = "G1"
    } else if(test$new_cluster[i] == 4){
      test$new_cluster[i] = "G2"
    } else if(test$new_cluster[i] == 6){
      test$new_cluster[i] = "G3"
    }else if(test$new_cluster[i] == 1){
      test$new_cluster[i] = "G4"
    }else if(test$new_cluster[i] == 5){
      test$new_cluster[i] = "G5"
    }else if(test$new_cluster[i] == 3){
      test$new_cluster[i] = "G6"
    }
  }
  test$new_cluster <- factor(test$new_cluster, levels = c("G1","G2","G3","G4","G5","G6"))
  colnames(test) <- paste0("Patient Group")
  breakList <- seq(-2, 4, by = .1)
  ann_colors2 <- list(`Patient Group` = c(G1 = "#BC3C29FF", G2 = "#0072B5FF", G3 = "#E18727FF",
                                          G4 = "#20854EFF", G5 = "#7876B1FF", G6 = "#6F99ADFF" ))
  palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(length(seq(-2, 4, by = .1)))
  p <- pheatmap(mat = norm_anno_matrix,
                color = palette,
                scale ="row",
                cellheight = 15,
                cellwidth = 15,
                annotation_col = test,
                annotation_colors = ann_colors2,
                border_color = "black",
                annotation_names_col = F,
                legend_labels = "Group",
                cluster_rows = T,
                cluster_cols = F)
  print(p)
}

getStats <- function(){

  library(reshape2)
  library(ggpubr)
  library(rstatix)
  
lnc_vector <- c("CYTOR", "VIM","PIK3CB")
count_file_name <- norm_anno[norm_anno$SYMBOL %in% lnc_vector,]
count_file_name <- count_file_name[,c(2:(ncol(count_file_name)-3))]
count_file_name_symbol <- count_file_name[!duplicated(count_file_name[,ncol(count_file_name)]), ]
count_file_name_symbol <- count_file_name_symbol[complete.cases(count_file_name_symbol), ]
row.names(count_file_name_symbol) <- count_file_name_symbol[,ncol(count_file_name)]
count_file_name <- count_file_name_symbol[,c(-(ncol(count_file_name)))]

ds = count_file_name
dd2 <- head(ds,topvar_genes)
dd2 = (t(dd2))

corresp_info = info_dataset[rownames(dd2)%in%rownames(info_dataset),]
voi_id = "new_cluster"

corresp_info$grpvar =purrr::pmap(corresp_info[intersect(voi_id,colnames(info_dataset))],
                                 paste, sep="-") %>% unlist()
norm_data_anno = merge(dd2, corresp_info["grpvar"], by="row.names",all.x=T)
norm_data_anno = norm_data_anno[,-1]
norm_data_anno <- norm_data_anno[ , c(ncol(norm_data_anno) , 1:(ncol(norm_data_anno)-1))]
for (i in 1:nrow(norm_data_anno)){
  if (norm_data_anno$grpvar[i] ==2){
    norm_data_anno$grpvar[i] = "G1"
  } else if(norm_data_anno$grpvar[i] == 4){
    norm_data_anno$grpvar[i] = "G2"
  } else if(norm_data_anno$grpvar[i] == 6){
    norm_data_anno$grpvar[i] = "G3"
  }else if(norm_data_anno$grpvar[i] == 1){
    norm_data_anno$grpvar[i] = "G4"
  }else if(norm_data_anno$grpvar[i] == 5){
    norm_data_anno$grpvar[i] = "G5"
  }else if(norm_data_anno$grpvar[i] == 3){
    norm_data_anno$grpvar[i] = "G6"
  }
}
norm_data_anno$grpvar <- factor(norm_data_anno$grpvar, levels = c("G1","G2","G3","G4","G5","G6"))
rownames(norm_data_anno) <- info_dataset$ID
norm_data_anno <- norm_data_anno[order(norm_data_anno$grpvar, decreasing = F),]
trans_norm <- setNames(data.frame(t(norm_data_anno[ , -1])) , norm_data_anno[,1])
if(data_in_log==TRUE){
  antilog <- function(lx , base) {
    
    lbx <- lx/log(exp(1) , base = base)
    result <- exp(lbx)
    result
  }
  trans_norm <- antilog(trans_norm , 2)
}

trans_norm <- t(apply(trans_norm , 1 , function(x) tapply(x , colnames(trans_norm) , mean)))
trans_norm <- cbind(trans_norm , rowMeans(trans_norm))
colnames(trans_norm)[ncol(trans_norm)] <- "group_mean"
grplist=colnames(trans_norm)[-(ncol(trans_norm))]


group_means= trans_norm[,"group_mean"]

gfc_calc=function(grp){
  df1= trans_norm[,grp]
  df2=gtools::foldchange(df1, group_means) %>% ifelse(.>2.0, 2.0,.) %>% ifelse(.< (-2.0), -2.0,.) %>% as.data.frame()
  colnames(df2) = paste0("GFC_",grp)
  return(df2)
}

GFC_all_genes=do.call("cbind", lapply(grplist, gfc_calc))
GFC_all_genes= round(GFC_all_genes,3)
GFC_all_genes$Gene = rownames(GFC_all_genes)
GFC_mat <- as.matrix(GFC_all_genes[,1:6])

GFC_oi <- melt(norm_data_anno, id.vars = c("grpvar"))
tmp <- list()

for (i in 1:length(lnc_vector)){
  GFC_oi_sub <- GFC_oi[GFC_oi$variable == lnc_vector[i],]
  tmp[[i]] <- aov(value ~ grpvar, data = GFC_oi_sub) %>% tukey_hsd()
  tmp[[i]] <- tmp[[i]][tmp[[i]]$group2 == "G6",]
}
names(tmp) <- lnc_vector
return(tmp)
}