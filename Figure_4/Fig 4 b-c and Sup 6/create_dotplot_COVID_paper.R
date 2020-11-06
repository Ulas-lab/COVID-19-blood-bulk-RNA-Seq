library(ggplot2)
library(magrittr)
require("GSVA")
require("reshape2")
library(clusterProfiler)
library(ggpubr)
library(data.table)
library(gtools)
library(readr)
library(stringr)


# load data and annotatin (data is normalized, filtered for the COVID-19 network genes, genes coding for ribosomal proteins as
# well as genes encoded by mitochondrial DNA have been removed) 

ext_data <- read.csv("INSERT PATH TO: norm_counts_covgenes_noRP_noMT_noMEN_noFAT.txt", sep = "\t", row.names = 1, 
                     check.names = F)

norm_counts_anno <- read_delim("INSERT PATH TO: norm_anno_covgenes_noRP_noMT_noMEN_noFAT.txt", 
                               "\t", escape_double = FALSE, trim_ws = TRUE)
norm_counts_anno <- norm_counts_anno[grepl("control", norm_counts_anno$merged) == F,]


# calculate mean expressions for each condition
meanedexp <- data.frame(gene = rownames(ext_data))
for(x in unique(norm_counts_anno$merged)){
  tmpanno <- norm_counts_anno[norm_counts_anno$merged == x, ]
  tmpcount <- ext_data[, colnames(ext_data) %in% tmpanno$ID]
  tmpdf <- data.frame(V1 = apply(tmpcount, 1, mean))
  colnames(tmpdf) <- c(x)
  meanedexp <- cbind(meanedexp, tmpdf)
}
rownames(meanedexp) <- meanedexp$gene
meanedexp$gene <- NULL




data_GSVA <- as.matrix(meanedexp)

# load selected HALLMARKs to be plotted:
HALLMARK_selected <- read.csv("INSERT PATH TO: DJplot_HALLMARK_selected only.csv")
HALLMARK_selected <- unlist(HALLMARK_selected)
HALLMARK_selected <- levels(HALLMARK_selected)
HALLMARK_selected <- HALLMARK_selected[!HALLMARK_selected == ""]

# load gene sets
gene_sets <- read.gmt("INSERT PATH TO: h.all.v6.1.symbols.gmt")
gene_sets <- data.frame(lapply(gene_sets, as.character), stringsAsFactors=FALSE)
gos <- unique(gene_sets$ont)
gene_sets_GSVA = list()
for (i in gos){
  gene_sets_GSVA[[i]] <- gene_sets[ gene_sets$ont == i ,2]
}


# GSVA:
es2 <- gsva(data_GSVA, gene_sets_GSVA, verbose=FALSE, parallel.sz=1, method = "zscore",abs.ranking = T, mx.diff=FALSE)

melted_cormat <- melt(es2)

colv <- c("Covid-19 Cluster 2", "Covid-19 Cluster 4", "Covid-19 Cluster 6", "Covid-19 Cluster 1", "Covid-19 Cluster 5", "Covid-19 Cluster 3",
          "Chikungunya acute", "Zika Early Acute", "Zika Late Acute", "Ebola vaccine",
          "HIV", "Tuberculosis HIV", "Tuberculosis HIV ART", 
          "Tuberculosis Lung", "Tuberculosis MTP", "Tuberculosis", "Latent tuberculosis", "uncomplicated sepsis", "severe sepsis", "septic shock", 
          "sepsis death", "SIRS",  "NLRC4-MAS anakinra treatment", "NLRC4-MAS No treatment", "NOMID anakinra treatment", 
          "NOMID No treatment", "Rheumatoid arthritis Unknown", "Autoimmune disease")


melted_cormat2 <- melted_cormat[melted_cormat$Var1 %in% HALLMARK_selected,]
melted_cormat2$Var2 <- factor(melted_cormat2$Var2, levels = unique(melted_cormat2$Var2[order(match(melted_cormat2$Var2, colv))]))
melted_cormat2$value[melted_cormat2$value > 2] <- 2
melted_cormat2$value[melted_cormat2$value < -2] <- -2

# format names nicely:
nice_names <- melted_cormat2$Var1 %>%
  str_replace(., "HALLMARK_", "")%>%
  str_replace_all(., "_", " ")%>%
  tolower()
nice_names <- str_replace_all(nice_names, "il6", "IL6")
nice_names <- str_replace_all(nice_names, "jak", "JAK")
nice_names <- str_replace_all(nice_names, "stat5", "STAT5")
nice_names <- str_replace_all(nice_names, "tnfa", "TNFA")
nice_names <- str_replace_all(nice_names, "stat3", "STAT3")
nice_names <- str_replace_all(nice_names, "il2", "IL2")
nice_names <- str_replace_all(nice_names, "mtorc1", "mTORC1")
nice_names <- str_replace_all(nice_names, "g2m", "G2M")
nice_names <- str_replace_all(nice_names, "nfkb", "NFKB")
melted_cormat2$Var1 <- nice_names


fucntional_order <- c("G2M checkpoint","mTORC1 signaling" ,"oxidative phosphorylation",  
                      "glycolysis","hypoxia",  "heme metabolism", "estrogen response early",  
                      "coagulation", "complement", "IL2 STAT5 signaling", "interferon alpha response",
                      "interferon gamma response", "TNFA signaling via NFKB","IL6 JAK STAT3 signaling",
                      "inflammatory response" )

# plot:
g <- ggplot(melted_cormat2, aes(x=Var2, y = Var1))+ geom_point(aes(size = value, col= value))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        axis.text = element_text(size = 10))+
  scale_color_gradient2(midpoint=0, low="blue", mid="white",
                        high="red", space ="Lab" ,
                        limits=c(-2, 2),
                        na.value = "black", labels = rev(c(">=2", "1", "0", "-1", "<= -2")))+
  theme(legend.position = "right",
        axis.title = element_blank())+
  scale_y_discrete(limits = fucntional_order, position = "right")



g

# save:

# Cairo::CairoPDF(file = "INSERT PATH: dotplot_with_legend.pdf", width = 15, 
#                 height = 10)
# g
# dev.off()
