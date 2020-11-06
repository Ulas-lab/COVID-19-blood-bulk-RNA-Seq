# Remove existing variables from R memory
rm(list=ls())

# Install packages
source("https://bioconductor.org/biocLite.R")
biocLite('GSVA')

# Load packages
library(GSVA)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(viridis)
library(dplyr)
library(readxl)

# Set working directory
setwd('D:/Theo/Documents/LIMES 2019-20+/Collaborations/Figure 3')
getwd()

# Load and prepare gene sets (Fig. S4D)
cluster.0<-read.csv('20200603_intersect_0_unique.csv')
cluster.2<-read.csv('20200603_intersect_2_unique.csv')
cluster.89<-read.csv('20200602_Seurat_markers_89.csv')

# Subset gene sets
cluster.0<-cluster.0$x
cluster.2<-cluster.2$x
cluster.89<-subset(cluster.89, cluster.89$p_val_adj<0.01 & cluster.89$avg_logFC>0)

list2 = list(as.character(cluster.0),as.character(cluster.2),as.character(cluster.89$X))
names(list2)<-c('Cluster.0','Cluster.2','Cluster.89')

# Load normalised data
norm_anno<-read.csv('normalized counts.csv')
norm_anno<-norm_anno[!duplicated(norm_anno$SYMBOL),]
rownames(norm_anno)<-norm_anno$SYMBOL
norm_anno<-norm_anno[,c(3:46)]
norm_anno<-as.matrix(norm_anno)

# Load annotation
sample_table <- read_xlsx("Fig3 201007 anno for granulos_final.xlsx")
sample_table <- sample_table[,c(1,2,5,6,7,8,9,10,11,13,14,15,16,17,19,20,26,28,30,33)]
sample_table$Donor<-substr(sample_table$Sample_Name_01, start=1, stop=4)

sample_table <- sample_table[which(!sample_table$ID %in% c(9221,9238,9242,9246,9247,9250,9251,9259,9260,9263,9266,9268,9270,9271,9272,9455)),]
sample_table <- sample_table[which(!sample_table$ID %in% sample_table[sample_table$Diagnosis=='COVID-19-neg',]$ID),]

sample_table$severity_phase<-factor(sample_table$severity_phase,levels=c("mild_early","mild_late","severe_early","severe_late"))
sample_table$ID<-paste0('X',sample_table$ID)

# Run GSVA
es<-gsva(norm_anno,
         list2,
         min.sz=1,
         max.sz=1081,
         verbose=T,
         method = 'gsva',
         #kcdf = 'Poisson',           # set to 'Poisson' when counts are integers
         mx.diff=F,                  # set to TRUE if you want to do downstream DE gene comparisons
         parallel.sz=1)

# Re-order count table
ordered<-sample_table %>% arrange(severity_phase, Days_post_1st_symptoms)
es<-es[,ordered$ID]
es<-data.frame(rowMeans(es[,ordered[ordered$severity_phase=="mild_early",]$ID]),
               rowMeans(es[,ordered[ordered$severity_phase=="mild_late",]$ID]),
               rowMeans(es[,ordered[ordered$severity_phase=="severe_early",]$ID]),
               rowMeans(es[,ordered[ordered$severity_phase=="severe_late",]$ID]))

es        
colnames(es)<-levels(sample_table$severity_phase)

# Heatmap annotation
df<-data.frame(severity_phase=colnames(es))
rownames(df)<-df$severity_phase

df[] <- lapply(df, factor)

col_severity_phase<-tableau_color_pal(palette='Tableau 10')(4) 
names(col_severity_phase)<-levels(sample_table$severity_phase)

pheatmap(es,
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         scale='row',
         show_rownames = T,
         show_colnames = F,
         cluster_rows = T,
         cluster_cols = F,
         cellheight=25,
         cellwidth = 35,
         main='',
         annotation_col = df,
         annotation_colors = list(severity_phase=col_severity_phase),
         color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(length(seq(-2, 2, by = .1))))

## Individual gene boxplots (Fig. 3D)
# Set working directory
setwd('D:/Theo/Documents/LIMES 2019-20+/Collaborations/Figure 3')
getwd()

# Load and prepare signatures
genes<-c('CD177','S100A6','MMP9')

# Load normalised data
norm_anno<-read.csv('normalized counts.csv')
norm_anno<-norm_anno[!duplicated(norm_anno$SYMBOL),]
rownames(norm_anno)<-norm_anno$SYMBOL
norm_anno<-norm_anno[,c(3:46)]
norm_anno<-norm_anno[genes,]
norm_anno<-as.matrix(norm_anno)

# Load annotation
sample_table <- read_xlsx("Fig3 201007 anno for granulos_final.xlsx")
sample_table <- sample_table[,c(1,2,5,6,7,8,9,10,11,13,14,15,16,17,19,20,26,28,30,33)]
sample_table$Donor<-substr(sample_table$Sample_Name_01, start=1, stop=4)

sample_table <- sample_table[which(!sample_table$ID %in% c(9221,9238,9242,9246,9247,9250,9251,9259,9260,9263,9266,9268,9270,9271,9272,9455)),]
sample_table <- sample_table[which(!sample_table$ID %in% sample_table[sample_table$Diagnosis=='COVID-19-neg',]$ID),]

sample_table$severity_phase<-factor(sample_table$severity_phase,levels=c("mild_early","mild_late","severe_early","severe_late"))
sample_table$ID<-paste0('X',sample_table$ID)

# Re-order count table
ordered<-sample_table %>% arrange(severity_phase, Days_post_1st_symptoms)

# Plot expression boxplots
table(ordered$severity_phase)
norm_anno<-melt(norm_anno)
norm_anno$grouping<-ordered$severity_phase[match(norm_anno$Var2,ordered$ID)]

ggplot(norm_anno,aes(x=grouping,y=value,color=grouping))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.2)+
  geom_boxplot()+
  geom_jitter(width=0.2)+
  facet_wrap(~Var1,scale='free',ncol=3)+
  scale_color_manual(values=c("#A49DC7","#A49DC7","#522685","#522685"))+
  ylab('Normalized counts')+
  guides(fill=F)+
  theme(plot.title=element_text(hjust=0.5,face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14,face='bold'),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=14,colour='black'),
        strip.text=element_text(size=12),
        strip.background = element_rect(fill='gray74',colour='black'),
        legend.title=element_text(size=13,face='bold',hjust=1),
        legend.key=element_rect(colour='black',fill='white'),
        panel.background=element_rect(fill='white'),
        panel.border=element_rect(colour='gray40',fill=NA))

## Individual gene boxplots (Fig. 3G, Fig. S5D)
# Set working directory
setwd('D:/Theo/Documents/LIMES 2019-20+/Collaborations/Figure 3')
getwd()

# Load and prepare signatures
genes<-c('ARG1','CD274','NLRC4')

# Load normalised data
norm_anno<-read.csv('normalized counts.csv')
norm_anno<-norm_anno[!duplicated(norm_anno$SYMBOL),]
rownames(norm_anno)<-norm_anno$SYMBOL
norm_anno<-norm_anno[,c(3:46)]
norm_anno<-norm_anno[genes,]
norm_anno<-as.matrix(norm_anno)

# Load annotation
sample_table <- read_xlsx("Fig3 201007 anno for granulos_final.xlsx")
sample_table <- sample_table[,c(1,2,5,6,7,8,9,10,11,13,14,15,16,17,19,20,26,28,30,33)]
sample_table$Donor<-substr(sample_table$Sample_Name_01, start=1, stop=4)

sample_table <- sample_table[which(!sample_table$ID %in% c(9221,9238,9242,9246,9247,9250,9251,9259,9260,9263,9266,9268,9270,9271,9272,9455)),]
sample_table <- sample_table[which(!sample_table$ID %in% sample_table[sample_table$Diagnosis=='COVID-19-neg',]$ID),]

sample_table$severity_phase<-factor(sample_table$severity_phase,levels=c("mild_early","mild_late","severe_early","severe_late"))
sample_table$ID<-paste0('X',sample_table$ID)
sample_table$`CoCena-G`<-c('G2','G2','G2','G2','G1','G1','G1','G5','G4','G2','G2','G1','G1','G1_','G4','G2','G5','G1','G2','G2','G2','G1_','G1','G1','G1',
                           'G5','G4','G2','G3','G3','G3','G3','G3','G3','G1','G2','G3','G2','G2','G1_','G1','G4','G1_','G2')

# Plot expression boxplots
norm_anno<-melt(norm_anno)
norm_anno$grouping<-sample_table$`CoCena-G`[match(norm_anno$Var2,sample_table$ID)]

ggplot(norm_anno,aes(x=grouping,y=value,color=grouping))+
  stat_boxplot(geom='errorbar', linetype=1, width=0.2)+
  geom_boxplot()+
  geom_jitter(width=0.2)+
  facet_wrap(~Var1,scale='free',ncol=3)+
  scale_color_nejm()+
  ylab('Normalized counts')+
  guides(fill=F)+
  theme(plot.title=element_text(hjust=0.5,face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=14,face='bold'),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=14,colour='black'),
        strip.text=element_text(size=12),
        strip.background = element_rect(fill='gray74',colour='black'),
        legend.title=element_text(size=13,face='bold',hjust=1),
        legend.key=element_rect(colour='black',fill='white'),
        panel.background=element_rect(fill='white'),
        panel.border=element_rect(colour='gray40',fill=NA))

## Individual genes from DECOI (Fig. 3F)
# Set working directory
setwd('D:/Theo/Documents/LIMES 2019-20+/Collaborations/Figure 3')
getwd()

# Load and prepare signatures
genes<-c('FUT4',
         'CD63',
         'MMP8',
         'S100A6',
         'S100A8',
         'S100A9',
         'GSN',
         'NLRC4',
         'PADI4',
         'CXCL8',
         'ITGA4',
         'SLC38A1',
         'MME',
         'PTPRC',
         'IFI6',
         'IFIT1',
         'IFIT3',
         'IFITM1',
         'IFITM2',
         'ISG15',
         'ARG1',
         'CD177',
         'CD274',
         'ZDHHC19')

#'IGSF2','IL3RA', 'MARC1' not expressed

# Load normalised data
norm_anno<-read.csv('normalized counts.csv')
norm_anno<-norm_anno[!duplicated(norm_anno$SYMBOL),]
rownames(norm_anno)<-norm_anno$SYMBOL
norm_anno<-norm_anno[,c(3:46)]
norm_anno<-norm_anno[genes,]
norm_anno<-as.matrix(norm_anno)

# Load annotation
sample_table <- read_xlsx("Fig3 201007 anno for granulos_final.xlsx")
sample_table <- sample_table[,c(1,2,5,6,7,8,9,10,11,13,14,15,16,17,19,20,26,28,30,33)]
sample_table$Donor<-substr(sample_table$Sample_Name_01, start=1, stop=4)

sample_table <- sample_table[which(!sample_table$ID %in% c(9221,9238,9242,9246,9247,9250,9251,9259,9260,9263,9266,9268,9270,9271,9272,9455)),]
sample_table <- sample_table[which(!sample_table$ID %in% sample_table[sample_table$Diagnosis=='COVID-19-neg',]$ID),]

sample_table$severity_phase<-factor(sample_table$severity_phase,levels=c("mild_early","mild_late","severe_early","severe_late"))
sample_table$ID<-paste0('X',sample_table$ID)

# Re-order count table
ordered<-sample_table %>% arrange(severity_phase, Days_post_1st_symptoms)
norm_anno<-norm_anno[,ordered$ID]
norm_anno<-data.frame(rowMeans(norm_anno[,ordered[ordered$severity_phase=="mild_early",]$ID]),
                      rowMeans(norm_anno[,ordered[ordered$severity_phase=="mild_late",]$ID]),
                      rowMeans(norm_anno[,ordered[ordered$severity_phase=="severe_early",]$ID]),
                      rowMeans(norm_anno[,ordered[ordered$severity_phase=="severe_late",]$ID]))

norm_anno        
colnames(norm_anno)<-c('mild_early', 'mild_late', 'severe_early', 'severe_late')

# Heatmap annotation
df<-data.frame(severity_grouping=colnames(norm_anno))
rownames(df)<-df$severity_grouping

df[] <- lapply(df, factor)

col_severity_grouping<-tableau_color_pal(palette='Tableau 10')(4) 
names(col_severity_grouping)<-c('mild_early', 'mild_late', 'severe_early', 'severe_late')

pheatmap(norm_anno[genes,c(4,3,2,1)],
         clustering_distance_rows = 'correlation',
         clustering_distance_cols = 'correlation',
         scale='row',
         show_rownames = T,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         cellheight=12,
         cellwidth=15,
         main='',
         annotation_col = df,
         fontsize = 7,
         annotation_colors = list(severity_grouping=col_severity_grouping),
         color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(length(seq(-2, 2, by = .1))))

## Overlay of DECOI genes and DE genes from PAX study
# Load and prepare gene sets
signature.6<-read.csv('Granulocyte-Cluster-Marker_UPDATE2.csv')
cocena<-read.delim('network_modules.txt',header=T)
upgenes<-read.csv('upgenes.csv')
upgenes<-upgenes[,c(3,4,7,8)]
downgenes<-read.csv('downgenes.csv')
downgenes<-downgenes[,c(3,4,7,8)]

upgenes$DECOI<-'No cluster'
upgenes$CoCena<-'No module'

for (i in intersect(upgenes$SYMBOL,signature.6$Marker)){
  upgenes[upgenes$SYMBOL==i,]$DECOI<-as.character(signature.6[signature.6$Marker==i,]$cluster)
}

for (i in intersect(upgenes$SYMBOL,cocena$gene_n)){
  upgenes[upgenes$SYMBOL==i,]$CoCena<-as.character(cocena[cocena$gene_n==i,]$color)
}

downgenes$DECOI<-'No cluster'
downgenes$CoCena<-'No module'

for (i in intersect(downgenes$SYMBOL,signature.6$Marker)){
  downgenes[downgenes$SYMBOL==i,]$DECOI<-as.character(signature.6[signature.6$Marker==i,]$cluster)
}

for (i in intersect(downgenes$SYMBOL,cocena$gene_n)){
  downgenes[downgenes$SYMBOL==i,]$CoCena<-as.character(cocena[cocena$gene_n==i,]$color)
}

write.csv(upgenes,'updated_upregulated neutrophils.csv')
write.csv(downgenes,'updated_downregulated neutrophils.csv')
