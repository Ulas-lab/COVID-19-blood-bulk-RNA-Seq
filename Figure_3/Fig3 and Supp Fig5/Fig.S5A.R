# Remove existing variables from R memory
rm(list=ls())

# Load packages
library(reshape2)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(RColorBrewer)
library(readxl)

# Set working directory
setwd("D:/Theo/Documents/LIMES 2019-20+/Collaborations/Figure 3/")
getwd()

## Fig. S4A
# Load the datasets
samples<-read.csv('CIBERSORTx_Job2_Results.csv')
anno<-read_xlsx('Fig3 201007 anno for granulos_final.xlsx')
anno$phase_days<-paste(anno$phase,sep='.',anno$Days_post_1st_symptoms)
anno<-anno[anno$ID %in% samples$Mixture,]

# Get cumulative results per cell type
data<-data.frame(patient=samples[,1], lymphocytes=rowSums(samples[,c(2:12)]),monocytes.dcs=rowSums(samples[,c(13:15)]),eosinophils=samples[,16],neutrophils=samples[,17])

data<-data[order(data$patient,decreasing=F),]
anno<-anno[order(anno$ID,decreasing=F),]

data<-cbind(data,anno[,11])

colnames(data)<-c('Patient.number','Lymphocytes','Monocytes/DCs','Eosinophils','Neutrophils','Phase.Days.post.onset')

# Prepare the data for lineage resolution
data<-melt(data,id.vars=c('Patient.number','Phase.Days.post.onset'))
data$Patient.number<-factor(data$Patient.number)
data$variable<-factor(data$variable, levels=c('Neutrophils','Eosinophils','Monocytes/DCs','Lymphocytes'))
data$value100<-data$value*100

# Plot boxplot
ggplot(data,aes(x=Patient.number,y=value100,fill=variable))+
  geom_bar(stat='identity',position=position_stack(rev=T),width=1)+
  scale_fill_manual(values=c('#C7B1D3','#20854E','#EF7901','#0072B5'))+
  ylab('Percentage of immune cells')+
  xlab('')+
  ggtitle('')+
  guides(fill=F)+
  theme(plot.title = element_text(size=16,colour='black',face='bold',hjust = 0.5),
        axis.text.y=element_text(size=14,colour='black'),
        axis.title.y=element_text(size=14,face='bold'),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_text(size=14,face='bold'),
        strip.text=element_text(size=12),
        strip.background = element_rect(fill='gray84',colour='black'),
        legend.title=element_blank(),
        legend.position = 'right',
        legend.key=element_rect(colour='black',fill='white'),
        legend.spacing.x = unit(0.3, 'cm'),
        panel.background=element_rect(fill='white'),
        panel.border=element_rect(colour='black',fill=NA))

ggplot(data,aes(x=variable,y=value100,fill=variable))+
  stat_summary(geom = "errorbar", fun.data = mean_se, position = position_dodge(0.9), width=0.7, cex=0.8)+
  stat_summary(geom = "bar", fun.y = mean, position = position_dodge(0.9), colour='black') +
  ylab('Percentage of immune cells')+
  xlab('')+
  ggtitle('')+
  guides(fill=F)+
  theme(plot.title = element_text(size=16,colour='black',face='bold',hjust = 0.5),
        axis.text.y=element_text(size=14,colour='black'),
        axis.title.y=element_text(size=14,face='bold'),
        axis.text.x=element_text(size=14,colour='black', angle = 270),
        axis.title.x=element_text(size=14,face='bold'),
        strip.text=element_text(size=12),
        strip.background = element_rect(fill='gray84',colour='black'),
        legend.title=element_blank(),
        legend.position = 'right',
        legend.key=element_rect(colour='black',fill='white'),
        legend.spacing.x = unit(0.3, 'cm'),
        panel.background=element_rect(fill='white'),
        panel.border=element_rect(colour='black',fill=NA))
