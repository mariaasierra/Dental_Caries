#picrust analysis
library(microbiome)
library(ggplot2)
library(ALDEx2)
library("readxl")
library("magrittr")
library(tidyverse)
#import abundance table
data<- read_xlsx("pathway.xlsx")
head(data)

#merge with meatdata
metadata=read.table("metadata.txt", header = T)
master_table=merge(data, metadata,by = "SampleID")
head(metadata)

result <- filter(master_table, type.plaque == "infected") #Select plaque infected
result2 <- result[,-1]
rownames(result2) <- result[,1]
result_trim <- result2 %>% dplyr::select(-sample.type, -treatment, -type.plaque)

metadata_order <- filter(metadata, type.plaque == "infected") #Select plaque infected
newdata <- metadata_order[order(metadata_order$treatment),] 
rownames(newdata) <- newdata[,1]

#PCA
library(devtools)
pathway_t<- t(result_trim)

data_t=cbind(pathway_t, newdata)

data_t$SampleID <- NULL

pca_res <- prcomp(data_t[1:256], scale. = TRUE)#choose only numeric variable
library(ggplot2)
library(grid)
library(gridExtra)
pdf("pcoa.path.pdf", height = 7, width = 7)
autoplot(pca_res, data = data_t,col="treatment",size = 2.8, frame = TRUE, frame.type = 'norm')+ theme_minimal()+scale_colour_manual(values=c("yellow4", "deepskyblue2", "hotpink1", "springgreen3"))+
  scale_fill_manual(values=c("yellow4", "deepskyblue2", "hotpink1", "springgreen3"))+
  scale_shape_manual(values=c(25,22,23))
dev.off()

#boxplot
#picked_pathways
picked_data<- read_excel("picked.xlsx")


#transforming data to long format
picked_long <- pivot_longer(data = picked_data, 
                            cols = -c(1:1), 
                            names_to = "Class", 
                            values_to = "Abundance")
head(picked_long)
#merge with matadata
picked_table=merge(picked_long, metadata,by = "SampleID")
pdf("picked_boxplot.pdf")
ggplot(picked_table, aes(x=treatment, y=Abundance, fill=treatment)) +geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3)+ylab("Abundance")+facet_wrap(~Class,scales = "free")+ theme(axis.text.x=element_blank(),
                      axis.ticks.x=element_blank())+geom_boxplot(alpha=0.6)+stat_summary(fun.y=mean, geom="point", shape=23, size=2)+
  theme_minimal()+scale_colour_manual(values=c("yellow4", "deepskyblue2", "hotpink1", "springgreen3"))+
  scale_fill_manual(values=c("yellow4", "deepskyblue2", "hotpink1", "springgreen3"))+theme(legend.position = "right",strip.text = element_text(size=10, face="bold")) 
print(pathway_plot)
dev.off()
