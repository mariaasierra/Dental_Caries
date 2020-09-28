library(compositions)
library(ggsignif)
library(ggpubr)
library(phyloseq)
library(data.table)
library(tidyr)
library(ggplot2)
library("ALDEx2")
library(microbiome)
library(ggtree)
library(vegan)
library(pheatmap)

#Import files
OTUs=read.table("../Data/table.frombiom.txt", header = T,sep="\t",stringsAsFactors = FALSE)
OTUs=t(OTUs)
colnames(OTUs) <- OTUs[1,]
OTUs <- OTUs[-1, ] 
class(OTUs) <- "numeric"
TAX=read.table("../Data/taxonomy.tsv", header = T,sep="\t")
row.names(TAX) = TAX$Feature.ID
TAX=TAX[,-1]
tax.clean = separate(TAX, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
tax.clean = tax.clean[,-which(names(tax.clean) %in% c("Confidence", "OTU"))]
OTU.UF = otu_table(OTUs, taxa_are_rows=F)
tax.UF = tax_table(as.matrix(tax.clean))
meta = read.table("../Data/metadata.txt", header=TRUE, row.names=1, sep="\t", dec = ".") #Sample names MUST be the same in both files, otherwise is going not going to pair data 
meta<-sample_data(meta)
tree <- read_tree ("../Data/tree.nwk")
physeq = phyloseq(OTU.UF,tax.UF, meta, tree)
sdf_phylo <- merge_phyloseq (OTU.UF,tax.UF, meta, tree)

#ALDEx2 comparisons by treatment

#Plaque Infected
OTU.UF = as.matrix(otu_table(plq_inf_phylo, taxa_are_rows=F))
tax= as.matrix(tax_table(plq_inf_phylo@tax_table@.Data))
gg.otus <- t(OTU.UF) #Transpose OTU
ann<- read.table("../Data/annotations.txt", header = T)
col_order <- row.names(ann)
OTU.UF.ord <- gg.otus[, col_order] #Reorder columns
sample_data(plq_inf_phylo)

conds <- c(rep("caries-active", 5), rep("caries-free", 5), rep("non-responders", 5),
                 rep("SDF",5)) #Treatment are the conditions to compare
x<- aldex.clr(OTU.UF.ord, conds, denom = "all")#clr: centred log-ratio transformation
x.al <- aldex.kw(x)#test glm: general linear model and Kruskal Wallace
OTU.aldex = x.al[which(x.al$kw.ep < 0.05),] #OTUs with p<0.05  #glm ANOVA

#Generate new phyloseq file with OTUs differentially abundant between treatment p<0.05 from ALDEx2
OTU.aldex$OTU <- rownames(OTU.aldex)
taxa.pvalues.aldex <- subset(tax, subset=rownames(tax) %in% OTU.aldex$OTU, select = c("Genus", "Species")) ##Extract Taxonomy for ALDEx2 OTUs
row.names(taxa.pvalues.aldex) = OTU.aldex$OTU
tax.UF = tax_table(as.matrix(taxa.pvalues.aldex))
OTU.clean<-OTU.UF.ord[row.names(OTU.aldex),] #Get abundance values of ALDEx2 OTUs only
OTU.aldex.uf = otu_table(as.matrix(OTU.clean), taxa_are_rows=F)
meta<-sample_data(ann)
physeq = phyloseq(OTU.aldex.uf, tax.UF, meta)
physeq.aldex = merge_phyloseq(physeq) ## New phyloseq file with ALDEx2 OTUs

#Plot significant OTUs

table_otus=t(OTU.clean)
row_order <- row.names(ann)
table_ord <- table_otus[row_order,] 
table_active_merged=cbind(table_ord, ann)

#Tannerella_sp
table_active_merged=table_active_merged %>% 
  rename(Tannerella_sp="9b37d7b16381c8c3b706e4ef7c3feeba")
compare_means(Tannerella_sp ~ Description,  data = table_active_merged)
ggplot(table_active_merged, aes(x=Description,y=Tannerella_sp, fill=Description))+
  geom_violin(alpha=0.5)+geom_dotplot(binaxis='y', stackdir='center',
                              position=position_dodge(1))+theme_minimal(base_size = 20, base_line_size = 0.5) + scale_fill_manual(values = c("yellow4","deepskyblue2",
                                                                                                                                              "hotpink1", "springgreen3"))

#Lachnospiraceae
table_active_merged=table_active_merged %>% 
  rename(Lachnospiraceae="ab4970076c2f42193e5c6dce685e7127")
compare_means(Lachnospiraceae ~ Description,  data = table_active_merged)
ggplot(table_active_merged, aes(x=Description,y=Lachnospiraceae, fill=Description))+
  geom_violin(alpha=0.5)+geom_dotplot(binaxis='y', stackdir='center',
                              position=position_dodge(1))+theme_minimal(base_size = 20, base_line_size = 0.5)+ scale_fill_manual(values = c("yellow4","deepskyblue2",
                                                                                                                                            "hotpink1", "springgreen3"))


#Granulicatella_adiacens
table_active_merged=table_active_merged %>% 
  rename(Granulicatella_adiacens="88286b80fadbc5732691bc1abd772563")
compare_means(Granulicatella_adiacens ~ Description,  data = table_active_merged)
ggplot(table_active_merged, aes(x=Description,y=Granulicatella_adiacens, fill=Description))+
  geom_violin(alpha=0.5)+geom_dotplot(binaxis='y', stackdir='center',
                              position=position_dodge(1))+theme_minimal(base_size = 20, base_line_size = 0.5)+ scale_fill_manual(values = c("yellow4","deepskyblue2",
                                                                                                                                                                                       "hotpink1", "springgreen3"))


