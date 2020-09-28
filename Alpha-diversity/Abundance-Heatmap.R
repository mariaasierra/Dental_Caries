##Heatmaps plaque

#Heatmap plaque infected

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

#Top40 taxa species
plq_phylo = subset_samples(sdf_phylo, sample.type == "Plaque") #Subset category
plq_inf_phylo= subset_samples(plq_phylo, type.plaque == "infected") #Subset TypePlaque
plq_inf_sp = tax_glom(plq_inf_phylo, "Species") #Glom at species level
top40 = prune_taxa(names(sort(taxa_sums(plq_inf_sp), TRUE))[1:40], plq_inf_sp) #Only Top 40
plq.trans <- microbiome::transform(top40, "log10") #Transform counts
otu_trans=t(as.data.frame(otu_table(plq.trans)))

ann<- read.table("../Data/annotations.txt", header = T)
col_order <- row.names(ann)
OTUdf.ord <- otu_trans[, col_order] #Rorder columns

tax.bin=tax_table(top40)
tax.bin=as.data.frame(tax.bin@.Data)
ann_plq=tax.bin[,6:7]
ann_plq$names <- paste(ann_plq$Genus, " ", ann_plq$Species)
ann_plq$names<-gsub("g__","",ann_plq$names)
ann_plq$names<-gsub("s__","",ann_plq$names)

pheatmap(as.matrix(OTUdf.ord), border_color = "black",fontsize_row=8, 
         fontsize_col = 8, scale = "none",
         cluster_cols = F, cluster_rows = T,
         labels_col = NULL, annotation_col = ann,
         labels_row=paste0(ann_plq$names))
         
#Top40 taxa genus glom

plq_inf_gen = tax_glom(plq_inf_phylo, "Genus") #Glom at genus level
top40 = prune_taxa(names(sort(taxa_sums(plq_inf_gen), TRUE))[1:40], plq_inf_gen)#Only Top 40
plq.trans <- microbiome::transform(top40, "log10") #Transform
otu_trans=t(as.data.frame(otu_table(plq.trans)))
data.dist <- vegdist(otu_trans, method = "bray")

ann<- read.table("../Data/annotations.txt", header = T)
col_order <- row.names(ann)
OTUdf.ord <- otu_trans[, col_order] #Rorder columns

tax.bin=tax_table(top40)
tax.bin=as.data.frame(tax.bin@.Data)
ann_plq=tax.bin[,6]
ann_plq=data.frame(ann_plq)
ann_plq$ann_plq<-gsub("g__","",ann_plq$ann_plq)
head(ann_plq)

pheatmap(as.matrix(OTUdf.ord), border_color = "black",fontsize_row=8, 
         fontsize_col = 8, scale = "none",
         cluster_cols = F, cluster_rows = T,
         labels_col = NULL, annotation_col = ann,
         labels_row=paste0(ann_plq$ann_plq))

#Top40 taxa Phylum glom
plq_inf_gen = tax_glom(plq_inf_phylo, "Phylum")
top40 = prune_taxa(names(sort(taxa_sums(plq_inf_gen), TRUE))[1:40], plq_inf_gen)
plq.trans <- microbiome::transform(top40, "log10")
otu_trans=t(as.data.frame(otu_table(plq.trans)))
data.dist <- vegdist(otu_trans, method = "bray")

ann<- read.table("../Data/annotations.txt", header = T)
col_order <- row.names(ann)
OTUdf.ord <- otu_trans[, col_order] #Rorder columns

tax.bin=tax_table(top40)
tax.bin=as.data.frame(tax.bin@.Data)
ann_plq=tax.bin[,2]
ann_plq=data.frame(ann_plq)

pheatmap(as.matrix(OTUdf.ord), border_color = "black",fontsize_row=8, 
         fontsize_col = 8, scale = "none",
         cluster_cols = F, cluster_rows = T,
         labels_col = NULL, annotation_col = ann,
         labels_row=paste0(ann_plq$ann_plq))

