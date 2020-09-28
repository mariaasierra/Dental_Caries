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


#PCOA saliva vs plaque

GP.ord <- ordinate(sdf_phylo, "PCoA", "Unifrac")
# New facet label names for Description variable
Description2<- c("caries-free" = "Caries Free", 
                 "caries-active" = "Caries Active",
                 "non-responders" = "Non responders", 
                 "SDF"="SDF")
#p value
dist <- phyloseq::distance(sdf_phylo, method = "Unifrac")
perma <- adonis(dist~sample.type, data = as(sample_data(sdf_phylo), "data.frame"), permutations = 1000)
perma

pcoa=plot_ordination(sdf_phylo, GP.ord,  color="sample.type", axes =1:2)+ 
  stat_ellipse(type = "t",linetype = 1,alpha=0.5)+theme_minimal(base_size = 30, base_line_size = 0.5)+ geom_point(size=5)+
  scale_color_brewer(palette="Set1",labels=c("Plaque", "Saliva")) +labs(color="Sample type")+
  annotate("text", x = -0.3, y = 0.5, 
           label = "paste(italic(p), \" = 0.00099\")", parse = TRUE, size=8)
pcoa

#PCoA saliva vs plaque By teatment

pcoa2=plot_ordination(sdf_phylo, GP.ord, color="sample.type", axes =1:2)+ 
  stat_ellipse(type = "t",linetype = 1,alpha=0.5)
pcoa3=pcoa2 + facet_wrap( ~ treatment, labeller = as_labeller(Description2)) + theme_minimal(base_size = 20, base_line_size = 0.5)+ geom_point(size=5)+
  scale_color_brewer(palette="Set1",labels=c("Plaque", "Saliva")) +labs(color="Sample type")

#p values

caries_act = subset_samples(sdf_phylo, treatment == "caries-active") #Subset TypePlaque
dist <- phyloseq::distance(caries_act, method = "bray")
adonis(dist~sample.type, data = as(sample_data(caries_act), "data.frame"), permutations = 1000)

caries_free = subset_samples(sdf_phylo, treatment == "caries-free") #Subset TypePlaque
dist <- phyloseq::distance(caries_free, method = "bray")
adonis(dist~sample.type, data = as(sample_data(caries_free), "data.frame"), permutations = 1000)

non_res = subset_samples(sdf_phylo, treatment == "non-responders") #Subset TypePlaque
dist <- phyloseq::distance(non_res, method = "bray")
adonis(dist~sample.type, data = as(sample_data(non_res), "data.frame"), permutations = 1000)

sdf = subset_samples(sdf_phylo, treatment == "SDF") #Subset TypePlaque
dist <- phyloseq::distance(sdf, method = "bray")
adonis(dist~sample.type, data = as(sample_data(sdf), "data.frame"), permutations = 1000)


dat_text <- data.frame(
  label = c("p = 0.0009", "p = 0.0019", "p = 0.0009", "p = 0.015"),
  treatment   = c("caries-active", "caries-free", "non-responders", "SDF"),
  Axis.1    = c(-.5, -.5, -.5,-.5),
  Axis.2    = c(.7, .7, .7,.7)
)

pcoa4=pcoa3 + geom_text(data = dat_text, aes(x = Axis.1,  y = Axis.2, label = label), color="black", size=5) 
pcoa4

#PCoA Saliva four treatments 
GP.ord <- ordinate(slv_phylo, "PCoA", "Unifrac")
pcoa_slv=plot_ordination(slv_phylo, GP.ord, color="treatment", axes =1:2)+ 
  stat_ellipse(type = "t",linetype = 1,alpha=0.5)
pcoa_slv + theme_minimal(base_size = 20, base_line_size = 0.5)+ geom_point(size=5)+
  scale_color_brewer(palette="Dark2") +labs(color="Sample type") + ggtitle("Saliva")

#PCoA Plaque 
#Infected
GP.ord <- ordinate(plq_inf_phylo, "PCoA", "Unifrac")
pcoa_infec=plot_ordination(plq_inf_phylo, GP.ord, color="treatment", axes =1:2)+ 
  stat_ellipse(type = "t",linetype = 1,alpha=0.5)
pcoa_infec + theme_minimal(base_size = 20, base_line_size = 0.5)+ geom_point(size=5)+
  scale_color_brewer(palette="Dark2") +labs(color="Sample type") + ggtitle("Plaque Infected")

GP.ord <- ordinate(plq_con_phylo, "PCoA", "Unifrac")
pcoa_con=plot_ordination(plq_con_phylo, GP.ord, color="treatment", axes =1:2)+ 
  stat_ellipse(type = "t",linetype = 1,alpha=0.5)
pdf("all_plq_contra_by_description_pcoa.pdf", height = 8, width = 8)

pcoa_con + theme_minimal(base_size = 20, base_line_size = 0.5)+ geom_point(size=5)+
  scale_color_brewer(palette="Dark2") +labs(color="Sample type") + ggtitle("Plaque Contralateral")


#Plaque comparisons
GP.ord <- ordinate(plq_phylo, "PCoA", "bray")
plot_ordination(plq_phylo, GP.ord,  color="type.plaque", axes =1:2)+ 
  stat_ellipse(type = "t",linetype = 1,alpha=0.5)+ 
  facet_wrap( ~ treatment,labeller = as_labeller(Description2)) + 
  theme_minimal(base_size = 20, base_line_size = 0.5)+ geom_point(size=5)+labs(color="Plaque type") + scale_color_manual(values=c("olivedrab3", "red3"), labels=c("Contralteral", "Infected"))

caries_act = subset_samples(plq_phylo, treatment == "caries-active") #Subset TypePlaque
dist <- phyloseq::distance(caries_act, method = "unifrac")
adonis(dist~type.plaque, data = as(sample_data(caries_act), "data.frame"), permutations = 1000)

caries_free = subset_samples(plq_phylo, treatment == "caries-free") #Subset TypePlaque
dist <- phyloseq::distance(caries_free, method = "unifrac")
adonis(dist~type.plaque, data = as(sample_data(caries_free), "data.frame"), permutations = 1000)

non_res = subset_samples(plq_phylo, treatment == "non-responders") #Subset TypePlaque
dist <- phyloseq::distance(non_res, method = "unifrac")
adonis(dist~type.plaque, data = as(sample_data(non_res), "data.frame"), permutations = 1000)

sdf = subset_samples(plq_phylo, treatment == "SDF") #Subset TypePlaque
dist <- phyloseq::distance(sdf, method = "unifrac")
adonis(dist~type.plaque, data = as(sample_data(sdf), "data.frame"), permutations = 1000)


#Plaque caries active vs caries free
plq_sdf_phylo = subset_samples(plq_phylo, treatment == "caries-active" | treatment == "caries-free"  ) #Subset TypePlaque
sample_data(plq_sdf_phylo)

GP.ord <- ordinate(plq_sdf_phylo, "PCoA", "unifrac")
plot_ordination(plq_sdf_phylo, GP.ord,  color="treatment", axes =1:2)+ 
  stat_ellipse(type = "t",linetype = 1,alpha=0.5)+
  theme_minimal(base_size = 20, base_line_size = 0.5)+ geom_point(size=5) +labs(color="Treatment")+ scale_color_manual(values=c("olivedrab3", "red3"), labels=c("Caries Active", "Caries Free"))+
  annotate("text", x = -0.4, y = 0.5, 
           label = "paste(italic(p), \" = 0.003\")", parse = TRUE, size=8)

dist <- phyloseq::distance(plq_sdf_phylo, method = "unifrac")
adonis(dist~treatment, data = as(sample_data(plq_sdf_phylo), "data.frame"), permutations = 1000)


#Plaque non-responders vs sdf
plq_sdf_phylo = subset_samples(plq_phylo, treatment == "non-responders" | treatment == "SDF"  ) #Subset TypePlaque
sample_data(plq_sdf_phylo)

GP.ord <- ordinate(plq_sdf_phylo, "PCoA", "unifrac")
plot_ordination(plq_sdf_phylo, GP.ord,  color="treatment", axes =1:2)+ 
  stat_ellipse(type = "t",linetype = 1,alpha=0.5)+
  theme_minimal(base_size = 20, base_line_size = 0.5)+ geom_point(size=5)+
  labs(color="Treatment")+ scale_color_manual(values=c("olivedrab3", "red3"), labels=c("Non-responders", "SDF"))

dist <- phyloseq::distance(plq_sdf_phylo, method = "bray")
adonis(dist~treatment, data = as(sample_data(plq_sdf_phylo), "data.frame"), permutations = 1000)





