#Plaque SDF
### Prevalence Total Caries data ###
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
TAX=read.table("taxonomy.tsv", header = T,sep="\t")
row.names(TAX) = TAX$Feature.ID
TAX=TAX[,-1]
tax.clean = separate(TAX, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
tax.clean = tax.clean[,-which(names(tax.clean) %in% c("Confidence", "OTU"))]
OTU.UF = otu_table(OTUs, taxa_are_rows=F)
tax.UF = tax_table(as.matrix(tax.clean))
meta = read.table("metadata.txt", header=TRUE, row.names=1, sep="\t", dec = ".") #Sample names MUST be the same in both files, otherwise is going not going to pair data 
meta<-sample_data(meta)
tree <- read_tree ("tree.nwk")
physeq = phyloseq(OTU.UF,tax.UF, meta, tree)
sdf_phylo <- merge_phyloseq (OTU.UF,tax.UF, meta, tree)

#Alpha diversity 
plot_richness=function (physeq, x = "samples", color = NULL, shape = NULL, fill=NULL,
                        title = NULL, scales = "free_y", nrow = 1, shsi = NULL, measures = NULL, 
                        sortby = NULL) 
{
  erDF = estimate_richness(physeq, split = TRUE, measures = measures)
  measures = colnames(erDF)
  ses = colnames(erDF)[grep("^se\\.", colnames(erDF))]
  measures = measures[!measures %in% ses]
  if (!is.null(sample_data(physeq, errorIfNULL = FALSE))) {
    DF <- data.frame(erDF, sample_data(physeq))
  }
  else {
    DF <- data.frame(erDF)
  }
  if (!"samples" %in% colnames(DF)) {
    DF$samples <- sample_names(physeq)
  }
  if (!is.null(x)) {
    if (x %in% c("sample", "samples", "sample_names", "sample.names")) {
      x <- "samples"
    }
  }
  else {
    x <- "samples"
  }
  mdf = reshape2::melt(DF, measure.vars = measures)
  mdf$se <- NA_integer_
  if (length(ses) > 0) {
    selabs = ses
    names(selabs) <- substr(selabs, 4, 100)
    substr(names(selabs), 1, 1) <- toupper(substr(names(selabs), 
                                                  1, 1))
    mdf$wse <- sapply(as.character(mdf$variable), function(i, 
                                                           selabs) {
      selabs[i]
    }, selabs)
    for (i in 1:nrow(mdf)) {
      if (!is.na(mdf[i, "wse"])) {
        mdf[i, "se"] <- mdf[i, (mdf[i, "wse"])]
      }
    }
    mdf <- mdf[, -which(colnames(mdf) %in% c(selabs, "wse"))]
  }
  if (!is.null(measures)) {
    if (any(measures %in% as.character(mdf$variable))) {
      mdf <- mdf[as.character(mdf$variable) %in% measures, 
                 ]
    }
    else {
      warning("Argument to `measures` not supported. All alpha-diversity measures (should be) included in plot.")
    }
  }
  if (!is.null(shsi)) {
    warning("shsi no longer supported option in plot_richness. Please use `measures` instead")
  }
  if (!is.null(sortby)) {
    if (!all(sortby %in% levels(mdf$variable))) {
      warning("`sortby` argument not among `measures`. Ignored.")
    }
    if (!is.discrete(mdf[, x])) {
      warning("`sortby` argument provided, but `x` not a discrete variable. `sortby` is ignored.")
    }
    if (all(sortby %in% levels(mdf$variable)) & is.discrete(mdf[, 
                                                                x])) {
      wh.sortby = which(mdf$variable %in% sortby)
      mdf[, x] <- factor(mdf[, x], levels = names(sort(tapply(X = mdf[wh.sortby, 
                                                                      "value"], INDEX = mdf[wh.sortby, x], mean, na.rm = TRUE, 
                                                              simplify = TRUE))))
    }
  }
  richness_map = aes_string(x = x, y = "value", colour = color, fill= fill,
                            shape = shape)
  p = ggplot(mdf, richness_map) + 
    if (any(!is.na(mdf[, "se"]))) {
      p = p + geom_errorbar(aes(ymax = value + se, ymin = value - 
                                  se), width = 0.1)
    }
  p = p + theme(axis.text.x = element_text(angle = -90, vjust = 0.5, 
                                           hjust = 0))
  p = p + ylab("Alpha Diversity Measure")
  p = p + facet_wrap(~variable, nrow = nrow, scales = scales)
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}

#Trim saliva samples

slv_phylo = subset_samples(sdf_phylo, sample.type == "Saliva") #Subset category
metadata_slv<-data.frame(sample_data(slv_phylo))
pr=plot_richness(slv_phylo, "treatment", color= "treatment", fill="treatment",
                 measures=c("Observed", "Shannon", "Simpson"))+
  geom_boxplot(color="gray40",size=0.3)+ scale_fill_brewer(palette="Blues")

pr2=pr+  theme(axis.text.x = element_text(size=15, 
                                          color="black"),
               legend.position='none', 
               axis.title.y = element_text(size = 15, color="black"),
               axis.title.x =element_blank(),
               strip.text.x = element_text(size = 15),
               panel.background = element_rect(fill = "white",colour = "gray50"),
               panel.grid.major = element_line(size = 0.1, 
                                               linetype = 'solid',
                                               colour = "gray70"), 
               panel.grid.minor = element_line(size = 0.1, 
                                               linetype = 'solid',
                                               colour = "gray70"))+ggtitle("Saliva")+
  scale_x_discrete(labels=c("caries_free" = "Caries Free", 
                            "caries_active" = "Caries Active",
                            "non_responder" = "Non responders", 
                            "sdf"="SDF"))

#Calculate p values saliva

rich=estimate_richness(slv_phylo, measures=c("Observed", "Simpson", "Shannon"))
rich$SampleID=rownames(rich)
metadata=data.frame(sample_data(slv_phylo))
metadata$SampleID=row.names(metadata)
joined_slv=merge(rich, metadata, by = "SampleID")

#shannon
TukeyHSD(aov(Shannon ~ treatment, data=joined_slv))
#simpson
TukeyHSD(aov(Simpson ~ treatment, data=joined_slv))
#observed
TukeyHSD(aov(Observed ~ treatment, data=joined_slv))

#Trim plaque samples
plq_phylo = subset_samples(sdf_phylo, sample.type == "Plaque") #Subset category
sample_data(plq_phylo)
plq_con_phylo = subset_samples(plq_phylo, type.plaque == "contra") #Subset TypePlaque
plq_inf_phylo= subset_samples(plq_phylo, type.plaque == "infected") #Subset TypePlaque

#Contralateral Site
con=plot_richness(plq_con_phylo, "treatment", color= "treatment", fill="treatment",
                  measures=c("Observed", "Shannon", "Simpson"))+
  geom_boxplot(color="gray40",size=0.3)+ scale_fill_brewer(palette="Oranges")

con2=con+  theme(axis.text.x = element_text(size=15, 
                                            color="black"),
                 legend.position='none', 
                 axis.title.y = element_text(size = 15, color="black"),
                 axis.title.x =element_blank(),
                 strip.text.x = element_text(size = 15),
                 panel.background = element_rect(fill = "white",colour = "gray50"),
                 panel.grid.major = element_line(size = 0.1, 
                                                 linetype = 'solid',
                                                 colour = "gray70"), 
                 panel.grid.minor = element_line(size = 0.1, 
                                                 linetype = 'solid',
                                                 colour = "gray70"))+ggtitle("Plaque Contralateral")+
  scale_x_discrete(labels=c("caries-free" = "Caries Free", 
                            "caries-active" = "Caries Active",
                            "non-responders" = "Non responders", 
                            "SDF"="SDF")) + geom_signif(comparisons =list( c("SDF", "non-responders")), map_signif_level = T, textsize = 5, color="gray30")
con2

#Calculate p values plaque contralateral
metadata=data.frame(sample_data(plq_con_phylo))
rich=estimate_richness(plq_con_phylo, measures=c("Observed", "Simpson", "Shannon"))
rich$SampleID=rownames(rich)
metadata$SampleID=row.names(metadata)
joined_plq_con=merge(rich, metadata, by = "SampleID")

#shannon
TukeyHSD(aov(Shannon ~ treatment, data=joined_plq_con))
#simpson
TukeyHSD(aov(Simpson ~ treatment, data=joined_plq_con))
#observed
TukeyHSD(aov(Observed ~ treatment, data=joined_plq_con))

#Infected Site
inf=plot_richness(plq_inf_phylo, "treatment", color= "treatment", fill="treatment",
                  measures=c("Observed", "Shannon", "Simpson"))+
  geom_boxplot(color="gray40",size=0.3)+ scale_fill_brewer(palette="Reds")

inf2=inf+  theme(axis.text.x = element_text(size=15, 
                                            color="black"),
                 legend.position='none', 
                 axis.title.y = element_text(size = 15, color="black"),
                 axis.title.x =element_blank(),
                 strip.text.x = element_text(size = 15),
                 panel.background = element_rect(fill = "white",colour = "gray50"),
                 panel.grid.major = element_line(size = 0.1, 
                                                 linetype = 'solid',
                                                 colour = "gray70"), 
                 panel.grid.minor = element_line(size = 0.1, 
                                                 linetype = 'solid',
                                                 colour = "gray70"))+ggtitle("Plaque Infected")+
  scale_x_discrete(labels=c("caries-free" = "Caries Free", 
                            "caries-active" = "Caries Active",
                            "non-responders" = "Non responders", 
                            "SDF"="SDF"))

#Calculate p values plaque infected
metadata=data.frame(sample_data(plq_inf_phylo))
rich=estimate_richness(plq_inf_phylo, measures=c("Observed", "Simpson", "Shannon"))
rich$SampleID=rownames(rich)
metadata$SampleID=row.names(metadata)
joined_plq_inf=merge(rich, metadata, by = "SampleID")

#shannon
TukeyHSD(aov(Shannon ~ treatment, data=joined_plq_inf))
#simpson
TukeyHSD(aov(Simpson ~ treatment, data=joined_plq_inf))
#observe
TukeyHSD(aov(Observed ~ treatment, data=joined_plq_inf))


#COMPARISONS BETWEEN SAME PATIENTS DIFFERENT TOOTH.

#Caries active
plq_phylo = subset_samples(sdf_phylo, sample.type == "Plaque") #Subset category
sample_data(plq_phylo)
plq_act_phylo = subset_samples(plq_phylo, treatment == "caries-active") #Subset TypePlaque

con=plot_richness(plq_act_phylo, "type.plaque", color= "type.plaque", fill="type.plaque",
                  measures=c("Observed", "Shannon", "Simpson"))+
  geom_boxplot(color="gray40",size=0.3)+ scale_fill_manual(values=c("olivedrab3", "red3"))

con2=con+  theme(axis.text.x = element_text(size=15, 
                                            color="black"),
                 legend.position='none', 
                 axis.title.y = element_text(size = 15, color="black"),
                 axis.title.x =element_blank(),
                 strip.text.x = element_text(size = 15),
                 panel.background = element_rect(fill = "white",colour = "gray50"),
                 panel.grid.major = element_line(size = 0.1, 
                                                 linetype = 'solid',
                                                 colour = "gray70"), 
                 panel.grid.minor = element_line(size = 0.1, 
                                                 linetype = 'solid',
                                                 colour = "gray70"))+ggtitle("Caries Active")+
  scale_x_discrete(labels=c("contra" = "Contralateral", 
                            "infected" = "Infected"))
#Calculate p values caries active
metadata=data.frame(sample_data(plq_act_phylo))
rich=estimate_richness(plq_act_phylo, measures=c("Observed", "Simpson", "Shannon"))
rich$SampleID=rownames(rich)
metadata$SampleID=row.names(metadata)
joined_plq_inf=merge(rich, metadata, by = "SampleID")

#shannon
t.test(Shannon~ type.plaque,data=joined_plq_inf)
#simpson
t.test(Simpson~ type.plaque,data=joined_plq_inf)
#observe
t.test(Observed~ type.plaque,data=joined_plq_inf)

#Non-responders
plq_phylo = subset_samples(sdf_phylo, sample.type == "Plaque") #Subset category
sample_data(plq_phylo)
plq_non_phylo = subset_samples(plq_phylo, treatment == "non-responders") #Subset TypePlaque

con=plot_richness(plq_non_phylo, "type.plaque", color= "type.plaque", fill="type.plaque",
                  measures=c("Observed", "Shannon", "Simpson"))+
  geom_boxplot(color="gray40",size=0.3)+ scale_fill_manual(values=c("olivedrab3", "red3"))

con2=con+  theme(axis.text.x = element_text(size=15, 
                                            color="black"),
                 legend.position='none', 
                 axis.title.y = element_text(size = 15, color="black"),
                 axis.title.x =element_blank(),
                 strip.text.x = element_text(size = 15),
                 panel.background = element_rect(fill = "white",colour = "gray50"),
                 panel.grid.major = element_line(size = 0.1, 
                                                 linetype = 'solid',
                                                 colour = "gray70"), 
                 panel.grid.minor = element_line(size = 0.1, 
                                                 linetype = 'solid',
                                                 colour = "gray70"))+ggtitle("Non-responders")+
  scale_x_discrete(labels=c("contra" = "Contralateral", 
                            "infected" = "Infected")) +geom_signif(comparisons = list(c("contra", "infected")), map_signif_level = T, color='gray30')

#Calculate p values non responders
metadata=data.frame(sample_data(plq_non_phylo))
rich=estimate_richness(plq_non_phylo, measures=c("Observed", "Simpson", "Shannon"))
rich$SampleID=rownames(rich)
metadata$SampleID=row.names(metadata)
joined_plq_inf=merge(rich, metadata, by = "SampleID")

#shannon
t.test(Shannon~ type.plaque,data=joined_plq_inf)
#simpson
t.test(Simpson~ type.plaque,data=joined_plq_inf)
#observe
t.test(Observed~ type.plaque,data=joined_plq_inf)

#SDF
plq_phylo = subset_samples(sdf_phylo, sample.type == "Plaque") #Subset category
sample_data(plq_phylo)
plq_sdf_phylo = subset_samples(plq_phylo, treatment == "SDF") #Subset TypePlaque

con=plot_richness(plq_sdf_phylo, "type.plaque", color= "type.plaque", fill="type.plaque",
                  measures=c("Observed", "Shannon", "Simpson"))+
  geom_boxplot(color="gray40",size=0.3)+ scale_fill_manual(values=c("olivedrab3", "red3"))

con2=con+  theme(axis.text.x = element_text(size=15, 
                                            color="black"),
                 axis.title.y = element_text(size = 15, color="black"),
                 axis.title.x =element_blank(),
                 strip.text.x = element_text(size = 15),
                 panel.background = element_rect(fill = "white",colour = "gray50"),
                 panel.grid.major = element_line(size = 0.1, 
                                                 linetype = 'solid',
                                                 colour = "gray70"), 
                 panel.grid.minor = element_line(size = 0.1, 
                                                 linetype = 'solid',
                                                 colour = "gray70"))+ggtitle("SDF")+
  scale_x_discrete(labels=c("contra" = "Contralateral", 
                            "infected" = "Infected"))

#Calculate p values non responders
metadata=data.frame(sample_data(plq_sdf_phylo))
rich=estimate_richness(plq_sdf_phylo, measures=c("Observed", "Simpson", "Shannon"))
rich$SampleID=rownames(rich)
metadata$SampleID=row.names(metadata)
joined_plq_inf=merge(rich, metadata, by = "SampleID")

#shannon
t.test(Shannon~ type.plaque,data=joined_plq_inf)
#simpson
t.test(Simpson~ type.plaque,data=joined_plq_inf)
#observe
t.test(Observed~ type.plaque,data=joined_plq_inf)






