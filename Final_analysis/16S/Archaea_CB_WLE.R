###########################################################
#########Bacterial community composition combined WLE and CB
##########################################################
##########using de-novo clustering at 99% identity #######
##########################################################
##########and using the new extractions done manually for WLE
##########################################################


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(grid)
library(reshape2)

#install.packages("vctrs")
#update.packages("tidyverse")
#update.packages("scales")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

library(decontam)

otu = read.csv("otu_table_dn_99.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("sample-metadata-cb-wle-decontam.csv", sep=",", row.names=1)
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 48637 taxa and 215 samples ]
sample_data() Sample Data:       [ 215 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 48637 taxa by 8 taxonomic ranks ]

sample_names(phyloseq_merged)

##check column names of taxonomy table
colnames(tax_table(phyloseq_merged))

##remove mito, chloroplast and archaea, eukarya
phyloseq_merged_clean <- phyloseq_merged %>%
  subset_taxa(
    Kingdom == "d__Archaea" &
      Family   != "f__Chloroplast" &
      Family  != "f__Mitochondria"
  )
phyloseq_merged_clean


##bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 46342 taxa and 215 samples ]
#sample_data() Sample Data:       [ 215 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 46342 taxa by 8 taxonomic ranks ]


#archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1454 taxa and 215 samples ]
#sample_data() Sample Data:       [ 215 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 1454 taxa by 8 taxonomic ranks ]

head(sample_data(phyloseq_merged_clean))

##inspect library sizes
df <- as.data.frame(sample_data(phyloseq_merged_clean)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phyloseq_merged_clean)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
plot<- ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point() ##Supplemental Figure S3

ggsave(filename = "librarySize_archaea.pdf", plot = plot,
       width = 16,
       height = 13, units = c("cm"),
       dpi = 300)


##identify contaminants by prevalence. in this method, the distribution 
#of the frequency of each sequence feature as a function of the 
#prevalence is used to identify contaminants. In this method, 
#the prevalence (presence/absence across samples) of each sequence feature 
#in true positive samples is compared to the prevalence in negative controls 
#to identify contaminants.

sample_data(phyloseq_merged_clean)$is.neg <- sample_data(phyloseq_merged_clean)$Sample_or_Control == "Control Sample"
contamdf.prev0.5 <- isContaminant(phyloseq_merged_clean, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev0.5$contaminant)

##bacteria
##FALSE  TRUE 
# 46192   150

##archaea
#FALSE  TRUE 
#1446     8


# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(phyloseq_merged_clean, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev0.5$contaminant)
plot<- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") ##Supplemental Figure S3

ggsave(filename = "prevalence_archaea.pdf", plot = plot,
       width = 16,
       height = 13, units = c("cm"),
       dpi = 300)

write.csv(df.pa, "contaminant-table-cb-wle-archaea.csv")

write.csv(contamdf.prev0.5, "contaminant-prev-0.5-cb-wle-archaea.csv")

##removing contaminants from phyloseq object
df.pa <- read.csv("contaminant-prev-0.5-cb-wle-archaea.csv")
View(df.pa)
subset.df <- subset(df.pa, contaminant== "FALSE")
View(subset.df)
keep.taxa <- as.vector(subset.df$X)
keep.taxa

phyloseq_merged_clean_decontam <- prune_taxa(keep.taxa, phyloseq_merged_clean)
phyloseq_merged_clean_decontam

##bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 46192 taxa and 215 samples ]
#sample_data() Sample Data:       [ 215 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 46192 taxa by 8 taxonomic ranks ]

##archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1446 taxa and 215 samples ]
#sample_data() Sample Data:       [ 215 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 1446 taxa by 8 taxonomic ranks ]

##export otu table out of phyloseq object 

OTU1 = as(otu_table(phyloseq_merged_clean_decontam), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

write.csv(OTUdf, "OTU_clean_archaea.csv")

#remove neg controls from otu table-- might add more rows which are rowsums=0, read back in to rarefy reads to 15K

otu = read.csv("OTU_clean_archaea.csv", sep=",", row.names=1)#with no contaminants, mito, chloroplast, neg control 
tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("sample-metadata-cb-wle-decontam.csv", sep=",", row.names=1)

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

##bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 46192 taxa and 205 samples ]
#sample_data() Sample Data:       [ 205 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 46192 taxa by 8 taxonomic ranks ]

#archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1446 taxa and 205 samples ]
#sample_data() Sample Data:       [ 205 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 1446 taxa by 8 taxonomic ranks ]


#############################################
#####rarefaction curves-- all samples#######
#############################################


otu1 = read.csv("OTU_clean_archaea.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
tax1 = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax1 = as.matrix(tax1)
metadata1 = read.csv("sample-metadata-cb-wle-decontam.csv", sep=",", row.names=1)
OTU1 = otu_table(otu1, taxa_are_rows = TRUE)
TAX1 = tax_table(tax1)
meta1 = sample_data(metadata1)

phyloseq_merged1 = phyloseq(OTU1, TAX1, meta1)
phyloseq_merged1


#BACTERIA
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 46192 taxa and 205 samples ]
#sample_data() Sample Data:       [ 205 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 46192 taxa by 8 taxonomic ranks ]

#archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1446 taxa and 205 samples ]
#sample_data() Sample Data:       [ 205 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 1446 taxa by 8 taxonomic ranks ]



#col<- cyan4
library(vegan)

mat<- as(t(otu_table(phyloseq_merged1)), "matrix")
system.time(rarecurve(mat, step=1000, col="cyan4", label=FALSE)) 
#rarecurve<-rarecurve(t(otu_table(phyloseq_merged1)), step=20, label = FALSE, col="cyan4") ## Supplemental Figure S3

###############################
sample_df <- data.frame(sample_data(phyloseq_merged))
View(sample_df)
write.csv(sample_df, "samples_pre_rarefied_archaea.csv")

library(vegan)
#rarecurve(t(otu_table(phyloseq_merged)), step=50)
phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged, sample.size = 2000, rngseed = TRUE, trimOTUs=FALSE)

##bacteria
#6 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 

#s20230217078s20230901075s20230901078s20230901081s20230901082

##archaea
#set.seed(TRUE)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(TRUE); .Random.seed` for the full vector
...
#131 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
  
  #s20230217001s20230217002s20230217003s20230217004s20230217005
...


sample_df <- data.frame(sample_data(phyloseq_rarefy))

View(sample_df)
write.csv(sample_df, "samples_rarefied_archaea.csv")

##use anti-join to select non-present samples in rarefied object
#This join is like eg: for df1-df2, it selects all rows from df1 that are not present in df2.
df1<- read.csv("samples_pre_rarefied_archaea.csv")
df2 <- read.csv("samples_rarefied_archaea.csv")
df= df1 %>% anti_join(df2,by="X")

View(df)

write.csv(df, "samples-lost-rarefaction-archaea.csv")

OTU2 = as(otu_table(phyloseq_rarefy), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, "OTU_clean_noneg_rarefied2k_archaea.csv")

########################
######## barplots and ordination plots
#########################



otu_final<- read.csv("OTU_clean_noneg_rarefied2k_archaea.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #467 otus with no counts for archaea
View(otu1)
max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
metadata <- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged


##bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 40216 taxa and 199 samples ]
#sample_data() Sample Data:       [ 199 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]

##archaea
#otu_table()   OTU Table:         [ 979 taxa and 74 samples ]
#sample_data() Sample Data:       [ 74 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 979 taxa by 8 taxonomic ranks ]

##barplot class level
phyloseq_class <- phyloseq_merged %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Genus) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."

View(phyloseq_class)


palette_new56 = c("#560d0d", "#dba4a4", "#cc1c1c","#111b77","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","green4", "greenyellow", "#b7ffdb","#825121",
                  "cyan",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#3f3e3d","#000000",
                  "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
                  "#599861", "#CBD588", "#508578","#FFD700", "coral",
                  "lightslateblue","green", "#AD6F3B","#FFA500", "#CD9BCD","darkslategrey","#a0fffc", "#DA5724","#fff0d6","#8A7C64",
                  "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4",
                  "#fa7efc","lightgoldenrod4","#D14285", "darkorchid1", "deepskyblue4","magenta",
                  "forestgreen","goldenrod",
                  "lightpink3", "lightgreen",
                  "lightgoldenrod1")

install.packages("ggh4x")
library(ggh4x)
phyloseq_class_new<- phyloseq_class %>% mutate(transect1= ifelse(transect=="wetland-center"|transect=="wetland-transition-edge", "wetland", transect))

phyloseq_class_new$transect1 <- factor(phyloseq_class_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)
View(phyloseq_class_new)


##plot by region/reansect/horizon/site
plot_class<- ggplot(phyloseq_class_new, aes(x= tree_number, y = Abundance, fill =Genus)) + 
  facet_nested(~region+site+transect1+horizon, scales="free", independent="x", labeller = labeller(horizon = function(x) {rep("", length(x))}))  + theme(strip.text.x = element_text(size = 20, face="bold"))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(palette_new56))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=20, face="bold"),
                                     axis.title.y=element_text(size=20, face="bold")) + 
  theme(legend.title = element_text(size=20, face="bold"))+ theme(legend.text = element_text(size=17))+ guides(fill=guide_legend(ncol=1))
plot_class

##Supplementary Figure
ggsave(filename = "barplot_genus_cb_wle_new_archaea_nohorizon_07.08.24.jpg", plot = plot_class,
       width = 80,
       height = 40, units = c("cm"),
       dpi = 300)



####ordination plots

##plot by region

metadata_new<- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)
sample_data_new<- metadata_new %>% mutate(transect1= ifelse(transect=="wetland-center"|transect=="wetland-transition-edge", "wetland", transect))

sample_data_new$transect1 <- factor(sample_data_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)

sample_data_phyloseq<- sample_data(sample_data_new)
phyloseq_merged_new<- merge_phyloseq(otu, tax, sample_data_phyloseq)
phyloseq_merged_new

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 979 taxa and 74 samples ]
#sample_data() Sample Data:       [ 74 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 979 taxa by 8 taxonomic ranks ]

phyloseq_pcoa <- ordinate(
  physeq = phyloseq_merged_new, ##considering using phyloseq_rarefy directly
  method = "PCoA", 
  distance = "bray"
)


pcoa_cb_wle<-plot_ordination(
  physeq = phyloseq_merged_new,
  ordination = phyloseq_pcoa,
  color = "region",
  shape = "transect1")+
  #title = "PCoA WLE CB archaea") + 
  #scale_color_manual(values = c("#CC79A7", "#E69F00", "#0072B2", "black"))+
  scale_color_manual(values=c("#009E73", "#D55E00"))+
  scale_shape_manual(values=c(15,16,17,18))+
  geom_point(aes(color = region), alpha = 0.7, size = 4) +
  geom_point(colour = "white", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))+labs(shape="transect")+theme_classic()

pcoa_cb_wle

ggsave(filename="pcoa-cb-wle-archaea-07.05.24.TIFF", plot=pcoa_cb_wle, width=8, height=6, units="in", dpi=300)


#pcoa separated by region

##cb
View(sample_data_new)
sample_data_new_cb<- sample_data_new %>% filter(region=="Chesapeake Bay")

sample_data_phyloseq_cb<- sample_data(sample_data_new_cb)
phyloseq_merged_new<- merge_phyloseq(otu, tax, sample_data_phyloseq_cb)
phyloseq_merged_new

##bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 40216 taxa and 117 samples ]
#sample_data() Sample Data:       [ 117 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]

##archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 979 taxa and 22 samples ]
#sample_data() Sample Data:       [ 22 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 979 taxa by 8 taxonomic ranks ]

##removing taxa with rows==0, testing the ordination with and without missing taxa

OTU2 = as(otu_table(phyloseq_merged_new), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, "OTU_CB_rarefied2k_archaea.csv")

##read back in new OTU csv for CB, remove rowsums and merge into phyloseq

otu_final<- read.csv("OTU_CB_rarefied2k_archaea.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #475 otus with no counts, this makes total otus less for cb than before because we have clustered otus at 99%
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)##504 otus

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)

metadata_new<- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)
sample_data_new<- metadata_new %>% mutate(transect1= ifelse(transect=="wetland-center", "wetland", transect))

sample_data_new$transect1 <- factor(sample_data_new$transect1,levels = c("upland", "transition", "wetland", "wetland-transition-edge"), ordered=TRUE)

View(sample_data_new)
sample_data_new_cb<- sample_data_new %>% filter(region=="Chesapeake Bay")

sample_data_phyloseq_cb<- sample_data(sample_data_new_cb)
tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)


colnames(tax)
tax
phyloseq_merged_new<- merge_phyloseq(otu, tax, sample_data_phyloseq_cb)
phyloseq_merged_new

##bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 24964 taxa and 117 samples ]
#sample_data() Sample Data:       [ 117 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 24964 taxa by 8 taxonomic ranks ]

##archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 504 taxa and 22 samples ]
#sample_data() Sample Data:       [ 22 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 504 taxa by 8 taxonomic ranks ]

  
  phyloseq_pcoa <- ordinate(
    physeq = phyloseq_merged_new, ##considering using phyloseq_rarefy directly
    method = "PCoA", 
    distance = "bray"
  )


pcoa_cb<-plot_ordination(
  physeq = phyloseq_merged_new,
  ordination = phyloseq_pcoa,
  color = "transect1",
  shape = "site")+
 # title = "PCoA CB archaea") + 
  scale_color_manual(values = c("#CC79A7", "#E69F00", "#0072B2"))+
  scale_shape_manual(values=c(15,16,17,18))+
  geom_point(aes(color = transect1), alpha = 0.7, size = 4) +
  geom_point(colour = "white", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))+labs(color="transect")+theme_classic()

pcoa_cb

##the plots dont change after removing rows of otus which are zero
#ggsave(filename="pcoa-cb-nozeroOTUs.TIFF", plot=pcoa_cb, width=8, height=6, units="in", dpi=300)
ggsave(filename="pcoa-cb-archaea-07.05.24.TIFF", plot=pcoa_cb, width=8, height=6, units="in", dpi=300)

##wle

otu_final<- read.csv("OTU_clean_noneg_rarefied2k_archaea.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #467 otus with no counts for archaea
View(otu1)
max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
metadata <- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged


metadata_new<- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)
sample_data_new<- metadata_new %>% mutate(transect1= ifelse(transect=="wetland-center", "wetland", transect))

sample_data_new$transect1 <- factor(sample_data_new$transect1,levels = c("upland", "transition", "wetland", "wetland-transition-edge"), ordered=TRUE)



View(sample_data_new)
sample_data_new_wle<- sample_data_new %>% filter(region=="Lake Erie ")##space added in csv file by mistake
View(sample_data_new_wle)
sample_data_phyloseq_wle<- sample_data(sample_data_new_wle)
phyloseq_merged_new<- merge_phyloseq(otu, tax, sample_data_phyloseq_wle)
phyloseq_merged_new

#bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 40216 taxa and 82 samples ]
#sample_data() Sample Data:       [ 82 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]

#archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 979 taxa and 52 samples ]
#sample_data() Sample Data:       [ 52 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 979 taxa by 8 taxonomic ranks ]


##removing taxa with rows==0, testing the ordination with and without missing taxa

OTU2 = as(otu_table(phyloseq_merged_new), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, "OTU_WLE_rarefied2k_archaea.csv")

##read back in new OTU csv for CB, remove rowsums and merge into phyloseq

otu_final<- read.csv("OTU_WLE_rarefied2k_archaea.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #484 otus with no counts
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)##495 otus, which is much higher than the 3000 ASVs we had before

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)

metadata_new<- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)
sample_data_new<- metadata_new %>% mutate(transect1= ifelse(transect=="wetland-center"|transect=="wetland-transition-edge", "wetland", transect))

sample_data_new$transect1 <- factor(sample_data_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)

View(sample_data_new)
sample_data_new_wle<- sample_data_new %>% filter(region=="Lake Erie ")

sample_data_phyloseq_wle<- sample_data(sample_data_new_wle)
tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)


colnames(tax)
tax
phyloseq_merged_new<- merge_phyloseq(otu, tax, sample_data_phyloseq_wle)
phyloseq_merged_new

#bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 19262 taxa and 82 samples ]
#sample_data() Sample Data:       [ 82 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 19262 taxa by 8 taxonomic ranks ]

#archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 495 taxa and 52 samples ]
#sample_data() Sample Data:       [ 52 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 495 taxa by 8 taxonomic ranks ]
  
  
  phyloseq_pcoa <- ordinate(
    physeq = phyloseq_merged_new, ##considering using phyloseq_rarefy directly
    method = "PCoA", 
    distance = "bray"
  )


pcoa_wle<-plot_ordination(
  physeq = phyloseq_merged_new,
  ordination = phyloseq_pcoa,
  color = "transect1",
  shape = "site")+
 # title = "PCoA WLE archaea") + 
  scale_color_manual(values = c("#CC79A7", "#E69F00", "#0072B2"))+
  scale_shape_manual(values=c(15,16,17,18))+
  geom_point(aes(color = transect1), alpha = 0.7, size = 4) +
  geom_point(colour = "white", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))+labs(color="transect")+theme_classic()

pcoa_wle

##exact same plot both times
#ggsave(filename="pcoa-wle-nozeroOTUs.TIFF", plot=pcoa_wle, width=8, height=6, units="in", dpi=300)
ggsave(filename="pcoa-wle-archaea-07.05.24.TIFF", plot=pcoa_wle, width=8, height=6, units="in", dpi=300)



####################
####PERMANOVA stats
##################

##by region
otu_final<- read.csv("OTU_clean_noneg_rarefied2k_archaea.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #5976 otus with no counts
View(otu1)
max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
metadata <- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

##bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 40216 taxa and 199 samples ]
#sample_data() Sample Data:       [ 199 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]

##archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 979 taxa and 74 samples ]
#sample_data() Sample Data:       [ 74 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 979 taxa by 8 taxonomic ranks ]


set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged, method = "bray") ##using merged object is fine here, gives same results as rarefy object
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged)) %>% mutate(transect1=ifelse(transect=="wetland-center"|transect=="wetland-transition-edge", "wetland", transect))
View(sample_df)
set.seed(1)
library(vegan)
adonis(phyloseq_bray ~ region*site*transect1, by="margin", data = sample_df)
adonis2(phyloseq_bray ~ region*site*transect1, by="terms", data = sample_df) #gives same result as above code

#CB AND WLE
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = phyloseq_bray ~ region * site * transect1, data = sample_df, by = "terms")
#Df SumOfSqs      R2       F Pr(>F)    
#region            1   6.6065 0.28217 78.4369  0.001 ***
  #site              4   4.4966 0.19206 13.3468  0.001 ***
 # transect1         2   3.7468 0.16003 22.2426  0.001 ***
  #region:transect1  2   2.7448 0.11724 16.2944  0.001 ***
#  site:transect1    3   0.6805 0.02906  2.6931  0.002 ** 
  #Residual         61   5.1378 0.21944                   
#Total            73  23.4131 1.00000                   

  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



##Chesapeake Bay
otu_final<- read.csv("OTU_CB_rarefied2k_archaea.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #15,252 otus with no counts, this makes total otus less for cb than before because we have clustered otus at 99%
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)##504 otus

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)

metadata_new<- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)
sample_data_new<- metadata_new %>% mutate(transect1= ifelse(transect=="wetland-center", "wetland", transect))

#sample_data_new$transect1 <- factor(sample_data_new$transect1,levels = c("upland", "transition", "wetland", "wetland-transition-edge"), ordered=TRUE)

View(sample_data_new)
sample_data_new_cb<- sample_data_new %>% filter(region=="Chesapeake Bay")

sample_data_phyloseq_cb<- sample_data(sample_data_new_cb)
tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)


colnames(tax)
tax
phyloseq_merged_new<- merge_phyloseq(otu, tax, sample_data_phyloseq_cb)
phyloseq_merged_new

#bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 24964 taxa and 117 samples ]
#sample_data() Sample Data:       [ 117 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 24964 taxa by 8 taxonomic ranks ]

#archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 504 taxa and 22 samples ]
#sample_data() Sample Data:       [ 22 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 504 taxa by 8 taxonomic ranks ]


set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_new, method = "bray") ##using merged object is fine here, gives same results as rarefy object
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged_new)) 
View(sample_df)
set.seed(1)
adonis2(phyloseq_bray ~ site*transect1, by="terms", data = sample_df)

## CB
#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = phyloseq_bray ~ site * transect1, data = sample_df, by = "terms")
#Df SumOfSqs      R2      F Pr(>F)    
#site       2   3.3187 0.42823 13.241  0.001 ***
 # transect1  2   2.3007 0.29686  9.179  0.001 ***
 # Residual  17   2.1305 0.27490                  
#Total     21   7.7499 1.00000                  

 # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



##Lake Erie

otu_final<- read.csv("OTU_WLE_rarefied2k_archaea.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #20954 otus with no counts
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)##495 otus

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)

metadata_new<- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)
sample_data_new<- metadata_new %>% mutate(transect1= ifelse(transect=="wetland-center"|transect=="wetland-transition-edge", "wetland", transect))

#sample_data_new$transect1 <- factor(sample_data_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)

View(sample_data_new)
sample_data_new_wle<- sample_data_new %>% filter(region=="Lake Erie ")

sample_data_phyloseq_wle<- sample_data(sample_data_new_wle)
tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)


colnames(tax)
tax
phyloseq_merged_new<- merge_phyloseq(otu, tax, sample_data_phyloseq_wle)
phyloseq_merged_new

#bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 19262 taxa and 82 samples ]
#sample_data() Sample Data:       [ 82 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 19262 taxa by 8 taxonomic ranks ]

#archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 495 taxa and 52 samples ]
#sample_data() Sample Data:       [ 52 samples by 16 sample variables ]
#tax_table()   Taxonomy Table:    [ 495 taxa by 8 taxonomic ranks ]


set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_new, method = "bray") 
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged_new)) 
View(sample_df)
set.seed(1)
library(vegan)
adonis2(phyloseq_bray ~ site*transect1, by="terms", data = sample_df)

#WLE

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#adonis2(formula = phyloseq_bray ~ site * transect1, data = sample_df, by = "terms")
#Df SumOfSqs      R2       F Pr(>F)    
#site            2   1.1779 0.13006  8.6168  0.001 ***
  #transect1       2   4.1910 0.46275 30.6590  0.001 ***
#  site:transect1  3   0.6805 0.07514  3.3187  0.003 ** 
#  Residual       44   3.0074 0.33206                   
#Total          51   9.0568 1.00000                   

#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##########################################
##load indicator species analysis package
##########################################

install.packages("indicspecies")
library(indicspecies)

###########

##phyloseq object
library(phyloseq)
library(tidyverse)
otu_final = read.csv("OTU_clean_noneg_rarefied2k_archaea.csv", sep=",", row.names=1) ##USING UNRAREFIED DATA, CAN ALSO CHECK WITH RAREFIED 

otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)



tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)

metadata = read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1) # use for crop indicators 


OTU = otu_table(otu_final, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

#bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 40216 taxa and 199 samples ]
#sample_data() Sample Data:       [ 199 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]

#archaea




##indicator species analysis

library(indicspecies)

GetIndicators <-function(dataframe, var){
  require(phyloseq); require(indicspecies); require(dplyr)
  #dataframe <- tax_glom(dataframe, taxrank = "Family") ##can be adjusted to indicate at family or class level, else default is OTU level
  otu <- as.data.frame(otu_table(dataframe, taxa_are_rows = TRUE))
  metadata = as(sample_data(dataframe), "data.frame")
  taxa <- as.data.frame(as.matrix(tax_table(dataframe)))
  multipatt <- multipatt(t(otu), metadata[,var], func = "r.g",
                         control=how(nperm=999), duleg=FALSE)
  multipatt -> multipatt_fdr
  multipatt_fdr$sign$p.value <- p.adjust(multipatt_fdr$sign$p.value, "fdr")
  #print(multipatt_fdr)
  multipatt_fdr$sign[which(
    multipatt_fdr$sign$p.value <= 0.05), ] -> indicator_taxa
  taxa$OTU <- rownames(taxa)
  data.frame(OTU = as.factor(row.names(indicator_taxa)), indicator_taxa) %>%
    dplyr::left_join(taxa, by="OTU") -> indicator_taxa
  rownames(indicator_taxa) <- indicator_taxa$OTU
  indicator_taxa <- arrange(indicator_taxa, desc(stat))
  return(indicator_taxa)
}



ind_region <-
  GetIndicators(dataframe=phyloseq_merged, "region")


##checking

View(ind_region)

ind_region_subset<- ind_region %>% subset(s.Chesapeake.Bay==1 & s.Lake.Erie.==1)
View(ind_region_subset) ##no common indicators across two regions, this hold true for this new analysis also


library(tidyverse)

write.csv(ind_region, "indicator_region_CB_WLE_archaea.csv")



###################################
###indicator for transect across sites in CB 
####################################

##load indicator species analysis package

install.packages("indicspecies")
library(indicspecies)

###########

##phyloseq object
library(phyloseq)
library(tidyverse)
otu = read.csv("OTU_CB_rarefied2k_archaea.csv", sep=",", row.names=1) ##USING RAREFIED DATA here unlike unrarefied data used for site specific indicators  

otu_final = otu[rowSums(otu[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) ##504 otus


tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)

metadata = read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1) 


OTU = otu_table(otu_final, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

#bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 24964 taxa and 117 samples ]
#sample_data() Sample Data:       [ 117 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 24964 taxa by 8 taxonomic ranks ]

#archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 504 taxa and 22 samples ]
#sample_data() Sample Data:       [ 22 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 504 taxa by 8 taxonomic ranks ]


##indicator species analysis

library(indicspecies)

GetIndicators <-function(dataframe, var){
  require(phyloseq); require(indicspecies); require(dplyr)
  #dataframe <- tax_glom(dataframe, taxrank = "Family") ##can be adjusted to indicate at family or class level, else default is OTU level
  otu <- as.data.frame(otu_table(dataframe, taxa_are_rows = TRUE))
  metadata = as(sample_data(dataframe), "data.frame")
  taxa <- as.data.frame(as.matrix(tax_table(dataframe)))
  multipatt <- multipatt(t(otu), metadata[,var], func = "r.g",
                         control=how(nperm=999), duleg=FALSE)
  multipatt -> multipatt_fdr
  multipatt_fdr$sign$p.value <- p.adjust(multipatt_fdr$sign$p.value, "fdr")
  #print(multipatt_fdr)
  multipatt_fdr$sign[which(
    multipatt_fdr$sign$p.value <= 0.05), ] -> indicator_taxa
  taxa$OTU <- rownames(taxa)
  data.frame(OTU = as.factor(row.names(indicator_taxa)), indicator_taxa) %>%
    dplyr::left_join(taxa, by="OTU") -> indicator_taxa
  rownames(indicator_taxa) <- indicator_taxa$OTU
  indicator_taxa <- arrange(indicator_taxa, desc(stat))
  return(indicator_taxa)
}

sample_data(phyloseq_merged)<- data.frame(sample_data(phyloseq_merged)) %>% mutate(transect1=ifelse(transect=="wetland-center", "wetland", transect))

View(sample_data(phyloseq_merged))

sample_data(phyloseq_merged)$transect1 <- factor(sample_data(phyloseq_merged)$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)


ind_transectCB <-
  GetIndicators(dataframe=phyloseq_merged, "transect1")


##checking

View(ind_transectCB)

ind_transectCB_subset<- ind_transectCB %>% subset(s.upland==0 & s.transition==0 & s.wetland==1)
View(ind_transectCB_subset) ##2115 Taxa uniquely associated with wetland

write.csv(ind_transectCB_subset, "ind_CB_wetland_unique_archaea.csv")



library(tidyverse)

write.csv(ind_transectCB, "indicator_transect_CB_archaea.csv")



################################
##indicator species by transect WLE 
################################

##load indicator species analysis package

install.packages("indicspecies")
library(indicspecies)

###########

##phyloseq object
library(phyloseq)
library(tidyverse)
otu = read.csv("OTU_WLE_rarefied2k_archaea.csv", sep=",", row.names=1) ##USING RAREFIED DATA here unlike unrarefied data used for site specific indicators  

otu_final = otu[rowSums(otu[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) ##495 otus


tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)

metadata = read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1) 


OTU = otu_table(otu_final, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

#bacteria
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 19262 taxa and 82 samples ]
#sample_data() Sample Data:       [ 82 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 19262 taxa by 8 taxonomic ranks ]

#archaea
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 495 taxa and 52 samples ]
#sample_data() Sample Data:       [ 52 samples by 15 sample variables ]
#tax_table()   Taxonomy Table:    [ 495 taxa by 8 taxonomic ranks ]



##indicator species analysis

library(indicspecies)

GetIndicators <-function(dataframe, var){
  require(phyloseq); require(indicspecies); require(dplyr)
  #dataframe <- tax_glom(dataframe, taxrank = "Family") ##can be adjusted to indicate at family or class level, else default is OTU level
  otu <- as.data.frame(otu_table(dataframe, taxa_are_rows = TRUE))
  metadata = as(sample_data(dataframe), "data.frame")
  taxa <- as.data.frame(as.matrix(tax_table(dataframe)))
  multipatt <- multipatt(t(otu), metadata[,var], func = "r.g",
                         control=how(nperm=999), duleg=FALSE)
  multipatt -> multipatt_fdr
  multipatt_fdr$sign$p.value <- p.adjust(multipatt_fdr$sign$p.value, "fdr")
  #print(multipatt_fdr)
  multipatt_fdr$sign[which(
    multipatt_fdr$sign$p.value <= 0.05), ] -> indicator_taxa
  taxa$OTU <- rownames(taxa)
  data.frame(OTU = as.factor(row.names(indicator_taxa)), indicator_taxa) %>%
    dplyr::left_join(taxa, by="OTU") -> indicator_taxa
  rownames(indicator_taxa) <- indicator_taxa$OTU
  indicator_taxa <- arrange(indicator_taxa, desc(stat))
  return(indicator_taxa)
}

sample_data(phyloseq_merged)<- data.frame(sample_data(phyloseq_merged)) %>% mutate(transect1=ifelse(transect=="wetland-center"|transect=="wetland-transition-edge", "wetland", transect))

View(sample_data(phyloseq_merged))

sample_data(phyloseq_merged)$transect1 <- factor(sample_data(phyloseq_merged)$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)


ind_transectWLE <-
  GetIndicators(dataframe=phyloseq_merged, "transect1")


##checking

View(ind_transectWLE)

ind_transectWLE_subset<- ind_transectWLE %>% subset(s.upland==0 & s.transition==0 & s.wetland==1)
View(ind_transectWLE_subset) ##28 Taxa uniquely associated with wetland

write.csv(ind_transectWLE_subset, "ind_WLE_wetland_unique_archaea.csv")



library(tidyverse)

write.csv(ind_transectWLE, "indicator_transect_WLE_archaea.csv")






