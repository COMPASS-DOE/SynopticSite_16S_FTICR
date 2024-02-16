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
    Domain == "d__Bacteria" &
      Family   != "f__Chloroplast" &
      Family  != "f__Mitochondria"
  )
phyloseq_merged_clean

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 46342 taxa and 215 samples ]
sample_data() Sample Data:       [ 215 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 46342 taxa by 8 taxonomic ranks ]

head(sample_data(phyloseq_merged_clean))

##inspect library sizes
df <- as.data.frame(sample_data(phyloseq_merged_clean)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(phyloseq_merged_clean)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
plot<- ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point() ##Supplemental Figure S3

ggsave(filename = "librarySize.pdf", plot = plot,
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
##FALSE  TRUE 
# 46192   150

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(phyloseq_merged_clean, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control Sample", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "True Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev0.5$contaminant)
plot<- ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") ##Supplemental Figure S3

ggsave(filename = "prevalence.pdf", plot = plot,
       width = 16,
       height = 13, units = c("cm"),
       dpi = 300)

write.csv(df.pa, "contaminant-table-cb-wle.csv")

write.csv(contamdf.prev0.5, "contaminant-prev-0.5-cb-wle.csv")

##removing contaminants from phyloseq object
df.pa <- read.csv("contaminant-prev-0.5-cb-wle.csv")
View(df.pa)
subset.df <- subset(df.pa, contaminant== "FALSE")
View(subset.df)
keep.taxa <- as.vector(subset.df$X)
keep.taxa

phyloseq_merged_clean_decontam <- prune_taxa(keep.taxa, phyloseq_merged_clean)
phyloseq_merged_clean_decontam


phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 46192 taxa and 215 samples ]
sample_data() Sample Data:       [ 215 samples by 16 sample variables ]
tax_table()   Taxonomy Table:    [ 46192 taxa by 8 taxonomic ranks ]

##export otu table out of phyloseq object 

OTU1 = as(otu_table(phyloseq_merged_clean_decontam), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU1)

write.csv(OTUdf, "OTU_clean.csv")

#remove neg controls from otu table-- might add more rows which are rowsums=0, read back in to rarefy reads to 15K

otu = read.csv("OTU_clean.csv", sep=",", row.names=1)#with no contaminants, mito, chloroplast, neg control 
tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)
metadata = read.csv("sample-metadata-cb-wle-decontam.csv", sep=",", row.names=1)

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 46192 taxa and 205 samples ]
sample_data() Sample Data:       [ 205 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 46192 taxa by 8 taxonomic ranks ]


#############################################
#####rarefaction curves-- all samples#######
#############################################


otu1 = read.csv("OTU_clean.csv", sep=",", row.names=1)##CHANGE OTU HEADER TO OTUID IN CSV FILE 
tax1 = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax1 = as.matrix(tax1)
metadata1 = read.csv("sample-metadata-cb-wle-decontam.csv", sep=",", row.names=1)
OTU1 = otu_table(otu1, taxa_are_rows = TRUE)
TAX1 = tax_table(tax1)
meta1 = sample_data(metadata1)

phyloseq_merged1 = phyloseq(OTU1, TAX1, meta1)
phyloseq_merged1

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 46192 taxa and 205 samples ]
sample_data() Sample Data:       [ 205 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 46192 taxa by 8 taxonomic ranks ]

#col<- cyan4
library(vegan)

mat<- as(t(otu_table(phyloseq_merged1)), "matrix")
system.time(rarecurve(mat, step=1000, col="cyan4", label=FALSE)) 
#rarecurve<-rarecurve(t(otu_table(phyloseq_merged1)), step=20, label = FALSE, col="cyan4") ## Supplemental Figure S3

###############################
sample_df <- data.frame(sample_data(phyloseq_merged))
View(sample_df)
write.csv(sample_df, "samples_pre_rarefied.csv")

library(vegan)
#rarecurve(t(otu_table(phyloseq_merged)), step=50)
phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged, sample.size = 20000, rngseed = TRUE, trimOTUs=FALSE)

#6 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
  
#s20230217078s20230901075s20230901078s20230901081s20230901082

sample_df <- data.frame(sample_data(phyloseq_rarefy))

View(sample_df)
write.csv(sample_df, "samples_rarefied.csv")

##use anti-join to select non-present samples in rarefied object
#This join is like eg: for df1-df2, it selects all rows from df1 that are not present in df2.
df1<- read.csv("samples_pre_rarefied.csv")
df2 <- read.csv("samples_rarefied.csv")
df= df1 %>% anti_join(df2,by="X")

View(df)

write.csv(df, "samples-lost-rarefaction.csv")

OTU2 = as(otu_table(phyloseq_rarefy), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, "OTU_clean_noneg_rarefied20k.csv")

########################
######## barplots and ordination plots
#########################



otu_final<- read.csv("OTU_clean_noneg_rarefied20k.csv", sep=",",row.names=1)

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

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 40216 taxa and 199 samples ]
sample_data() Sample Data:       [ 199 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]



##barplot class level
phyloseq_class <- phyloseq_merged %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Class) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.05] <- "< 5% abund."

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
phyloseq_class_new<- phyloseq_class %>% mutate(transect1= ifelse(transect=="wetland-center", "wetland", transect))

phyloseq_class_new$transect1 <- factor(phyloseq_class_new$transect1,levels = c("upland", "transition", "wetland", "wetland-transition-edge"), ordered=TRUE)
View(phyloseq_class_new)


##plot by region/reansect/horizon/site
plot_class<- ggplot(phyloseq_class_new, aes(x= tree_number, y = Abundance, fill =Class)) + 
  facet_nested(~region+site+transect1+horizon, scales="free", independent="x", labeller = label_wrap_gen(width=10)) + theme(strip.text.x = element_text(size = 20, face="bold"))+ 
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

ggsave(filename = "barplot_class_cb_wle_new.jpg", plot = plot_class,
       width = 80,
       height = 40, units = c("cm"),
       dpi = 300)

####ordination plots

##plot by region

metadata_new<- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)
sample_data_new<- metadata_new %>% mutate(transect1= ifelse(transect=="wetland-center", "wetland", transect))

sample_data_new$transect1 <- factor(sample_data_new$transect1,levels = c("upland", "transition", "wetland", "wetland-transition-edge"), ordered=TRUE)

sample_data_phyloseq<- sample_data(sample_data_new)
phyloseq_merged_new<- merge_phyloseq(otu, tax, sample_data_phyloseq)
phyloseq_merged_new

phyloseq_pcoa <- ordinate(
  physeq = phyloseq_merged_new, ##considering using phyloseq_rarefy directly
  method = "PCoA", 
  distance = "bray"
)


pcoa_cb_wle<-plot_ordination(
  physeq = phyloseq_merged_new,
  ordination = phyloseq_pcoa,
  color = "region",
  shape = "transect1",
  title = "PCoA WLE CB") + 
  scale_color_manual(values = c("#CC79A7", "#E69F00", "#0072B2", "black"))+
  scale_shape_manual(values=c(15,16,17,18))+
  geom_point(aes(color = region), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))

pcoa_cb_wle

ggsave(filename="pcoa-cb-wle.TIFF", plot=pcoa_cb_wle, width=8, height=6, units="in", dpi=300)


#pcoa separated by region

##cb
View(sample_data_new)
sample_data_new_cb<- sample_data_new %>% filter(region=="Chesapeake Bay")

sample_data_phyloseq_cb<- sample_data(sample_data_new_cb)
phyloseq_merged_new<- merge_phyloseq(otu, tax, sample_data_phyloseq_cb)
phyloseq_merged_new


phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 40216 taxa and 117 samples ]
sample_data() Sample Data:       [ 117 samples by 16 sample variables ]
tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]


##removing taxa with rows==0, testing the ordination with and without missing taxa

OTU2 = as(otu_table(phyloseq_merged_new), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, "OTU_CB_rarefied20k.csv")

##read back in new OTU csv for CB, remove rowsums and merge into phyloseq

otu_final<- read.csv("OTU_CB_rarefied20k.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #15,252 otus with no counts, this makes total otus less for cb than before because we have clustered otus at 99%
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)##24964 otus

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

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 24964 taxa and 117 samples ]
sample_data() Sample Data:       [ 117 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 24964 taxa by 8 taxonomic ranks ]
> 

phyloseq_pcoa <- ordinate(
  physeq = phyloseq_merged_new, ##considering using phyloseq_rarefy directly
  method = "PCoA", 
  distance = "bray"
)


pcoa_cb<-plot_ordination(
  physeq = phyloseq_merged_new,
  ordination = phyloseq_pcoa,
  color = "site",
  shape = "transect1",
  title = "PCoA CB") + 
  scale_color_manual(values = c("#CC79A7", "#E69F00", "#0072B2", "black"))+
  scale_shape_manual(values=c(15,16,17,18))+
  geom_point(aes(color = site), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))

pcoa_cb

##the plots dont change after removing rows of otus which are zero
ggsave(filename="pcoa-cb-nozeroOTUs.TIFF", plot=pcoa_cb, width=8, height=6, units="in", dpi=300)
ggsave(filename="pcoa-cb.TIFF", plot=pcoa_cb, width=8, height=6, units="in", dpi=300)

##wle
View(sample_data_new)
sample_data_new_wle<- sample_data_new %>% filter(region=="Lake Erie ")##space added in csv file by mistake
View(sample_data_new_wle)
sample_data_phyloseq_wle<- sample_data(sample_data_new_wle)
phyloseq_merged_new<- merge_phyloseq(otu, tax, sample_data_phyloseq_wle)
phyloseq_merged_new

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 40216 taxa and 82 samples ]
sample_data() Sample Data:       [ 82 samples by 16 sample variables ]
tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]

##removing taxa with rows==0, testing the ordination with and without missing taxa

OTU2 = as(otu_table(phyloseq_merged_new), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, "OTU_WLE_rarefied20k.csv")

##read back in new OTU csv for CB, remove rowsums and merge into phyloseq

otu_final<- read.csv("OTU_WLE_rarefied20k.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #20954 otus with no counts
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)##19262 otus, which is much higher than the 3000 ASVs we had before

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)

metadata_new<- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)
sample_data_new<- metadata_new %>% mutate(transect1= ifelse(transect=="wetland-center", "wetland", transect))

sample_data_new$transect1 <- factor(sample_data_new$transect1,levels = c("upland", "transition", "wetland", "wetland-transition-edge"), ordered=TRUE)

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

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19262 taxa and 82 samples ]
sample_data() Sample Data:       [ 82 samples by 16 sample variables ]
tax_table()   Taxonomy Table:    [ 19262 taxa by 8 taxonomic ranks ]
> 



phyloseq_pcoa <- ordinate(
  physeq = phyloseq_merged_new, ##considering using phyloseq_rarefy directly
  method = "PCoA", 
  distance = "bray"
)


pcoa_wle<-plot_ordination(
  physeq = phyloseq_merged_new,
  ordination = phyloseq_pcoa,
  color = "site",
  shape = "transect1",
  title = "PCoA WLE") + 
  scale_color_manual(values = c("#CC79A7", "#E69F00", "#0072B2", "black"))+
  scale_shape_manual(values=c(15,16,17,18))+
  geom_point(aes(color = site), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))

pcoa_wle

##exact same plot both times
ggsave(filename="pcoa-wle-nozeroOTUs.TIFF", plot=pcoa_wle, width=8, height=6, units="in", dpi=300)
ggsave(filename="pcoa-wle.TIFF", plot=pcoa_wle, width=8, height=6, units="in", dpi=300)


####################
####PERMANOVA stats
##################

##by region
otu_final<- read.csv("OTU_clean_noneg_rarefied20k.csv", sep=",",row.names=1)

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

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 40216 taxa and 199 samples ]
sample_data() Sample Data:       [ 199 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]



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

$aov.tab
Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
region             1    11.994 11.9940  60.867 0.15072  0.001 ***
  site               4     8.390  2.0975  10.644 0.10543  0.001 ***
  transect1          2     7.557  3.7785  19.175 0.09496  0.001 ***
  region:transect1   2     5.785  2.8924  14.678 0.07269  0.001 ***
  site:transect1     8    10.187  1.2734   6.462 0.12801  0.001 ***
  Residuals        181    35.666  0.1971         0.44819           
Total            198    79.579                 1.00000           
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

$call
adonis(formula = phyloseq_bray ~ region * site * transect1, data = sample_df, 
       by = "margin")


##Chesapeake Bay
otu_final<- read.csv("OTU_CB_rarefied20k.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #15,252 otus with no counts, this makes total otus less for cb than before because we have clustered otus at 99%
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)##24964 otus

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

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 24964 taxa and 117 samples ]
sample_data() Sample Data:       [ 117 samples by 16 sample variables ]
tax_table()   Taxonomy Table:    [ 24964 taxa by 8 taxonomic ranks ]


set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_new, method = "bray") ##using merged object is fine here, gives same results as rarefy object
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged_new)) 
View(sample_df)
set.seed(1)
adonis2(phyloseq_bray ~ site*transect1, by="terms", data = sample_df)

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ site * transect1, data = sample_df, by = "terms")
Df SumOfSqs      R2       F Pr(>F)    
site             2    6.210 0.13491 13.5090  0.001 ***
  transect1        2    7.698 0.16722 16.7453  0.001 ***
  site:transect1   4    7.302 0.15861  7.9413  0.001 ***
  Residual       108   24.825 0.53926                   
Total          116   46.035 1.00000                   
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



##Lake Erie

otu_final<- read.csv("OTU_WLE_rarefied20k.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #20954 otus with no counts
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)##19262 otus, which is much higher than the 3000 ASVs we had before

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

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19262 taxa and 82 samples ]
sample_data() Sample Data:       [ 82 samples by 16 sample variables ]
tax_table()   Taxonomy Table:    [ 19262 taxa by 8 taxonomic ranks ]

set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged_new, method = "bray") 
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged_new)) 
View(sample_df)
set.seed(1)
library(vegan)
adonis2(phyloseq_bray ~ site*transect1, by="terms", data = sample_df)

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ site * transect1, data = sample_df, by = "terms")
Df SumOfSqs      R2       F Pr(>F)    
site            2   2.1796 0.10114  7.3382  0.001 ***
  transect1       2   5.6437 0.26189 19.0007  0.001 ***
  site:transect1  4   2.8853 0.13389  4.8569  0.001 ***
  Residual       73  10.8414 0.50308                   
Total          81  21.5499 1.00000                   
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1





##################################################################
###new analysis for CB and WLE samples from compass poster########
### code adapted from May 25, 2023#################################
###################################################################


##Venn diagram  noting shared ASVs across transects

otu_final<- read.csv("OTU_clean_noneg_rarefied20k.csv", row.names=1) ##change name 1k--- depending on rarefaction curve
View(otu_final)
any(rowSums(otu_final[])<1)

##check
otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #5590 otus with no counts
View(otu1)
max(rowSums(otu1))

##subset
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) #40216 otus 

write.csv(otu_final, "otu_rarefiedCBWLE_nozerorow.csv")


library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
metadata <- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)


##subset to cb and wle as needed for the venn diagram

metadata<- metadata %>% filter(region=="Chesapeake Bay")

metadata<- metadata %>% filter(region=="Lake Erie ")

tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

##WLE samples
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 40216 taxa and 82 samples ]
sample_data() Sample Data:       [ 82 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]

##all samples
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 40216 taxa and 199 samples ]
sample_data() Sample Data:       [ 199 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]


##CB samples
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 40216 taxa and 117 samples ]
sample_data() Sample Data:       [ 117 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]

phyloseq <- phyloseq_merged %>%
  #tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  #transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt()                                       # Melt to long format
#filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
#arrange(Class) # Sort data frame alphabetically by phylum/family etc

#phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."

View(phyloseq)


#Chesapeake bay
phyloseq_gwi<- phyloseq %>% filter(site=="GWI")
phyloseq_gcrew<- phyloseq %>% filter(site=="GCREW")
phyloseq_msm<- phyloseq %>% filter(site=="MSM")

##GWI

View(phyloseq_gwi)
phyloseq_gwi_wc<- phyloseq_gwi %>% filter(transect=="wetland-center")
View(phyloseq_gwi_wc)

phyloseq_gwi_wc_nozero<- phyloseq_gwi_wc %>% filter(Abundance!=0)

View(phyloseq_gwi_wc_nozero)

phyloseq_gwi_upland<- phyloseq_gwi %>% filter(transect=="upland")
View(phyloseq_gwi_upland)

phyloseq_gwi_upland_nozero<- phyloseq_gwi_upland %>% filter(Abundance!=0)
View(phyloseq_gwi_upland_nozero)

phyloseq_gwi_transition<- phyloseq_gwi %>% filter(transect=="transition")
View(phyloseq_gwi_transition)

phyloseq_gwi_transition_nozero<- phyloseq_gwi_transition %>% filter(Abundance!=0)
View(phyloseq_gwi_transition_nozero)

#Define sets for diagram
SET1 <- unique(phyloseq_gwi_wc_nozero$OTU)
SET2 <- unique(phyloseq_gwi_upland_nozero$OTU)
SET3 <- unique(phyloseq_gwi_transition_nozero$OTU)



#Draw the diagram 


library(VennDiagram)
venn.diagram(list(Wetland=SET1, Upland=SET2, Transition=SET3),main.cex=20, 
             sub.cex=10,
             fill = c("red", "green", "blue"),
             alpha = c(0.5, 0.5, 0.5),
             filename='venn_gwi_cb.tiff', height = 6000 , 
             width = 6000 , resolution = 300, units="px", 
             imagetype="tiff", compression = "lzw", output=TRUE)

##GCREW

View(phyloseq_gcrew)
phyloseq_gcrew_wc<- phyloseq_gcrew %>% filter(transect=="wetland")
View(phyloseq_gcrew_wc)

phyloseq_gcrew_wc_nozero<- phyloseq_gcrew_wc %>% filter(Abundance!=0)

View(phyloseq_gcrew_wc_nozero)

phyloseq_gcrew_upland<- phyloseq_gcrew %>% filter(transect=="upland")
View(phyloseq_gcrew_upland)

phyloseq_gcrew_upland_nozero<- phyloseq_gcrew_upland %>% filter(Abundance!=0)
View(phyloseq_gcrew_upland_nozero)

phyloseq_gcrew_transition<- phyloseq_gcrew %>% filter(transect=="transition")
View(phyloseq_gcrew_transition)

phyloseq_gcrew_transition_nozero<- phyloseq_gcrew_transition %>% filter(Abundance!=0)
View(phyloseq_gcrew_transition_nozero)

#Define sets for diagram
SET1 <- unique(phyloseq_gcrew_wc_nozero$OTU)
SET2 <- unique(phyloseq_gcrew_upland_nozero$OTU)
SET3 <- unique(phyloseq_gcrew_transition_nozero$OTU)



#Draw the diagram 

library(VennDiagram)
venn.diagram(list(Wetland=SET1, Upland=SET2, Transition=SET3),main.cex=20, 
             sub.cex=10,
             fill = c("red", "green", "blue"),
             alpha = c(0.5, 0.5, 0.5),
             filename='venn_gcrew_cb.tiff', height = 6000 , 
             width = 6000 , resolution = 300, units="px", 
             imagetype="tiff", compression = "lzw", output=TRUE)




##MSM

View(phyloseq_msm)
phyloseq_msm_wc<- phyloseq_msm %>% filter(transect=="wetland-center")
View(phyloseq_msm_wc)

phyloseq_msm_wc_nozero<- phyloseq_msm_wc %>% filter(Abundance!=0)

View(phyloseq_msm_wc_nozero)

phyloseq_msm_upland<- phyloseq_msm %>% filter(transect=="upland")
View(phyloseq_msm_upland)

phyloseq_msm_upland_nozero<- phyloseq_msm_upland %>% filter(Abundance!=0)
View(phyloseq_msm_upland_nozero)

phyloseq_msm_transition<- phyloseq_msm %>% filter(transect=="transition")
View(phyloseq_msm_transition)

phyloseq_msm_transition_nozero<- phyloseq_msm_transition %>% filter(Abundance!=0)
View(phyloseq_msm_transition_nozero)

#Define sets for diagram
SET1 <- unique(phyloseq_msm_wc_nozero$OTU)
SET2 <- unique(phyloseq_msm_upland_nozero$OTU)
SET3 <- unique(phyloseq_msm_transition_nozero$OTU)



#Draw the diagram 

library(VennDiagram)
venn.diagram(list(Wetland=SET1, Upland=SET2, Transition=SET3),main.cex=20, 
             sub.cex=10,
             fill = c("red", "green", "blue"),
             alpha = c(0.5, 0.5, 0.5),
             filename='venn_msm_cb.tiff', height = 6000 , 
             width = 6000 , resolution = 300, units="px", 
             imagetype="tiff", compression = "lzw", output=TRUE)

#new venn 

venn.diagram(x=list(Wetland=SET1, Upland=SET2, Transition=SET3),
             category.names = c("Wetland" , "Upland" , "Transition"),
             filename = 'venn_msm_cb.png',
             output = TRUE ,
             imagetype="png" ,
             height = 1000 , 
             width = 1000 , 
             resolution = 1000,
             compression = "lzw",
             lwd = 1,
             fill=c("#f2ac79", "#5ba640", "#eed580"),
             col = c(alpha("#f2ac79",0.3), alpha("#5ba640",0.3), alpha("#eed580",0.3)),
             cex = 0.5,
             fontfamily = "sans",
             cat.cex = 0.3,
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             #cat.col = c("#f2ac79", "#5ba640", "#eed580"),
             rotation = 1
)


##Lake Erie

phyloseq_pr<- phyloseq %>% filter(site=="PR")
phyloseq_cc<- phyloseq %>% filter(site=="CC")
phyloseq_owc<- phyloseq %>% filter(site=="OWC")

##PR

View(phyloseq_pr)
phyloseq_pr_wc<- phyloseq_pr %>% filter(transect=="wetland-center")
View(phyloseq_pr_wc)

phyloseq_pr_wc_nozero<- phyloseq_pr_wc %>% filter(Abundance!=0)

View(phyloseq_pr_wc_nozero)

phyloseq_pr_upland<- phyloseq_pr %>% filter(transect=="upland")
View(phyloseq_pr_upland)

phyloseq_pr_upland_nozero<- phyloseq_pr_upland %>% filter(Abundance!=0)
View(phyloseq_pr_upland_nozero)

phyloseq_pr_transition<- phyloseq_pr %>% filter(transect=="transition")
View(phyloseq_pr_transition)

phyloseq_pr_transition_nozero<- phyloseq_pr_transition %>% filter(Abundance!=0)
View(phyloseq_pr_transition_nozero)

#Define sets for diagram
SET1 <- unique(phyloseq_pr_wc_nozero$OTU)
SET2 <- unique(phyloseq_pr_upland_nozero$OTU)
SET3 <- unique(phyloseq_pr_transition_nozero$OTU)


## CC

View(phyloseq_cc)
phyloseq_cc_wc<- phyloseq_cc %>% filter(transect=="wetland-center")
View(phyloseq_cc_wc)

phyloseq_cc_wc_nozero<- phyloseq_cc_wc %>% filter(Abundance!=0)

View(phyloseq_cc_wc_nozero)

phyloseq_cc_upland<- phyloseq_cc %>% filter(transect=="upland")
View(phyloseq_cc_upland)

phyloseq_cc_upland_nozero<- phyloseq_cc_upland %>% filter(Abundance!=0)
View(phyloseq_cc_upland_nozero)

phyloseq_cc_transition<- phyloseq_cc %>% filter(transect=="transition")
View(phyloseq_cc_transition)

phyloseq_cc_transition_nozero<- phyloseq_cc_transition %>% filter(Abundance!=0)
View(phyloseq_cc_transition_nozero)

#Define sets for diagram
SET1 <- unique(phyloseq_cc_wc_nozero$OTU)
SET2 <- unique(phyloseq_cc_upland_nozero$OTU)
SET3 <- unique(phyloseq_cc_transition_nozero$OTU)



##OWC

View(phyloseq_owc)
phyloseq_owc_wte<- phyloseq_owc %>% filter(transect=="wetland-transition-edge")
View(phyloseq_owc_wte)

phyloseq_owc_wte_nozero<- phyloseq_owc_wte %>% filter(Abundance!=0)

View(phyloseq_owc_wte_nozero)

phyloseq_owc_upland<- phyloseq_owc %>% filter(transect=="upland")
View(phyloseq_owc_upland)

phyloseq_owc_upland_nozero<- phyloseq_owc_upland %>% filter(Abundance!=0)
View(phyloseq_owc_upland_nozero)

phyloseq_owc_transition<- phyloseq_owc %>% filter(transect=="transition")
View(phyloseq_owc_transition)

phyloseq_owc_transition_nozero<- phyloseq_owc_transition %>% filter(Abundance!=0)
View(phyloseq_owc_transition_nozero)

#Define sets for diagram
SET1 <- unique(phyloseq_owc_wte_nozero$OTU)
SET2 <- unique(phyloseq_owc_upland_nozero$OTU)
SET3 <- unique(phyloseq_owc_transition_nozero$OTU)


#new venn 

venn.diagram(x=list(Wetland=SET1, Upland=SET2, Transition=SET3),
             category.names = c("Wetland" , "Upland" , "Transition"),
             filename = 'venn_owc_wle.png',
             output = TRUE ,
             imagetype="png" ,
             height = 1000 , 
             width = 1000 , 
             resolution = 1000,
             compression = "lzw",
             lwd = 1,
             fill=c("#f2ac79", "#5ba640", "#eed580"),
             col = c(alpha("#f2ac79",0.3), alpha("#5ba640",0.3), alpha("#eed580",0.3)),
             cex = 0.5,
             fontfamily = "sans",
             cat.cex = 0.3,
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             #cat.col = c("#f2ac79", "#5ba640", "#eed580"),
             rotation = 1
)


#########################################################
################Alpha Diversity #########################
#########################################################

otu_final<- read.csv("OTU_clean_noneg_rarefied20k.csv", row.names=1) ##change name 1k--- depending on rarefaction curve
View(otu_final)
any(rowSums(otu_final[])<1)

##check
otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #5590 otus with no counts
View(otu1)
max(rowSums(otu1))

##subset
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) #40216 otus 

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
metadata <- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)


##subset to cb and wle as needed for the venn diagram

metadata<- metadata %>% filter(region=="Chesapeake Bay")

metadata<- metadata %>% filter(region=="Lake Erie ")

tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 40216 taxa and 199 samples ]
sample_data() Sample Data:       [ 199 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]


##richness (observed OTUs) estimates, inverse simpson

library(ggplot2)
library(dplyr)
library(vegan)
library(reshape2)
library(scales)
library(grid)
library(readxl)


# Initialize matrices to store richness and evenness estimates
nsamp = nsamples(phyloseq_merged)
trials = 100
richness <- matrix(nrow = nsamp, ncol = trials)
row.names(richness) <- sample_names(phyloseq_merged)
evenness <- matrix(nrow = nsamp, ncol = trials)
row.names(evenness) <- sample_names(phyloseq_merged)
# It is always important to set a seed when you subsample so your result is replicable, note of caution that this will yield different results of alpha diversity depending on the version of R being used. This is because set.seed function can vary across R versions. Here I am reporting results from R v.3.4.0 
set.seed(3)
for (i in 1:100) {
  # Subsample
  #r <- rarefy_even_depth(phyloseq_merged, sample.size = 40000, verbose = FALSE, replace = TRUE) #no rarefaction a second time.
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(phyloseq_merged, measures = "Observed"))) ##changed r to phyloseq_merged, no second rarefaction
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(phyloseq_merged, measures = "InvSimpson"))) ##changed r to phyloseq_merged, no second rarefaction
  evenness[ ,i] <- even
}
# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)
# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd) 
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)
alpha <- rbind(rich_stats, even_stats)

s <- read.csv("sample-metadata-cb-wle-alphadiv.csv")
alphadiv <- merge(alpha, s, by = "SampleID") 
write.csv(alphadiv,file = "alphadiv_CB_WLE.csv")


##subset if needed
alphadiv_wle<-alphadiv[alphadiv$region=="Lake Erie ",]
dim(alphadiv_wle)
alphadiv_cb<-alphadiv[alphadiv$region=="Chesapeake Bay",]
dim(alphadiv_cb)
View(alphadiv_wle)
View(alphadiv_cb)
write.csv(alphadiv_wle, file="alphadiv_wle.csv")
write.csv(alphadiv_cb, file="alphadiv_cb.csv")


##average if needed
group_wle<-alphadiv_wle%>%
  group_by(measure, transect, site) %>%
  summarise(mean_update = mean(mean))
write.csv(group_wle, file="alphadivmean_wle.csv")
group_cb<-alphadiv_cb%>%
  group_by(measure, transect, site) %>%
  summarise(mean_update = mean(mean))
write.csv(group_cb, file="alphadivmean_cb.csv")




### Inverse Simpson CB
alphadiv_CB <- read.csv("alphadiv_cb.csv")
pd<-position_dodge(0.7)
alphadiv_CB_InvS =  alphadiv_CB %>% filter(measure=="Inverse Simpson")

View(alphadiv_CB_InvS)


#### Richness CB
alphadiv_CB <- read.csv("alphadiv_cb.csv")
pd<-position_dodge(0.7)
alphadiv_CB_rich =  alphadiv_CB %>% filter(measure=="Richness")

View(alphadiv_CB_rich)



##combine richness and evenness together in panel


View(alphadiv_CB_InvS)


View(alphadiv_CB_rich)

alphadiv_combined= rbind(alphadiv_CB_rich, alphadiv_CB_InvS)
View(alphadiv_combined)



write.csv(alphadiv_combined, "alphadiv_combined.csv")
alphadiv_combined<- read.csv("alphadiv_combined.csv")

alphadiv_new<- alphadiv %>% mutate(transect1=ifelse(transect=="WC", "wetland", transect))

alphadiv_new$transect1 <- factor(alphadiv_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)

plot_CB<- ggplot(alphadiv_new, aes(x= transect1,y=mean,color=site)) + theme_classic()+ scale_color_manual(values=c("#39568CFF", "#55C667FF", "orange")) +
  # scale_shape_manual(values=c(15,17, 19))+ #scale_size(range = c(1,10))+ # point size range
  #label = scales::dollar)+
  #expand_limits(shape = seq(2, 6, by = 1))+
  #scale_alpha_manual(values = c(2,3,4,5,6))+
  geom_boxplot() + geom_jitter()+ 
  facet_wrap(~measure)+theme(
    strip.text.x = element_text(
      size = 16
    ))+
  theme(plot.title = element_text(size= 18, hjust = 0.5))+
  labs(y = "measure")+ xlab("transect")+theme(axis.text.x = element_text(colour = "black", size = 17,angle = 90, vjust = 0.5, hjust=1))+
  theme( axis.text.y = element_text(size = 16, colour = "black"), 
         axis.title.y= element_text(size=16, color="black"),
         axis.title.x= element_text(size=16, color="black"),
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =14, colour = "black"),
         legend.text = element_text(size = 12, colour = "black")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")), color= guide_legend(order = 1)) 


plot_CB


##Figure
ggsave(filename="alphadiv_rich_even_CB_1.TIFF", plot=plot_CB, width=25, height=15, units="cm", dpi=300)

##alpha diversity ANOVA stats
library(ggplot2)
library(grid)
library(lattice)
library(multcompView)
library(tidyverse)
library(dplyr)
alphadiv_rich <- read.csv("alphadiv_cb.csv") ## combined richness and evenness estimates of all samples


pd<-position_dodge(0.7)


alphadiv_rich =  alphadiv_rich %>% filter(measure=="Inverse Simpson")
alphadiv_rich =  alphadiv_rich %>% filter(measure=="Richness")

##Conduct the following code for each of the datasets selected above
View(alphadiv_rich)

alphadiv_rich<- alphadiv_rich %>% mutate(transect1=ifelse(transect=="wetland-center", "wetland", transect))
options(contrasts = c("contr.sum", "contr.poly")) ### must set this before running the model for type III tests
model = lm(mean~ transect1+site + transect1:site +transect1 + transect1:site,
           data=alphadiv_rich)

shapiro.test(resid(model)) #W statistic > 0.9 means normality okay #check normality to flag outliers 

##richness cb

Shapiro-Wilk normality test

data:  resid(model)
W = 0.96946, p-value = 0.009037


##inverse simpson without horizon in model ##CB
Shapiro-Wilk normality test

data:  resid(model)
W = 0.89682, p-value = 1.811e-07


##richness without horizon in model ##WLE

Shapiro-Wilk normality test

data:  resid(model)
W = 0.98289, p-value = 0.3457

#inverse simpson WLE

Shapiro-Wilk normality test

data:  resid(model)
W = 0.98605, p-value = 0.5201



ee= as.matrix(resid(model))
qqnorm (ee) #look for outliers (far away from curve) and test for normality 
ee

library(car)


#options(contrasts = c("contr.sum", "contr.poly"))
Anova(model, type=c("III"))

####Results inverse simpson CB, works without horizon as horizon levels are unequal among transects and sites

Anova Table (Type III tests)

Response: mean
Sum Sq  Df  F value    Pr(>F)    
(Intercept)    1647162   1 503.5234 < 2.2e-16 ***
  transect1       226999   2  34.6959 2.304e-12 ***
  site             43586   2   6.6619  0.001870 ** 
  transect1:site   59530   4   4.5495  0.001958 ** 
  Residuals       353297 108                       
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##cb richness

Anova Table (Type III tests)

Response: mean
Sum Sq  Df   F value    Pr(>F)    
(Intercept)    93988748   1 1982.6360 < 2.2e-16 ***
  transect1       6990289   2   73.7280 < 2.2e-16 ***
  site             240632   2    2.5380   0.08373 .  
transect1:site  1768344   4    9.3255 1.638e-06 ***
  Residuals       5119843 108                        
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> 



##Richness WLE
Anova Table (Type III tests)

Response: mean
Sum Sq Df   F value    Pr(>F)    
(Intercept)    119795971  1 2026.8848 < 2.2e-16 ***
  transect1        6894182  2   58.3230 7.351e-16 ***
  site             1038481  2    8.7853 0.0003813 ***
  transect1:site   2348455  4    9.9337 1.791e-06 ***
  Residuals        4314555 73                        
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##Inverse simpson WLE

Anova Table (Type III tests)

Response: mean
Sum Sq Df  F value    Pr(>F)    
(Intercept)    3130132  1 548.3830 < 2.2e-16 ***
  transect1       377590  2  33.0759 5.942e-11 ***
  site             56154  2   4.9189  0.009907 ** 
  transect1:site  174906  4   7.6607 3.279e-05 ***
  Residuals       416679 73                       
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




##add letters on graph 
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(multcompView)
library(multcomp)


#set up model ##letters will be evaluated depending on what factor is different based on the anova above
model <- lm(mean~transect1, data = alphadiv_rich)

# get (adjusted) weight means per group
library(emmeans)
model_means <- emmeans(object = model,
                       specs = "transect1")

cld <- cld(object = model_means,
           adjust = "Tukey",
           Letters = letters,
           alpha = 0.05)

View(cld)

combined<- full_join(alphadiv_rich,cld, by=c("transect1"))
View(combined)
###

#plot<-combined%>%
#group_by("Soil")%>% unique(combined$.group)>1

library(dplyr)
combined<- as.tibble(combined)
class(combined)

data<- combined %>% mutate(group=.group)
View(data)
data$group <-str_remove(data$group, " ") ##run thrice

View(data)

df_new<- as.tibble(data)

class(df_new)
View(df_new)
as.tibble(df_new)
#data1<- df_new %>% dplyr::group_by(Site)%>%
#mutate(count= length(unique(group)), new= ifelse(count!=1, df_new$Soil, "combined"))
#df_new$new <- ifelse(length(unique(df_new$group))>1 %in% df_new$Site, df_new$Soil)

View(df_new)

df_new<- as.data.frame(df_new)

length(unique(df_new$group))

colnames(df_new)

#data1<- data1%>% unite(soil_site, Soil,Site, sep="_", remove=FALSE)

#df_new<- df_mew %>% mutate(transect1=ifelse(transect=="WC", "wetland", transect))
summarized <- df_new %>% dplyr::select(c('transect1', 'mean', 'group')) %>% group_by(transect1)%>% summarise(max_mean=max(mean), across(group))%>% distinct(.keep_all=TRUE)
View(summarized)

#summarized<- summarized %>% mutate(transect1=ifelse(transect=="WC", "wetland", transect))

summarized$transect1 <- factor(summarized$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)



#alphadiv_rich$site<-factor(alphadiv_rich$site,levels = c("CC", "PR", "OWC"), ordered=TRUE)      

#alphadiv_rich$transect<-factor(alphadiv_rich$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)         


df_new$transect1<- factor(df_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)         

p<- ggplot(df_new, aes(transect1, 
                       mean, 
                       colour = site)) +
  theme_classic()+ scale_color_manual(values=c("#39568CFF", "#55C667FF", "orange"))+
  geom_boxplot(outlier.shape = NA,
               fill = "white", 
               outlier.colour = NA, 
               position = position_dodge(width = 0.9)) + 
  geom_point(position = position_jitterdodge()) +
  scale_x_discrete(name = 'transect1', 
                   breaks = c('upland', 'transition', 'wetland'), 
                   labels = c('upland', 'transition', 'wetland'))+
  geom_text(data=summarized, aes(x=summarized$transect1, y =1890, 
                                 group=summarized$transect1, label = summarized$group, vjust=-0.3), 
            check_overlap = T, 
            position = position_dodge(width = 0.9), col = "black") + 
  ylab("Richness (Number of Observed OTUs)") +
  theme(axis.text.x = element_text(colour = "black", size = 17,angle = 90, vjust = 0.5, hjust=1))+
  theme( axis.text.y = element_text(size = 16, colour = "black"), 
         axis.title.y= element_text(size=16, color="black"),
         axis.title.x= element_text(size=16, color="black"),
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =14, colour = "black"),
         legend.text = element_text(size = 12, colour = "black")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")), color= guide_legend(order = 1)) +
  xlab("Transect") 


p

ggsave(filename = "alphadiv_letters_CB_richness.tiff", plot = p,
       width = 15,
       height = 15, units = c("cm"),
       dpi = 300)

ggsave(filename = "alphadiv_letters_CB_InvSimp.tiff", plot = p,
       width = 15,
       height = 15, units = c("cm"),
       dpi = 300)

ggsave(filename = "alphadiv_letters_WLE_Richness_ASV.tiff", plot = p,
       width = 15,
       height = 15, units = c("cm"),
       dpi = 300)

ggsave(filename = "alphadiv_letters_WLE_inversesimpson.tiff", plot = p,
       width = 15,
       height = 15, units = c("cm"),
       dpi = 300)



#######################################
###CAP analysis with soil chemistry data
#######################################


##sample names check with microbiome and soil chemistry data

data<- read.csv("allhorizon_soilchem.csv")
data_cb<- data %>% filter(region=="Chesapeake Bay")

unique(data_cb$SampleID) ##125 SAMPLES

data<- read.csv("allhorizon_soilchem.csv")
data_wle<- data %>% filter(region=="Lake Erie " & horizon=="A")

unique(data_wle$SampleID) ##95 SAMPLES


##reading in new phyloseq object with all samples in CB
otu_final<- read.csv("OTU_CB_rarefied20k.csv", row.names=1) ##change name to WLE or CB as intended
View(otu_final)
any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #20954 otus with no counts (wle), 15252 otus with no counts (cb)
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) #19262 otus for WLE, 24964 otus for CB

data_chem<- read.csv("allhorizon_soilchem.csv")
View(data_chem)
metadata_cb_wle<- read.csv("sample-metadata-cb-wle-cap.csv")
View(metadata_cb_wle)
metadata_cap<- merge(data_chem, metadata_cb_wle, by="SampleID")

View(metadata_cap)

unique(metadata_cap$SampleID) #205 samples

samples_lost<- anti_join(metadata_cb_wle, metadata_cap, by = "SampleID")
View(samples_lost) ##negative controls lost only


##subset to CB samples

cap_cb<- metadata_cap %>% filter(region.x=="Chesapeake Bay")

View(cap_cb)
##subset to WLE samples
cap_wle<- metadata_cap %>% filter(region.x=="Lake Erie ")


View(cap_wle)


##melt to wide format 

cap_cb_new<- cap_cb%>% dplyr::select(!analysis) %>% pivot_wider(names_from="name", values_from="value")  
View(cap_cb_new) ##118 samples
cap_wle_new<- cap_wle%>% dplyr::select(!analysis) %>% pivot_wider(names_from="name", values_from="value")  
View(cap_wle_new) #87 samples

write.csv(cap_cb_new, "cap_cb_new.csv") ##did not save this a second time

write.csv(cap_wle_new, "cap_wle_new.csv") 

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
metadata <- read.csv("cap_cb_new.csv", sep=",", row.names=8) ##change as needed for cb and wle

##removing columns that are >50% NAs, adjust separately for CB and WLE
##WLE
metadata_new<- metadata %>% dplyr::select(!c(region.y, Sample.description, #Nitrate_meq100g, #Sulfate_meq100g, Phosphate_meq100g, #Bromide_meq100g, Ammonia_meq100g
                                             tree_number.y, transect.y, site.y, horizon.y, project, date, DNA.Well, DNA.Plate, PCR.Plate, PCR.Well, Index.ID, Barcode, Potassium_meq100g, dic_ugg, Sodium_meq100g, Fluoride_meq100g, Chloride_meq100g, Magnesium_meq100g, Calcium_meq100g,
                                             percentOM)) 
#CB
metadata_new<- metadata %>% dplyr::select(!c(region.y, Sample.description, #Nitrate_meq100g, #Sulfate_meq100g, #Bromide_meq100g, Ammonia_meq100g, Phosphate_meq100g,
                                             tree_number.y, transect.y, site.y, horizon.y, project, date, DNA.Well, DNA.Plate, PCR.Plate, PCR.Well, Index.ID, Barcode,Potassium_meq100g, dic_ugg, Sodium_meq100g, Fluoride_meq100g, Chloride_meq100g, Magnesium_meq100g, Calcium_meq100g,
                                             percentOM)) 


View(metadata_new)

##removing rows that are NAs, need to do this because ordinate function does not work with NAs
metadata_new=na.omit(metadata_new)

#metadata_new<- na.pass(metadata)
View(metadata_new)
tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata_new) ## 76 samples for wle because omit na removed some, 67 samples remain for CB

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged


  
  ##after NA are removed, for WLE

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19262 taxa and 71 samples ]
sample_data() Sample Data:       [ 71 samples by 29 sample variables ]
tax_table()   Taxonomy Table:    [ 19262 taxa by 8 taxonomic ranks ]
 
## after NA are removed for cb

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 24964 taxa and 66 samples ]
sample_data() Sample Data:       [ 66 samples by 30 sample variables ]
tax_table()   Taxonomy Table:    [ 24964 taxa by 8 taxonomic ranks ]
  
 
##weighted unifrac distance matrix needs to be computed for CAP plot to work below
set.seed(1)
phyloseq_bray <- phyloseq::distance(phyloseq_merged, method = "bray")
phyloseq_bray

library("ape")
set.seed(1)
random_tree = rtree(ntaxa(phyloseq_merged), rooted=TRUE, tip.label=taxa_names(phyloseq_merged))
plot(random_tree)
phyloseq_tree = merge_phyloseq(phyloseq_merged, random_tree)
phyloseq_tree

##CB 
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 24964 taxa and 66 samples ]
sample_data() Sample Data:       [ 66 samples by 30 sample variables ]
tax_table()   Taxonomy Table:    [ 24964 taxa by 8 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 24964 tips and 24963 internal nodes ]


##WLE


phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19262 taxa and 71 samples ]
sample_data() Sample Data:       [ 71 samples by 29 sample variables ]
tax_table()   Taxonomy Table:    [ 19262 taxa by 8 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 19262 tips and 19261 internal nodes ]


#calculate weighted unifrac distance matrix
phyloseq_wunifrac <- phyloseq::distance(phyloseq_tree, method = "wunifrac" )


##CAP plot

# First select the variables that significantly explain variation in community composition.
## To do this I'm using the ordistep function which steps through the variables removing ones that don't add to the significance. 
## There are better, more thurough explanations of that online.
### First set new names for all variables that fit nicer on the final figure.
var_goodnames = data.frame(labels = c("gwc_perc", "nh4n_ug_g", "no3n_ug_g",
                                      "spConductance_mscm", "pH", 
                                      "TC_perc", "TN_perc", "TS_perc", "cec_meq100g",
                                      "Na_meq100g", "Nitrate_meq100g", "Sulfate_meq100g", "Phosphate_meq100g", 
                                      "Ammonia_meq100g"), ##chloride only present for may samples in CB
                           goodnames = c("Grav.Moisture", "NH4","NO3", "sp.Cond",
                                         "pH", "TC", "TN", "TS", "CEC", "Na",
                                         "NO3", "SO4","PO4", "NH3"))

#WLE only
var_goodnames = data.frame(labels = c("gwc_perc", "nh4n_ug_g", "no3n_ug_g",
                                      "spConductance_mscm", "pH", 
                                      "TC_perc", "TN_perc", "TS_perc", "cec_meq100g",
                                      "Na_meq100g", "Nitrate_meq100g", "Sulfate_meq100g", "Phosphate_meq100g", 
                                      "Ammonia_meq100g"), ##chloride only present for may samples in CB
                           goodnames = c("Grav.Moisture", "NH4","NO3", "sp.Cond",
                                         "pH", "TC", "TN", "TS", "CEC", "Na", 
                                         "NO3", "SO4","PO4", "NH3"))



### Full model with all variables
set.seed(4242)
cap_ord.full <- ordinate(physeq = phyloseq_merged, method = "CAP", distance = phyloseq_wunifrac, 
                         formula = ~ gwc_perc + nh4n_ug_g + no3n_ug_g +
                           spConductance_mscm + pH + TC_perc + TN_perc + TS_perc + cec_meq100g  + 
                           Na_meq100g + Nitrate_meq100g + Sulfate_meq100g + Phosphate_meq100g +
                           Ammonia_meq100g)
### Null model with no variables
set.seed(4242)
cap_ord.null <- ordinate(physeq = phyloseq_merged, method = "CAP", distance = phyloseq_wunifrac, 
                         formula = ~ 1)

class(cap_ord.null)

class(cap_ord.full)

### Model selection to get just significant variables
set.seed(4242)
library(vegan)
ordistep.res = ordistep(cap_ord.null, scope = formula(cap_ord.full), perm.max = 1000, trace=F)
goodform = ordistep.res$call$formula

# Get main CAP ordination
## This uses phyloseq package. The formula here is generated by that ordistep function.
set.seed(4242)
cap_ord <- ordinate(physeq = phyloseq_merged, method = "CAP", distance = phyloseq_wunifrac, formula = goodform)

# CAP plot CB
## For this I use my own code rather than phyloseq plotting for more control. To do this I have to pull out the x (CAP1) and y (CAP2) axes from the ordination and add in my metadata
key_match<- read.csv("key_match_cap.csv")
View(key_match)
cap.ord.df1 <- data.frame(vegan::scores(cap_ord, display="sites")) %>%
  tibble::rownames_to_column(var="sample.id") %>%
  dplyr::select(sample.id,CAP1, CAP2) %>% 
  left_join(key_match, by="sample.id")
cap.ord.df = cap.ord.df1 %>%
  left_join(sample_data(phyloseq_merged), by = "SampleID")

# CAP plot WLE
## For this I use my own code rather than phyloseq plotting for more control. To do this I have to pull out the x (CAP1) and y (CAP2) axes from the ordination and add in my metadata

key_match<- read.csv("key_match_cap.csv")
View(key_match)
cap.ord.df1 <- data.frame(vegan::scores(cap_ord, display="sites")) %>%
tibble::rownames_to_column(var="sample.id") %>%
  dplyr::select(sample.id,CAP1, CAP2) %>% 
  left_join(key_match, by="sample.id")
cap.ord.df = cap.ord.df1 %>%
  left_join(sample_data(phyloseq_merged), by = "SampleID")

View(cap.ord.df)
## Calculate eignevalues and fraction of variation explained by each CAP axis. I use this for the axis labels in the plot
eigvec = vegan::eigenvals(cap_ord)
fracvar = round(eigvec/sum(eigvec)*100, 2)

## Plot initial figure of points
library(ggplot2)
cap.ord.df$transect.x <- factor(cap.ord.df$transect.x,levels = c("upland", "transition", "wetland"), ordered=TRUE)



cap.ord.df<- cap.ord.df %>% filter(!SampleID=="COMPASS.Dec2021.016")##filter out the outlier for wle COMPASS.Dec2021.016
cap_plot<- ggplot(data=cap.ord.df, aes(x=CAP1, y=CAP2)) +
  geom_point(aes(shape=site.x, fill=transect.x), size=3)+
  scale_shape_manual(values=c("MSM"=24, "GWI"=21))+##change site names for CB
  scale_fill_manual(values=c("#5ba640", "#eed580", "#f2ac79", "red"))+
  labs(x=paste("CAP1 (", fracvar[1], "%)", sep=""),y=paste("CAP2 (", fracvar[2], "%)", sep=""))+
  labs(shape="site", fill="transect")


#cap_plot$labels$site.x <- "site"
#cap_plot$labels$transect.x <- "transect"

cap_plot


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

## Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat) %>%
  mutate(labels = gsub("\\.", ":", labels))
colnames(arrowdf) = c("labels", "xend", "yend")
arrowdf = arrowdf %>%
  left_join(var_goodnames, by = "labels") %>%
  rename(old_labels = labels) %>%
  rename(labels = goodnames)

## Define the arrow aesthetic mapping
arrow_map <- aes(xend = xend, yend = yend, x = 0, y = 0, 
                 color = NULL)

label_map <- aes(x = xend + 0.02*xend/abs(xend), y = yend, 
                 color = NULL, label = labels)

arrowhead = arrow(length = unit(0.02, "npc"), type = "closed")

# Make a new graphic including labels and axes for variables.
cap.plot = cap_plot + 
  geom_segment(mapping = arrow_map, linewidth = 1.2, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_segment(mapping = arrow_map, linewidth = 0.5, data = arrowdf, color = "orange", arrow = arrowhead) + 
  geom_label(mapping = label_map, data = filter(arrowdf, xend < 0), show.legend = FALSE, size=6*5/14, hjust=1, fill="orange", color="black") +
  geom_label(mapping = label_map, data = filter(arrowdf, xend > 0), show.legend = FALSE, size=6*5/14, hjust=0, fill="orange", color="black") +
  labs(title = paste("Variables explain ", round(100*RsquareAdj(cap_ord)$r.squared, 3), "% of the variation", sep="")) +
  theme_classic() +
  theme(legend.position="bottom",
        legend.direction = "vertical") +
  guides(shape = guide_legend(order = 1, nrow=1),
         fill = guide_legend(order = 2, override.aes=list(shape=22), nrow=1))

cap.plot
ggsave(filename = "cap_cb.jpg", plot = cap.plot,
       width = 19,
       height = 19, units = c("cm"),
       dpi = 700)

##DESeq and Lefse analysis for differentially abundant taxa

##LEFSE in R, similar to DeSeq


otu_final<- read.csv("OTU_clean_noneg_rarefied20k.csv", row.names=1) ##change depending on region, all or cb or wle
View(otu_final)

otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) #24964 otus CB, 19262 for WLE, 40216 for cb and wle


library(phyloseq)


tax<- read.csv("taxonomy-dn-99-lefse.csv", sep=",", row.names=1)
metadata <- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)
metadata_new<- metadata %>% mutate(transect1= ifelse(transect=="wetland-center", "wetland", transect))
#CB
metadata_new$transect1 <- factor(metadata_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)
#WLE
metadata_new$transect1 <- factor(metadata_new$transect1,levels = c("upland", "transition", "wetland", "wetland-transition-edge"), ordered=TRUE)


View(metadata_new)

tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata_new)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

##CB

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 24964 taxa and 117 samples ]
sample_data() Sample Data:       [ 117 samples by 16 sample variables ]
tax_table()   Taxonomy Table:    [ 24964 taxa by 7 taxonomic ranks ]


##WLE
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19262 taxa and 82 samples ]
sample_data() Sample Data:       [ 82 samples by 16 sample variables ]
tax_table()   Taxonomy Table:    [ 19262 taxa by 7 taxonomic ranks ]

##both cb and wle

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 40216 taxa and 199 samples ]
sample_data() Sample Data:       [ 199 samples by 16 sample variables ]
tax_table()   Taxonomy Table:    [ 40216 taxa by 7 taxonomic ranks ]


BiocManager::install("microbiomeMarker")

library(microbiomeMarker)

##CB
lefse_transect_cb<-run_lefse(phyloseq_merged, group="transect1", subgroup="site", taxa_rank = "Genus", transform="identity", norm="none") ##no marker was identified when grouping by site and subgroup is transect
##markers were identified with group as transect and subgroup as site

#WLE
lefse_transect_wle<-run_lefse(phyloseq_merged, group="transect1", subgroup="site", taxa_rank = "Genus", transform="identity", norm="none")

##both regions
lefse_transect_cbwle<-run_lefse(phyloseq_merged, group="transect1", subgroup="site", taxa_rank = "Genus", transform="identity", norm="none")

marker_table(lefse_transect_cb) ##11 markers at class level for CB, 23 markers at genus level for CB

marker_table(lefse_transect_wle) ##23 markers at genus level for wle

marker_table(lefse_transect_cbwle) #10 markers at genus level for cb and wle

plot_lefse<- plot_abundance(lefse_transect_cbwle, group="transect1") ##change cb to wle as needed
plot_lefse

ggsave(filename = "lefse_plot_abund_transect_cbwle_GENUS.tiff", plot = plot_lefse,
       width = 15,
       height = 19, units = c("cm"),
       dpi = 300)



#Bar plot or dot plot for effect size
# bar plot
plot_effect_size<-plot_ef_bar(lefse_transect_cbwle)

ggsave(filename = "lefse_effect_size_transect_cbwle_GENUS.tiff", plot = plot_effect_size,
       width = 15,
       height = 19, units = c("cm"),
       dpi = 300)


##DESeq to compare with lefse results , location differences only

library(DESeq2)
packageVersion("DESeq2")


phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 40216 taxa and 199 samples ]
sample_data() Sample Data:       [ 199 samples by 16 sample variables ]
tax_table()   Taxonomy Table:    [ 40216 taxa by 7 taxonomic ranks ]


##convert phyloseq format to DESEq dataset with dispersions estimated using the experimental design formula
diagdds = phyloseq_to_deseq2(phyloseq_merged, ~ region)
diagdds$region
diagdds

diagdds$region <-factor(diagdds$region, levels = c("Chesapeake Bay", "Lake Erie "))
diagdds$region #make sure that Control is the first level in the treatment factor, so that the
#default log2 fold changes are calculated as treatment over control and not the other way around.

cts <- counts(diagdds)

##if all rows have at least one zero then the DeSeq function fails because it cant compute geometric means 
##f only a few rows have this issue then the function still works by using the other rows, otherwise it does not.
##since in this case all rows have at least one zero, mostly due to low read counts in WLE I am using the below code 
##suggested by Michael Love (developer DeSeq2)
##this code basically removes all rows where all entries are zero, 
##else it subsets to those cells in a row which are NOT zero, and does geometric mean on log of those values. 
##this is because log of Zero is undefined

##DeSeq error

##estimating size factors
##Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : 
##every gene contains at least one zero, cannot compute log geometric means


geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
diagdds <- estimateSizeFactors(diagdds, geoMeans=geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

##The default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function.

#The following results function call creates a table of the results of the tests. Very fast. 
#The hard work was already stored with the rest of the DESeq2-related data in our latest version 
#of the diagdds object (see above). I then order by the adjusted p-value, removing the entries 
#with an NA value. The rest of this example is just formatting the results table with taxonomic information for nice(ish) display.


res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseq_merged)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

##dont do this if wanting both positive and negatives
sigtab_positives<- sigtab %>% filter(!log2FoldChange<0)
View(sigtab_positives) 


#Let's look at the OTUs that were positively enriched in CB AND WLE. The following makes a nice ggplot2 summary of the results.

library("ggplot2")
install.packages("Polychrome")
library(Polychrome)
P52<- createPalette(52,  c("#010101", "#ff0000"))
names(P52) <- NULL
library(viridis)


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "P52", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

View(sigtab) 


write.csv(sigtab, "sigtab_CB_WLE_enriched.depleted.inCB.csv")

## Figure
sigtab<- read.csv("sigtab_CB_WLE_enriched.depleted.inCB.csv")

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "P52", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

View(sigtab)
sigtab_thres<- sigtab %>% filter(log2FoldChange>=10|log2FoldChange<=-10)

plot<- ggplot(sigtab_thres, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  #scale_x_reordered() +
  geom_point(size=3) + scale_color_manual(values=as.vector(P52))+ 
  ggtitle("DeSeq analysis")+
  theme(plot.title = element_text(hjust = 0.5, size=18))+
  theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_text(size=18))+
  theme(axis.text.y=element_text(size=13), axis.title.y=element_text(size=18)) + 
  theme(legend.title = element_text(size=18))+ theme(legend.text = element_text(size=13)) +
  guides(col=guide_legend(ncol=1)) +theme(strip.text.x = element_text(size = 30)) +geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  ylab("Log2FoldChange")+coord_flip() 
plot

#1171 taxa at genus level are enriched in CB compared to CB. Thats too many to plot
ggsave(filename = "Deseq2_CB_WLE_enriched_depleted_inCB.tiff", plot = plot,
       width = 90,
       height = 100, units = c("cm"),
       dpi = 300)
ggsave(filename = "Deseq2_CB_WLE_enriched_depleted_inCB_zoomed.tiff", plot = plot,
       width = 40,
       height = 100, units = c("cm"),
       dpi = 300)

##64 taxa at genus level when log 2 fold change threshold set to -10 to +10
ggsave(filename = "Deseq2_CB_WLE_enriched_depleted_inCB_sigtab_thres.tiff", plot = plot,
       width = 25,
       height = 25, units = c("cm"),
       dpi = 300)

##load indicator species analysis package

install.packages("indicspecies")
library(indicspecies)

###########

##phyloseq object
library(phyloseq)
library(tidyverse)
otu_final = read.csv("OTU_clean_noneg_rarefied20k.csv", sep=",", row.names=1) ##USING UNRAREFIED DATA, CAN ALSO CHECK WITH RAREFIED 

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

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 40216 taxa and 199 samples ]
sample_data() Sample Data:       [ 199 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 40216 taxa by 8 taxonomic ranks ]

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

write.csv(ind_region, "indicator_region_CB_WLE.csv")



###################################
###indicator for transect site in CB 
####################################

##load indicator species analysis package

install.packages("indicspecies")
library(indicspecies)

###########

##phyloseq object
library(phyloseq)
library(tidyverse)
otu = read.csv("OTU_CB_rarefied20k.csv", sep=",", row.names=1) ##USING RAREFIED DATA here unlike unrarefied data used for site specific indicators  

otu_final = otu[rowSums(otu[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) ##24964 otus


tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)

metadata = read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1) 


OTU = otu_table(otu_final, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 24964 taxa and 117 samples ]
sample_data() Sample Data:       [ 117 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 24964 taxa by 8 taxonomic ranks ]

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

write.csv(ind_transectCB_subset, "ind_CB_wetland_unique.csv")



library(tidyverse)

write.csv(ind_transectCB, "indicator_transect_CB.csv")

##creating a dataframe for CB only indicators and passing on as a join to the WLE count table. 

##need to use unrarefied count tables for WLE here for the join

ind_transect_CB<- read.csv("indicator_transect_CB.csv", row.names=1)
View(ind_transect_CB)
ind_transect_CB_wetland<- ind_transect_CB %>% subset(s.wetland==1)
View(ind_transect_CB_wetland) ##2803 ASVs

merge_cbind_wetland<- merge(ind_transect_CB_wetland, otu_final, by=0)
View(merge_cbind_wetland)

write.csv(merge_cbind_wetland, "CB_wetlandindAll_count.csv")

##remove extra columns 

merge_cbind_wetland <- as.data.frame(merge_cbind_wetland)
rnames <-merge_cbind_wetland[,1] 
rownames(merge_cbind_wetland) <- rnames  # assign row names
View(merge_cbind_wetland)
merge_cbind_wetland_new <- merge_cbind_wetland %>% select(!c(Row.names, OTU, s.upland, s.transition, s.wetland,index, stat, sequencefeature, p.value, Domain, Phylum, Class, Order, Family, Genus, Species)) ##add other columns
View(merge_cbind_wetland_new)

write.csv(merge_cbind_wetland_new, "indicatorCB_wetland_count_final.csv")

##get these otus from WLE to assess treatment variation 
##merge to phyloseq object

##max_standardize 
max_standardize_CBwetland <-decostand(merge_cbind_wetland_new, method = "max", MARGIN = 1) ##from Jackson's phil trans paper
View(max_standardize_CBwetland)

tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("combined_metadataCB.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(max_standardize_CBwetland, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 2803 taxa and 116 samples ]
sample_data() Sample Data:       [ 116 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 2803 taxa by 8 taxonomic ranks ]


phyloseq_otu <- phyloseq_merged %>%
  #tax_glom(taxrank = "Class")   %>% 
  #transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()
View(phyloseq_otu)

##combining class and otu information

phyloseq_otu$ASVgenus <- paste(phyloseq_otu$Genus, phyloseq_otu$OTU) ##trying with genus names here
View(phyloseq_otu)




palette_new72 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
                  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
                  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "firebrick2",
                  "darkslategrey", "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4", "lightgreen", "magenta",
                  "lightgoldenrod1", "lightgoldenrod4", "lightpink3", "darkgoldenrod", "dodgerblue", "darkblue",
                  "hotpink", "hotpink4","khaki4","ivory3","ivory4","greenyellow","gold","aquamarine","coral4",
                  "coral", "sienna4", "springgreen", "rosybrown4","red")

install.packages("ggh4x")
library(ggh4x)
phyloseq_otu_new<- phyloseq_otu %>% mutate(transect1= ifelse(transect=="WC", "wetland", transect))
View(phyloseq_otu_new)

phyloseq_otu_new$transect1 <- factor(phyloseq_otu_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)
plot<- ggplot(phyloseq_otu_new, aes(x= transect1, y = Abundance, fill =transect1)) + 
  facet_grid2(site~ASVclass, scales="free", independent = "x", labeller = label_wrap_gen(width=25)) + theme(strip.text.x = element_text(size = 23, face="bold"))+ 
  theme(strip.text.y = element_text(size = 28, face="bold"))+
  geom_boxplot() +
  scale_fill_manual(values=as.vector(palette_new72))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=28,face="bold"),
                                     axis.title.y=element_text(size=28, face="bold")) + 
  theme(legend.title = element_text(size=28))+ theme(legend.text = element_text(size=28))+ guides(fill=guide_legend(ncol=1))
plot



ggsave(filename = "ASV_CBind_wetland.jpg", plot = plot,
       width = 120,
       height = 40, units = c("cm"),
       dpi = 300)


##plot the most abundant ones only 

merge_cbind_wetland_new<- read.csv("indicatorCB_wetland_count_final.csv", row.names=1)
View(merge_cbind_wetland_new)
merge_cbind_wetland_new<-merge_cbind_wetland_new[order(rowSums(merge_cbind_wetland_new), decreasing = TRUE),]
View(merge_cbind_wetland_new)
merge_cbind_wetland_new_20 <- merge_cbind_wetland_new[c(1:20),] ##select top 20 most abundant Classes
View(merge_cbind_wetland_new_20)

write.csv(merge_cbind_wetland_new_20, "CBind_wetland_top20.csv")

##max_standardize 
max_standardize_CBwetland <-decostand(merge_cbind_wetland_new_20, method = "max", MARGIN = 1) ##from Jackson's phil trans paper
View(max_standardize_CBwetland)

tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("combined_metadataCB.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(max_standardize_CBwetland, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 20 taxa and 116 samples ]
sample_data() Sample Data:       [ 116 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 20 taxa by 8 taxonomic ranks ]


phyloseq_otu <- phyloseq_merged %>%
  #tax_glom(taxrank = "Class")   %>% 
  #transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()
View(phyloseq_otu)

##combining class and otu information

phyloseq_otu$ASVgenus <- paste(phyloseq_otu$Genus, phyloseq_otu$OTU)
View(phyloseq_otu)




palette_new72 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
                  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
                  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "firebrick2",
                  "darkslategrey", "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4", "lightgreen", "magenta",
                  "lightgoldenrod1", "lightgoldenrod4", "lightpink3", "darkgoldenrod", "dodgerblue", "darkblue",
                  "hotpink", "hotpink4","khaki4","ivory3","ivory4","greenyellow","gold","aquamarine","coral4",
                  "coral", "sienna4", "springgreen", "rosybrown4","red")

install.packages("ggh4x")
library(ggh4x)
phyloseq_otu_new<- phyloseq_otu %>% mutate(transect1= ifelse(transect=="WC", "wetland", transect))
View(phyloseq_otu_new)

phyloseq_otu_new$transect1 <- factor(phyloseq_otu_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)
plot<- ggplot(phyloseq_otu_new, aes(x= transect1, y = Abundance, fill =transect1)) + 
  facet_grid2(site~ASVgenus, scales="free", independent = "x", labeller = label_wrap_gen(width=25)) + theme(strip.text.x = element_text(size = 18, face="bold"))+ 
  theme(strip.text.y = element_text(size = 28, face="bold"))+
  geom_boxplot() +
  scale_fill_manual(values=as.vector(palette_new72))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=28,face="bold"),
                                     axis.title.y=element_text(size=28, face="bold")) + 
  theme(legend.title = element_text(size=28))+ theme(legend.text = element_text(size=15))+ guides(fill=guide_legend(ncol=1))
plot



ggsave(filename = "ASV_CBind_wetland_top20.jpg", plot = plot,
       width =48,
       height = 10, units = c("in"),
       dpi = 300)


###############################
###plotting clade with colors of treatments
###############################

##Chesapeake Bay

library(phyloseq)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(vegan)

indicator_CB<- read.csv("indicator_transect_CB.csv" , header=TRUE)
View(indicator_CB)
indicator_filter<- indicator_CB %>% filter(s.upland==1 & s.transition==0 & s.wetland==0 | s.upland==0 & s.transition==1 & s.wetland==0|s.upland==0 & s.transition==0 & s.wetland==1)
View(indicator_filter)
indicator_filter$indicator_type<-  ifelse(indicator_filter$s.upland==1 & indicator_filter$s.transition==0 & indicator_filter$s.wetland==0, "upland_indicator", ifelse(indicator_filter$s.upland==0 & indicator_filter$s.transition==1 & indicator_filter$s.wetland==0, "transition_indicator", ifelse(indicator_filter$s.upland==0 & indicator_filter$s.transition==0 & indicator_filter$s.wetland==1, "wetland_indicator", NA)))


##instead of merging melt_simple with indicator taxa, lets merge the OTU table in phyloseq_Class object with indicator taxa first to select the indicator taxa, then order the table as decreasing rowsums. Pick the top 20 indicators. Then use psmelt for long format and use max standardization to account for relabund



otu_final<- read.csv("OTU_CB_rarefied20k.csv", sep=",", row.names=1)
View(otu_final)


otu_final1 = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final1) 

OTU <- rownames(otu_final1)
#rownames(d) <- NULL
otu_final1 <- cbind(OTU,otu_final1)


View(otu_final1)
join_abs_abun<- left_join(indicator_filter, otu_final1, by = "OTU")
View(join_abs_abun)
rnames <-join_abs_abun[,1] 
rownames(join_abs_abun) <- rnames  # assign row names

join_abs_abun1<- join_abs_abun %>% select(!c(X, OTU, s.upland, s.transition, s.wetland,index, stat, Confidence, p.value, Kingdom, Phylum, Class, Order, Family, Genus, Species, indicator_type)) ##add other columns
View(join_abs_abun1)

join_abs_abun2<-join_abs_abun1[order(rowSums(join_abs_abun1), decreasing = TRUE),]
View(join_abs_abun2)
join_abs_abun2_20 <- join_abs_abun2[c(1:20),] ##select top 20 most abundant Classes
View(join_abs_abun2_20)

##DO MAX STANDARDIZATION HERE decostand vegan
library(vegan)
max_standardize <-decostand(join_abs_abun2_20, method = "max", MARGIN = 1) ##from Jackson's phil trans paper
View(max_standardize)

##read new phyloseq object with the new joinabsabun object

OTU = otu_table(max_standardize, taxa_are_rows = TRUE)

tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)
TAX = tax_table(tax)


metadata = read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1) 
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged



phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 20 taxa and 117 samples ]
sample_data() Sample Data:       [ 117 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 20 taxa by 8 taxonomic ranks ]



melt_simple_new_1 <- phyloseq_merged %>% #transform_sample_counts(function(x) {x/sum(x)} )%>% ##dont transfrom to relabun yet, we need max standardization
  psmelt() %>%
  select(OTU, value=Abundance, transect, site, horizon, tree_number)
View(melt_simple_new)

melt_simple_new<- melt_simple_new_1 %>% mutate(transect1= ifelse(transect=="wetland-center", "wetland", transect))
join_data_new<- left_join(melt_simple_new, indicator_filter, by="OTU")%>% select (OTU, transect1, site, horizon, tree_number, indicator_type)

##reducing taxa to only indicators for purpose of plotting
#join_data_new<- right_join(melt_simple, indicator_filter, by="OTU") 

View(join_data_new)

library("ape")
set.seed(1)
random_tree = rtree(ntaxa(phyloseq_merged), rooted=TRUE, tip.label=taxa_names(phyloseq_merged))
plot(random_tree)
phyloseq_tree = merge_phyloseq(phyloseq_merged, random_tree)
phyloseq_tree


phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 20 taxa and 117 samples ]
sample_data() Sample Data:       [ 117 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 20 taxa by 8 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 20 tips and 19 internal nodes ]


##saving node and phylum inofrnation together in a dataframe

phylo<- phy_tree(phyloseq_tree)
phylo

write.tree(phylo, file="phylo.tre")

phylo_tibble<- as_tibble(phylo)
View(phylo_tibble)

phylo_tibble<- phylo_tibble %>% mutate(OTU=label)
View(phylo_tibble)

merge_phylo<- phylo_tibble %>% left_join(indicator_filter, by="OTU")
View(merge_phylo)

write.csv(merge_phylo, "merge_phylo_node_cb_transect_ind.csv")

write.csv(merge_phylo, "merge_phylo_node_all_cb_transect_ind.csv")

#merge_phylo_new<- read.csv("merge_phylo_node_all_cb_transect_ind.csv", check.names = F) #deleted the 19 internal nodes and modified as per the package github code
#merge_phylo_new_tbl<- as.tibble(merge_phylo_new)

#phylo_tree<- read.tree("phylo.tre")

#p <- ggtree(phylo_tree) + geom_tiplab2(aes(label=node), hjust=-.3)

#merge_phylo_new$id<-nodeid(phylo_tree,merge_phylo_new$id) #gives error
#merge_phylo$node<- nodeid(phylo_tree, merge_phylo$node) ## node is becoming NA - not working
#merge_phylo_new$id<- nodeid(phylo_tree,merge_phylo_new$id[nchar(merge_phylo_new$id)>=1])
#merge_phylo_new$transect<-factor(merge_phylo_new$transect,levels=c("upland", "transition", "wetland")) #gives error

merge_phylo$indicator_type<-factor(merge_phylo$indicator_type,levels=c("upland_indicator", "transition_indicator", "wetland_indicator"))

p <- ggtree(phylo, layout="fan", open.angle=10)#+ geom_text(aes(label=node), hjust=-.3)

#geom_tippoint(mapping=aes(color=Phylum), 
#  size=0.1,
# show.legend=FALSE)
p<- p %<+% join_data_new + geom_tippoint(pch=16, aes(col=indicator_type)) +scale_color_manual(values=c("#b2df8a","#33a02c","#fb9a99"))
p
p <- rotate_tree(p, -90)
p

#ggtree(phyloseq_tree) + geom_text(aes(label=node), hjust=-.3)


##if we remove NAs in join_data, the below code fails because now when the abundance is mapped to the tree there are missing ASVs with no "val"

##i pruned the to only indicator taxa in the phyloseq tree object as well as the join_data 
##the error with object "val" not being found was corrected by changing the name "val" to "value" when making melt_simple object, and subsequently creating a new join_data. not sure why it did that
##seems liek the ggtree object has matching variables to join_data namely indicator_type and "val".
##this is most likley why the "val" object was not being recognized in the below code.
##i figured it out: the tree has matching variables as it used join_data to color by indicator taxa above in the ggtree code
##so the variables in join data are doubling up when the below code is executed and it cannot recognize "val"
## changing val to value in melt_simple object and making a new joindata object, while keeping the tree the same with older join_data object, reduced one problem but the variable "indicator_taxa" is still doubled up, so is the other variables in join_data


join_data_new<- left_join(melt_simple_new, indicator_filter, by="OTU")%>% select (OTU, value, transect1, site, horizon, tree_number, indicator_type, Phylum, Class)

library(ggtreeExtra)
library(ggtree) 
library(treeio) 
library(tidytree) 
library(ggnewscale) 
library(ggpmisc) 
library(ggplot2) 
library(ggimage) 
library(patchwork)


if(!require("tidytree", quietly=TRUE))
  install_version("tidytree", version = "0.4.1")
library(tidytree)

p <- p +
  geom_fruit(
    data=join_data_new,
    geom=geom_boxplot,
    mapping = aes(
      y=OTU,
      x=value,
      #group=transect, ##change to site and others and see the best fit 
      fill=Phylum,
    ),
    size=0.5,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 2,
    ),
    grid.params=list()
  )  

p

write.csv(merge_phylo, "merge_phylo_cb_all_ind.csv")
merge_phylo<- read.csv("merge_phylo_cb_all_ind.csv")

##remove extra rows from phylo tibble instead of saving as csv

merge_phylo_20<- merge_phylo %>% filter(label!="NA")
View(merge_phylo_20)

#merge_phylo_20_tbl<- as_tibble(merge_phylo_20)
#merge_phylo_treedata<-as.treedata(merge_phylo_new)
#merge_phylo_treedata
#merge_phylo_treedata_tibble<- as_tibble(merge_phylo_treedata)

p1<- p+ 
  geom_cladelab(data=merge_phylo_20,
                mapping=aes(node=node,label=transect1,colour=transect1),
                textcolour=NA, barsize=2, extend=0.4, offset=1, align=TRUE)+ 
  scale_colour_manual(name="Transect", values=c("#b2df8a","#33a02c","#fb9a99", "#e31a1c","#EACB47","#6a3d9a"), 
                      guide=guide_legend(keywidth=0.3, keyheight=0.3, order=2, override.aes=list(size=1.5,alpha=1)))
p1

p2 <- p1 +
  scale_fill_discrete(
    name="Phylum",
    guide=guide_legend(keywidth=0.8, keyheight=0.8, ncol=2)
  ) +
  theme(
    legend.title=element_text(size=9), 
    legend.text=element_text(size=7) 
  ) 
p2



p3<- p2 %<+%merge_phylo_20 + geom_tiplab2(aes(label=Phylum), align=T, linetype=NA, 
                     size=2, offset=2.5, hjust=0.5) 
p3
  
  
  

  

##NEED TO ADD PHYLUM NAMES TO TREE
ggsave(filename = "tree_indicator_taxa_revised_cb_fullversion2.tiff", plot = p3,
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)

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
otu = read.csv("OTU_WLE_rarefied20k.csv", sep=",", row.names=1) ##USING RAREFIED DATA here unlike unrarefied data used for site specific indicators  

otu_final = otu[rowSums(otu[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) ##19262 otus


tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)

metadata = read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1) 


OTU = otu_table(otu_final, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19262 taxa and 82 samples ]
sample_data() Sample Data:       [ 82 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 19262 taxa by 8 taxonomic ranks ]

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
View(ind_transectWLE_subset) ##1292 Taxa uniquely associated with wetland

write.csv(ind_transectWLE_subset, "ind_WLE_wetland_unique.csv")



library(tidyverse)

write.csv(ind_transectWLE, "indicator_transect_WLE.csv")





##plot clade/tree WLE 

library(phyloseq)
library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(vegan)

indicator_WLE<- read.csv("indicator_transect_WLE.csv" , header=TRUE)
View(indicator_WLE)
indicator_filter<- indicator_WLE %>% filter(s.upland==1 & s.transition==0 & s.wetland==0 | s.upland==0 & s.transition==1 & s.wetland==0|s.upland==0 & s.transition==0 & s.wetland==1)
View(indicator_filter)
indicator_filter$indicator_type<-  ifelse(indicator_filter$s.upland==1 & indicator_filter$s.transition==0 & indicator_filter$s.wetland==0, "upland_indicator", ifelse(indicator_filter$s.upland==0 & indicator_filter$s.transition==1 & indicator_filter$s.wetland==0, "transition_indicator", ifelse(indicator_filter$s.upland==0 & indicator_filter$s.transition==0 & indicator_filter$s.wetland==1, "wetland_indicator", NA)))


##instead of merging melt_simple with indicator taxa, lets merge the OTU table in phyloseq_Class object with indicator taxa first to select the indicator taxa, then order the table as decreasing rowsums. Pick the top 20 indicators. Then use psmelt for long format and use max standardization to account for relabund



otu_final<- read.csv("OTU_WLE_rarefied20k.csv", sep=",", row.names=1)
View(otu_final)


otu_final1 = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final1) ##19262 OTUS

OTU <- rownames(otu_final1)
#rownames(d) <- NULL
otu_final1 <- cbind(OTU,otu_final1)


View(otu_final1)
join_abs_abun<- left_join(indicator_filter, otu_final1, by = "OTU")
View(join_abs_abun)
rnames <-join_abs_abun[,1] 
rownames(join_abs_abun) <- rnames  # assign row names

join_abs_abun1<- join_abs_abun %>% select(!c(X, OTU, s.upland, s.transition, s.wetland,index, stat, Confidence, p.value, Kingdom, Phylum, Class, Order, Family, Genus, Species, indicator_type)) ##add other columns
View(join_abs_abun1)

join_abs_abun2<-join_abs_abun1[order(rowSums(join_abs_abun1), decreasing = TRUE),]
View(join_abs_abun2)
join_abs_abun2_20 <- join_abs_abun2[c(1:20),] ##select top 20 most abundant Classes
View(join_abs_abun2_20)

##DO MAX STANDARDIZATION HERE decostand vegan
library(vegan)
max_standardize <-decostand(join_abs_abun2_20, method = "max", MARGIN = 1) ##from Jackson's phil trans paper
View(max_standardize)

##read new phyloseq object with the new joinabsabun object

OTU = otu_table(max_standardize, taxa_are_rows = TRUE)

tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)
TAX = tax_table(tax)


metadata = read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1) 
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged



phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 20 taxa and 82 samples ]
sample_data() Sample Data:       [ 82 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 20 taxa by 8 taxonomic ranks ]


melt_simple_new_1 <- phyloseq_merged %>% #transform_sample_counts(function(x) {x/sum(x)} )%>% ##dont transfrom to relabun yet, we need max standardization
  psmelt() %>%
  select(OTU, value=Abundance, transect, site, horizon, tree_number)


melt_simple_new<- melt_simple_new_1 %>% mutate(transect1= ifelse(transect=="wetland-center"|transect=="wetland-transition-edge", "wetland", transect))
join_data_new<- left_join(melt_simple_new, indicator_filter, by="OTU")%>% select (OTU, transect1, site, horizon, tree_number, indicator_type)

##reducing taxa to only indicators for purpose of plotting
#join_data_new<- right_join(melt_simple, indicator_filter, by="OTU") 

View(join_data_new)

library("ape")
set.seed(1)
random_tree = rtree(ntaxa(phyloseq_merged), rooted=TRUE, tip.label=taxa_names(phyloseq_merged))
plot(random_tree)
phyloseq_tree = merge_phyloseq(phyloseq_merged, random_tree)
phyloseq_tree


phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 20 taxa and 82 samples ]
sample_data() Sample Data:       [ 82 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 20 taxa by 8 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 20 tips and 19 internal nodes ]


##saving node and phylum inofrnation together in a dataframe

phylo<- phy_tree(phyloseq_tree)
phylo

#write.tree(phylo, file="phylo.tre")

phylo_tibble<- as_tibble(phylo)
View(phylo_tibble)

phylo_tibble<- phylo_tibble %>% mutate(OTU=label)
View(phylo_tibble)

merge_phylo<- phylo_tibble %>% left_join(indicator_filter, by="OTU")
View(merge_phylo)



write.csv(merge_phylo, "merge_phylo_node_all_wle_transect_ind.csv")

#merge_phylo_new<- read.csv("merge_phylo_node_all_cb_transect_ind.csv", check.names = F) #deleted the 19 internal nodes and modified as per the package github code
#merge_phylo_new_tbl<- as.tibble(merge_phylo_new)

#phylo_tree<- read.tree("phylo.tre")

#p <- ggtree(phylo_tree) + geom_tiplab2(aes(label=node), hjust=-.3)

#merge_phylo_new$id<-nodeid(phylo_tree,merge_phylo_new$id) #gives error
#merge_phylo$node<- nodeid(phylo_tree, merge_phylo$node) ## node is becoming NA - not working
#merge_phylo_new$id<- nodeid(phylo_tree,merge_phylo_new$id[nchar(merge_phylo_new$id)>=1])
#merge_phylo_new$transect<-factor(merge_phylo_new$transect,levels=c("upland", "transition", "wetland")) #gives error

merge_phylo$indicator_type<-factor(merge_phylo$indicator_type,levels=c("upland_indicator", "transition_indicator", "wetland_indicator"))

p <- ggtree(phylo, layout="fan", open.angle=10)#+ geom_text(aes(label=node), hjust=-.3)

#geom_tippoint(mapping=aes(color=Phylum), 
#  size=0.1,
# show.legend=FALSE)
p<- p %<+% join_data_new + geom_tippoint(pch=16, aes(col=indicator_type)) +scale_color_manual(values=c("#b2df8a","#33a02c","#fb9a99"))
p
p <- rotate_tree(p, -90)
p

#ggtree(phyloseq_tree) + geom_text(aes(label=node), hjust=-.3)


##if we remove NAs in join_data, the below code fails because now when the abundance is mapped to the tree there are missing ASVs with no "val"

##i pruned the to only indicator taxa in the phyloseq tree object as well as the join_data 
##the error with object "val" not being found was corrected by changing the name "val" to "value" when making melt_simple object, and subsequently creating a new join_data. not sure why it did that
##seems liek the ggtree object has matching variables to join_data namely indicator_type and "val".
##this is most likley why the "val" object was not being recognized in the below code.
##i figured it out: the tree has matching variables as it used join_data to color by indicator taxa above in the ggtree code
##so the variables in join data are doubling up when the below code is executed and it cannot recognize "val"
## changing val to value in melt_simple object and making a new joindata object, while keeping the tree the same with older join_data object, reduced one problem but the variable "indicator_taxa" is still doubled up, so is the other variables in join_data


join_data_new<- left_join(melt_simple_new, indicator_filter, by="OTU")%>% select (OTU, value, transect1, site, horizon, tree_number, indicator_type, Phylum, Class)

library(ggtreeExtra)
library(ggtree) 
library(treeio) 
library(tidytree) 
library(ggnewscale) 
library(ggpmisc) 
library(ggplot2) 
library(ggimage) 
library(patchwork)


if(!require("tidytree", quietly=TRUE))
  install_version("tidytree", version = "0.4.1")
library(tidytree)

p <- p +
  geom_fruit(
    data=join_data_new,
    geom=geom_boxplot,
    mapping = aes(
      y=OTU,
      x=value,
      #group=transect, ##change to site and others and see the best fit 
      fill=Phylum,
    ),
    size=0.5,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 2,
    ),
    grid.params=list()
  )  

p

write.csv(merge_phylo, "merge_phylo_wle_all_ind.csv")
merge_phylo<- read.csv("merge_phylo_wle_all_ind.csv")

##remove extra rows from phylo tibble instead of saving as csv

merge_phylo_20<- merge_phylo %>% filter(label!="NA")
View(merge_phylo_20)

#merge_phylo_20_tbl<- as_tibble(merge_phylo_20)
#merge_phylo_treedata<-as.treedata(merge_phylo_new)
#merge_phylo_treedata
#merge_phylo_treedata_tibble<- as_tibble(merge_phylo_treedata)

p1<- p+ 
  geom_cladelab(data=merge_phylo_20,
                mapping=aes(node=node,label=transect1,colour=transect1),
                textcolour=NA, barsize=2, extend=0.4, offset=1, align=TRUE)+ 
  scale_colour_manual(name="Transect", values=c("#b2df8a","#33a02c","#fb9a99", "#e31a1c","#EACB47","#6a3d9a"), 
                      guide=guide_legend(keywidth=0.3, keyheight=0.3, order=2, override.aes=list(size=1.5,alpha=1)))
p1

p2 <- p1 +
  scale_fill_discrete(
    name="Phylum",
    guide=guide_legend(keywidth=0.8, keyheight=0.8, ncol=2)
  ) +
  theme(
    legend.title=element_text(size=9), 
    legend.text=element_text(size=7) 
  ) 
p2



p3<- p2 %<+%merge_phylo_20 + geom_tiplab2(aes(label=Phylum), align=T, linetype=NA, 
                                          size=2, offset=2.5, hjust=0.5) 
p3






##NEED TO ADD PHYLUM NAMES TO TREE
ggsave(filename = "tree_indicator_taxa_revised_wle_fullversion2.tiff", plot = p3,
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)







###abundance occupancy relationship 
##Chesapeake bay
library(phyloseq)
library(tidyverse)

otu_final<- read.csv("OTU_CB_rarefied20k.csv", sep=",", row.names=1)
View(otu_final)


otu_final1 = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final1) 


OTU = otu_table(otu_final1, taxa_are_rows = TRUE)

tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)
TAX = tax_table(tax)


metadata = read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1) 
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 24964 taxa and 117 samples ]
sample_data() Sample Data:       [ 117 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 24964 taxa by 8 taxonomic ranks ]

##the psmelt and agglom is reducing the read count per sample slightly, total reads is not 117*20000=2340000 but 2327516
physeq_phylum_subset <- phyloseq_merged %>%  tax_glom(taxrank = "Phylum") %>% psmelt() %>% #filter(Abundance > 0) ##since the reads numbers are being lost I am removing the  filter option to see what happens
  arrange(Phylum)
write.csv(physeq_phylum_subset, 'physeq_phylum_subset_abs_nozero_cb.csv')
write.csv(physeq_phylum_subset, 'physeq_phylum_subset_abs_cb.csv') ##no change, same as above, which makes sense

##reducing the phylum to plot only the most abundant ones, this below code is wrong,because everything here is absolute abundance, there is no fraction. so better to sort the final relabund table in decreasing order of relabund and plot the highest abundant ones that way

physeq_phylum_subset <- phyloseq_merged %>%  tax_glom(taxrank = "Phylum") %>% psmelt() %>% filter(Abundance > 0.2) %>%
  arrange(Phylum)
write.csv(physeq_phylum_subset, 'physeq_phylum_subset_abs_cb_0.2filter.csv')

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5"
)


library(tidyverse)
data<- read.csv("physeq_phylum_subset_abs_nozero_cb.csv")##removed excess columns to clean up and be able to format to wide format
data<- read.csv("physeq_phylum_subset_abs_cb.csv")##removed excess columns to clean up and be able to format to wide format
data<- read.csv("physeq_phylum_subset_abs_cb_0.2filter.csv")##removed excess columns to clean up and be able to format to wide format


df.wide<-pivot_wider(data, names_from = SampleID, values_from = Abundance, values_fill = 0)
write.csv(df.wide, 'abs_nozero_wide_cb.csv')
write.csv(df.wide, 'abs_wide_cb.csv')
write.csv(df.wide, 'abs_wide_cb_0.2filter.csv')

otu_abs<- read.csv("physeq_phylum_subset_abs_cb_0.2filter.csv") #change this file as needed
otu_abs.PA<- mutate(otu_abs, PA = if_else(Abundance>0, 1, 0))
write.csv(otu_abs.PA, 'nozero_PA_cb.csv') ##remove all columns except phylum, PA, Sample ID
write.csv(otu_abs.PA, 'PA_cb_0.2filter.csv') ##remove all columns except phylum, PA, Sample ID

otu_PA<-read.csv("nozero_PA_cb.csv")
otu_PA<-read.csv("PA_cb_0.2filter.csv")
df.wide.PA<-pivot_wider(otu_PA, names_from = SampleID, values_from = PA, values_fill = 0)
write.csv(df.wide.PA, 'abs_nozero_wide_PA_cb.csv') #remove ID numbers from first column
write.csv(df.wide.PA, 'abs_wide_PA_cb_0.2filter.csv') #remove ID numbers from first column


##occupancy
##calculating occupancy of all taxa
otu_PA<-read.csv("abs_nozero_wide_PA_cb.csv", row.names = 1)
otu_PA<-read.csv("abs_wide_PA_cb_0.2filter.csv", row.names = 1)
otu_PA.df<-as.data.frame(otu_PA)
#checking if the factor is class or numeric
str(otu_PA)
#need to make numeric
otu_PA_transform<-transform(otu_PA.df, occupancy=rowSums(sapply(otu_PA.df, as.numeric))/length(colnames(otu_PA.df)))
#View(otu_PA_transform)          
write.csv(otu_PA_transform,"otu_PA_occupancy_cb.csv")
write.csv(otu_PA_transform,"otu_PA_occupancy_cb_0.2filter.csv")

otu_PA_transform<- read.csv("otu_PA_occupancy_cb.csv")
otu_PA_transform<- read.csv("otu_PA_occupancy_cb_0.2filter.csv")
dim(otu_PA_transform)

##filtering rel abund to greater than 0.5%

# edit abs_nozero_wide_cb.csv by adding the rel abund column row-wise
#abs_nozero_wide_cb<- read.csv("abs_nozero_wide_cb.csv")


data<- read.csv("physeq_phylum_subset_abs_nozero_cb.csv")##removed excess columns to clean up and be able to format to wide format
data<- read.csv("physeq_phylum_subset_abs_cb.csv")
data<- read.csv("physeq_phylum_subset_abs_cb_0.2filter.csv")

df.wide<-pivot_wider(data, names_from = SampleID, values_from = Abundance, values_fill = 0)
View(df.wide)
class(df.wide)
df.wide<- as.data.frame(df.wide)
rnames <-df.wide[,1] 
rownames(df.wide) <- rnames  # assign row names

df.wide1<- df.wide %>% select(!c(Phylum)) ##add other columns
View(df.wide1)

df.wide1.transform<- transform(df.wide1, sum=rowSums(sapply(df.wide1, as.numeric)))
View(df.wide1.transform)

df.wide1.transform.new<- transform(df.wide1.transform, relabund=rowSums(sapply(df.wide1.transform, as.numeric))/sum(sum))
View(df.wide1.transform.new)

write.csv(df.wide1.transform.new, "otu_RA_cb.csv")
write.csv(df.wide1.transform.new, "otu_RA_cb_0.2filter.csv")

##testing original read count of CB otu table 

otu_final2<- transform(otu_final1, rowsum=rowSums(sapply(otu_final1, as.numeric)))

sum(otu_final2$rowsum)


#read in abund_occ 0.2 filter file and arrange in deceasing order of relabund
abund_occ<- read.csv("occ_abund_0.2filter.csv")

abund_occ_new<-abund_occ[order(abund_occ$relabund, decreasing=TRUE), ]

abund_occ_new_20 <- abund_occ_new[c(1:20),] ##select top 20 most abundant Classes
View(abund_occ_new_20)


##ggplot
#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
library(ggplot2)
library(viridis)
library(rcartocolor)
library(pals)
#data<- read.csv("otu_RA_PA_cb.csv") ##obtained by combining occupancy and rel abundance data in csv file.
gg <- ggplot(abund_occ_new_20, aes(x=relabund, y=occupancy)) +
  geom_point(aes(col=Phylum), size=6) + 
  #geom_errorbar(aes(xmin = relabund - stderror, 
  #xmax = relabund + stderror), width = 0.2, size = 1, color="white")+
  #geom_text(aes(label=Phylum), hjust = 0, nudge_y = 0.25)+
  #scale_fill_manual(values="#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD","#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5") +
  #geom_smooth(method="loess", se=F) + 
  #xlim(c(0, 100)) + 
  #ylim(c(85, 100)) + 
  #scale_color_viridis(discrete=TRUE)+
  labs(
    y="Occupancy", 
    x="Mean Relative Abundance") +theme(legend.justification = c("right"))+theme_classic()+guides(color=guide_legend(ncol=1))

gg 

#gg+scale_color_carto_d(palette="Safe")


palette_new72 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
                  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
                  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "firebrick2",
                  "darkslategrey", "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4", "lightgreen", "magenta",
                  "lightgoldenrod1", "lightgoldenrod4", "lightpink3", "darkgoldenrod", "dodgerblue", "darkblue",
                  "hotpink", "hotpink4","khaki4","ivory3","ivory4","greenyellow","gold","aquamarine","coral4",
                  "coral", "sienna4", "springgreen", "rosybrown4","red")

plot<-gg+scale_color_manual(values=as.vector(palette_new72))


#plot<-gg+scale_color_manual(values=c("black","forestgreen", "red2", "orange", "cornflowerblue", 
#"magenta", "tan3", "darkblue", 
# "mediumorchid3","firebrick4",  "yellowgreen", "lightsalmon"))
plot

ggsave(plot, file="abundocc_cb_final.TIFF", dpi=300, height=7, width=7, unit=c("in"))


###not needed below code
data<- read.csv("occ_abund2.csv")
gg <- ggplot(data, aes(x=relabund, y=occupancy_percent)) + 
  geom_point(aes(col=Phylum, size=occupancy_percent)) + 
  #geom_smooth(method="loess", se=F) + 
  xlim(c(0, 100)) + 
  ylim(c(0, 100)) + 
  labs(
    y="Occupancy(%)", 
    x="Mean Relative Abundance")

plot(gg)

#############################
##testing whether the uncultured taxa are being removed in the agglomeration above since there are multiple uncultured taxa named similarly belonging to multiple phyla
##############################

taxonomy<- read.csv("taxonomy-dn-99.csv")
taxonomy_duplicated<-subset(taxonomy,duplicated(Species))
taxonomy_duplicated_species<- as.tibble(taxonomy_duplicated) %>% select (c('OTUID', 'Species')) %>% filter(Species!="")

taxonomy_know_duplicated<- unique(taxonomy_duplicated_species[c("Species")])
write.csv(taxonomy_know_duplicated, "taxonomy_know_duplicated_species.csv")

##merge the above taxonomy file with the original taxonomy dn 99 file to get the otus I need, then do rowsums 

taxonomy1<- read.csv("taxonomy-dn-99.csv")
taxonomy2<- read.csv("taxonomy_know_duplicated_species.csv")

merge_tax<- taxonomy2 %>% left_join(taxonomy1, by="Species")

otu_final<- read.csv("OTU_CB_rarefied20k.csv", sep=",", row.names=1) ##read CB otu table
otu_final1 = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final1) 

otu_cb_new<- otu_final1 %>% mutate(OTUID=rownames(otu_final1))

##merge otu cb table with the merge tax to keep those duplicated species OTUs

merge_otu_tax_cb<- inner_join(merge_tax, otu_cb_new, by="OTUID")

merge_otu_tax_cb_final<- merge_otu_tax_cb %>% select(!c('Species', 'Kingdom','Phylum','Class','Order', 'Family', 'Genus','Confidence', 'X'))

merge_otu_tax_cb_final<- as.data.frame(merge_otu_tax_cb_final)
rnames <-merge_otu_tax_cb_final[,1] 
rownames(merge_otu_tax_cb_final) <- rnames  # assign row names

merge_otu_tax_cb_final<- merge_otu_tax_cb_final %>% select(!c(OTUID)) ##add other columns
View(merge_otu_tax_cb_final)

merge_otu_tax_cb_final_transform<- transform(merge_otu_tax_cb_final, sum=rowSums(sapply(merge_otu_tax_cb_final, as.numeric)))
sum(merge_otu_tax_cb_final_transform$sum) ##1179081, does not match the depleted numbers in agglom. something else is off
##solution: this is possibly because we are removing all species duplicates, some duplicates may belong to the same phylum and those should be counted when agglomerating. only those
#which are duplicated but belong to different phyla are likely being removed

#############
#####abund-occupancy WLE
#############


library(phyloseq)
library(tidyverse)

otu_final<- read.csv("OTU_WLE_rarefied20k.csv", sep=",", row.names=1)
View(otu_final)


otu_final1 = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final1) ##19262 OTUs


OTU = otu_table(otu_final1, taxa_are_rows = TRUE)

tax = read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
tax = as.matrix(tax)
TAX = tax_table(tax)


metadata = read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1) 
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19262 taxa and 82 samples ]
sample_data() Sample Data:       [ 82 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 19262 taxa by 8 taxonomic ranks ]

##the psmelt and agglom is reducing the read count per sample slightly, total reads is not 117*20000=2340000 but 2327516
physeq_phylum_subset <- phyloseq_merged %>%  tax_glom(taxrank = "Phylum") %>% psmelt() %>% filter(Abundance > 0)%>% ##since the reads numbers are being lost I am removing the  filter option to see what happens
  arrange(Phylum)
write.csv(physeq_phylum_subset, 'physeq_phylum_subset_abs_nozero_wle.csv')


phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5"
)


library(tidyverse)
data<- read.csv("physeq_phylum_subset_abs_nozero_wle.csv")##removed excess columns to clean up and be able to format to wide format


df.wide<-pivot_wider(data, names_from = SampleID, values_from = Abundance, values_fill = 0)
write.csv(df.wide, 'abs_nozero_wide_wle.csv')

otu_abs<- read.csv("physeq_phylum_subset_abs_nozero_wle.csv") #change this file as needed
otu_abs.PA<- mutate(otu_abs, PA = if_else(Abundance>0, 1, 0))
write.csv(otu_abs.PA, 'nozero_PA_wle.csv') ##remove all columns except phylum, PA, Sample ID


otu_PA<-read.csv("nozero_PA_wle.csv")
df.wide.PA<-pivot_wider(otu_PA, names_from = SampleID, values_from = PA, values_fill = 0)
write.csv(df.wide.PA, 'abs_nozero_wide_PA_wle.csv') #remove ID numbers from first column



##occupancy
##calculating occupancy of all taxa
otu_PA<-read.csv("abs_nozero_wide_PA_wle.csv", row.names = 1)
otu_PA.df<-as.data.frame(otu_PA)
#checking if the factor is class or numeric
str(otu_PA)
#need to make numeric
otu_PA_transform<-transform(otu_PA.df, occupancy=rowSums(sapply(otu_PA.df, as.numeric))/length(colnames(otu_PA.df)))
#View(otu_PA_transform)          
write.csv(otu_PA_transform,"otu_PA_occupancy_wle.csv")

otu_PA_transform<- read.csv("otu_PA_occupancy_wle.csv")
dim(otu_PA_transform)

##filtering rel abund to greater than 0.5%

# edit abs_nozero_wide_cb.csv by adding the rel abund column row-wise
#abs_nozero_wide_cb<- read.csv("abs_nozero_wide_cb.csv")


data<- read.csv("physeq_phylum_subset_abs_nozero_wle.csv")##removed excess columns to clean up and be able to format to wide format

df.wide<-pivot_wider(data, names_from = SampleID, values_from = Abundance, values_fill = 0)
View(df.wide)
class(df.wide)
df.wide<- as.data.frame(df.wide)
rnames <-df.wide[,1] 
rownames(df.wide) <- rnames  # assign row names

df.wide1<- df.wide %>% select(!c(Phylum)) ##add other columns
View(df.wide1)

df.wide1.transform<- transform(df.wide1, sum=rowSums(sapply(df.wide1, as.numeric)))
View(df.wide1.transform)
sum(df.wide1.transform$sum) #1631248 reads which is less than 1640000



df.wide1.transform.new<- transform(df.wide1.transform, relabund=rowSums(sapply(df.wide1.transform, as.numeric))/sum(sum))
View(df.wide1.transform.new)

write.csv(df.wide1.transform.new, "otu_RA_wle.csv")


##testing original read count of WLE otu table 

otu_final2<- transform(otu_final1, rowsum=rowSums(sapply(otu_final1, as.numeric)))

sum(otu_final2$rowsum) ##1640000 reads 


#read in abund_occ file and arrange in deceasing order of relabund
abund_occ<- read.csv("occ_abund_wle.csv")

abund_occ_new<-abund_occ[order(abund_occ$relabund, decreasing=TRUE), ]

abund_occ_new_20 <- abund_occ_new[c(1:20),] ##select top 20 most abundant Classes
View(abund_occ_new_20)


##ggplot
#color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
library(ggplot2)
library(viridis)
library(rcartocolor)
library(pals)
#data<- read.csv("otu_RA_PA_cb.csv") ##obtained by combining occupancy and rel abundance data in csv file.
gg <- ggplot(abund_occ_new_20, aes(x=relabund, y=occupancy)) +
  geom_point(aes(col=Phylum), size=6) + 
  #geom_errorbar(aes(xmin = relabund - stderror, 
  #xmax = relabund + stderror), width = 0.2, size = 1, color="white")+
  #geom_text(aes(label=Phylum), hjust = 0, nudge_y = 0.25)+
  #scale_fill_manual(values="#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD","#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5") +
  #geom_smooth(method="loess", se=F) + 
  #xlim(c(0, 100)) + 
  #ylim(c(85, 100)) + 
  #scale_color_viridis(discrete=TRUE)+
  labs(
    y="Occupancy", 
    x="Mean Relative Abundance") +theme(legend.justification = c("right"))+theme_classic()+guides(color=guide_legend(ncol=1))

gg 

#gg+scale_color_carto_d(palette="Safe")


palette_new72 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
                  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
                  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "firebrick2",
                  "darkslategrey", "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4", "lightgreen", "magenta",
                  "lightgoldenrod1", "lightgoldenrod4", "lightpink3", "darkgoldenrod", "dodgerblue", "darkblue",
                  "hotpink", "hotpink4","khaki4","ivory3","ivory4","greenyellow","gold","aquamarine","coral4",
                  "coral", "sienna4", "springgreen", "rosybrown4","red")

plot<-gg+scale_color_manual(values=as.vector(palette_new72))


#plot<-gg+scale_color_manual(values=c("black","forestgreen", "red2", "orange", "cornflowerblue", 
#"magenta", "tan3", "darkblue", 
# "mediumorchid3","firebrick4",  "yellowgreen", "lightsalmon"))
plot

ggsave(plot, file="abundocc_wle_final.TIFF", dpi=300, height=7, width=7, unit=c("in"))


###not needed code below
data<- read.csv("occ_abund2.csv")
gg <- ggplot(data, aes(x=relabund, y=occupancy_percent)) + 
  geom_point(aes(col=Phylum, size=occupancy_percent)) + 
  #geom_smooth(method="loess", se=F) + 
  xlim(c(0, 100)) + 
  ylim(c(0, 100)) + 
  labs(
    y="Occupancy(%)", 
    x="Mean Relative Abundance")

plot(gg)


#combining the abund occupancy data between CB and WLE

abund_occ_wle<- read.csv("occ_abund_wle.csv")
abund_occ_cb<- read.csv("occ_abund_cb_0.2filter.csv")

abund_occ_new_wle<-abund_occ_wle[order(abund_occ_wle$relabund, decreasing=TRUE), ]

abund_occ_new_wle_20 <- abund_occ_new_wle[c(1:20),] ##select top 20 most abundant Classes
abund_occ_new_wle_20<- abund_occ_new_wle_20 %>% mutate(region="Lake Erie")
View(abund_occ_new_wle_20)

abund_occ_new_cb<-abund_occ_cb[order(abund_occ_cb$relabund, decreasing=TRUE), ]

abund_occ_new_cb_20 <- abund_occ_new_cb[c(1:20),] ##select top 20 most abundant Classes
abund_occ_new_cb_20<- abund_occ_new_cb_20 %>% mutate(region="Chesapeake Bay")
View(abund_occ_new_cb_20)

abund_occ_cb_wle_20<-rbind(abund_occ_new_wle_20, abund_occ_new_cb_20)

gg <- ggplot(abund_occ_cb_wle_20, aes(x=relabund, y=occupancy)) +
  facet_grid(~region)+
  geom_point(aes(col=Phylum), size=6) + 
  #geom_errorbar(aes(xmin = relabund - stderror, 
  #xmax = relabund + stderror), width = 0.2, size = 1, color="white")+
  #geom_text(aes(label=Phylum), hjust = 0, nudge_y = 0.25)+
  #scale_fill_manual(values="#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD","#AD6F3B", "#673770","#D14285", "#652926", "#C84248", "#8569D5") +
  #geom_smooth(method="loess", se=F) + 
  #xlim(c(0, 100)) + 
  #ylim(c(85, 100)) + 
  #scale_color_viridis(discrete=TRUE)+
  labs(
    y="Occupancy", 
    x="Mean Relative Abundance") +theme(legend.justification = c("right"))+theme_classic()+guides(color=guide_legend(ncol=2))

gg 

#gg+scale_color_carto_d(palette="Safe")


palette_new72 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
                  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
                  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "firebrick2",
                  "darkslategrey", "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4", "lightgreen", "magenta",
                  "lightgoldenrod1", "lightgoldenrod4", "lightpink3", "darkgoldenrod", "dodgerblue", "darkblue",
                  "hotpink", "hotpink4","khaki4","ivory3","ivory4","greenyellow","gold","aquamarine","coral4",
                  "coral", "sienna4", "springgreen", "rosybrown4","red")

plot<-gg+scale_color_manual(values=as.vector(palette_new72))


plot

ggsave(plot, file="abundocc_cb_wle_20phyla_final.TIFF", dpi=300, height=5, width=11, unit=c("in"))






