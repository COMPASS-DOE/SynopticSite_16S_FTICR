library(tidyverse)
library(reshape2)
data<- read.csv("COMPASS_CB_WLE_count.csv", sep=",",row.names=1)
View(data)
data_WLE<- data %>% select(!starts_with("s"))
View(data_WLE)


library(phyloseq)


tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("WLE_metadata.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(data_WLE, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 50421 taxa and 85 samples ]
sample_data() Sample Data:       [ 85 samples by 5 sample variables ]
tax_table()   Taxonomy Table:    [ 50421 taxa by 8 taxonomic ranks ]


sample_sum_df <- data.frame(sum = sample_sums(phyloseq_merged))
sum(sample_sum_df) ##426275 reads

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

smin <- min(sample_sums(phyloseq_merged))
smin # 79
smean <- mean(sample_sums(phyloseq_merged))
smean # 5015
smax <- max(sample_sums(phyloseq_merged))
smax # 43018


taxa_are_rows(phyloseq_merged)
#### Rarefaction curves

library(scales)
library(vegan)
library(dplyr)

mat<- as(t(otu_table(phyloseq_merged)), "matrix")
system.time(rarecurve(mat, step=1000, label = FALSE, col="cyan4")) 

##rarefaction

sample_df <- data.frame(sample_data(phyloseq_merged))
View(sample_df)
write.csv(sample_df, "samples_pre_rarefied_WLE.csv")

library(vegan)
#rarecurve(t(otu_table(phyloseq_merged)), step=50)
phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged, sample.size = 1000, rngseed = TRUE, trimOTUs=FALSE)

#19 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
  
 # COMPASS.Dec2021.003COMPASS.Dec2021.004COMPASS.Dec2021.006COMPASS.Dec2021.008COMPASS.Dec2021.010

sample_df <- data.frame(sample_data(phyloseq_rarefy))

View(sample_df)
write.csv(sample_df, "samples_rarefied_WLE.csv")

##use anti-join to select non-present samples in rarefied object
#This join is like eg: for df1-df2, it selects all rows from df1 that are not present in df2.
df1<- read.csv("samples_pre_rarefied_WLE.csv")
df2 <- read.csv("samples_rarefied_WLE.csv")
df= df1 %>% anti_join(df2,by="X")

View(df)

write.csv(df, "samples-lost-rarefaction-WLE.csv")

OTU2 = as(otu_table(phyloseq_rarefy), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, "OTU_rarefied1k.csv")

otu_final<- read.csv("OTU_rarefied1k.csv", row.names=1) 
View(otu_final)
any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #47164 otus with no counts
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) #3257 otus 

asv_wle<- read.csv("ASV_WLE.csv")
View(asv_wle)

asv_long<- asv_wle %>% pivot_longer(!X.OTU.ID, names_to="sample", values_to="abund")
View(asv_long)

samples_lost<- read.csv("samples-lost-rarefaction-WLE.csv")
View(samples_lost)

asv_long_mutate<- asv_long%>% mutate(X=X.OTU.ID)
View(asv_long_mutate)


samples_lost_mutate<- samples_lost %>% mutate(sample=X)
View(samples_lost_mutate)

merge_bysample<- anti_join(asv_long_mutate, samples_lost_mutate, by="sample")
View(merge_bysample)

merge_bysample_wide<- merge_bysample %>% pivot_wider(names_from="sample", values_from="abund", values_fill=0)
View(merge_bysample_wide)
merge_bysample_wide_nozero = merge_bysample_wide[rowSums(merge_bysample_wide[])>0, ,drop=FALSE]
merge_bysample_wide<- merge_bysample_wide %>% select(!X.OTU.ID)
merge_bysample_wide<- merge_bysample_wide %>% column_to_rownames(var="X")
merge_bysample_wide_nozero = merge_bysample_wide[rowSums(merge_bysample_wide[])>0, ,drop=FALSE]
View(merge_bysample_wide_nozero) ##this wil not match with 3257 ASVs as the samples have not been "rarefied" per se, so we see miore ASVs 3307 rather than 3257 


##lets try to do the phyloseq merged with the ASV table from scratch

data<- read.csv("ASV_WLE.csv", sep=",",row.names=1) ##should be rarefied one but works for testing purposes
View(data)


library(phyloseq)


tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("WLE_metadata.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(data, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

##Not working because the taxonomy table has updated asv calls produced after merging using dada2


##continuing with analysis using prior rarefied data

otu_final<- read.csv("OTU_rarefied1k.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #47164 otus with no counts
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final)

library(phyloseq)


tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("WLE_metadata.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 3257 taxa and 66 samples ]
sample_data() Sample Data:       [ 66 samples by 5 sample variables ]
tax_table()   Taxonomy Table:    [ 3257 taxa by 8 taxonomic ranks ]


#metadata<- read.csv("sample_metadata_dna_cdna.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones
#keep.samples <- as.vector(metadata$SampleID)
#keep.samples

#phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
#phyloseq_merged_final

##barplot class level
phyloseq_class <- phyloseq_merged %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Class) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."

View(phyloseq_class)


##Desulfobacteria class plot only

phyloseq_class_desulfo<- phyloseq_class %>% filter(Class==" c__Desulfobacteria")
View(phyloseq_class_desulfo)

library(ggh4x)
plot_class<- ggplot(phyloseq_class_desulfo, aes(x= site, y = Abundance, fill =Class)) + 
  facet_nested(.~transect+horizon, scales="free", independent = "x", labeller = label_wrap_gen(width=20)) + theme(strip.text.x = element_text(size = 14, face="bold"))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_boxplot() +
  scale_fill_manual(values=c("#39568CFF", "#55C667FF", "orange"))+
  theme(axis.title.x = element_text(size=17)) + theme(axis.text.x = element_text(size=15))+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=20, face="bold"),
                                     axis.title.y=element_text(size=20, face="bold")) + 
  ggtitle("Class Desulfobacteria")+
  theme(legend.title = element_blank())+ theme(legend.text = element_blank()) #guides(fill=guide_legend(ncol=1))
plot_class

ggsave(filename = "boxplot_desulfo_WLE_filtered0.01.tiff", plot = plot_class,
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)

##kruskal test for wetland desulfo

phyloseq_class_desulfo_wetland<- phyloseq_class_desulfo %>% filter(transect=="wetland/wetland-transition edge")

kruskal.test(phyloseq_class_desulfo_wetland$Abundance ~ phyloseq_class_desulfo_wetland$site)

##not significant

Kruskal-Wallis rank sum test

data:  phyloseq_class_desulfo_wetland$Abundance by phyloseq_class_desulfo_wetland$site
Kruskal-Wallis chi-squared = 2.0661, df = 2, p-value = 0.3559


##desulfo by transect
kruskal.test(phyloseq_class_desulfo$Abundance ~ phyloseq_class_desulfo$transect)

Kruskal-Wallis rank sum test

data:  phyloseq_class_desulfo$Abundance by phyloseq_class_desulfo$transect
Kruskal-Wallis chi-squared = 1.6206, df = 1, p-value = 0.203

##desulfo by site
kruskal.test(phyloseq_class_desulfo$Abundance ~ phyloseq_class_desulfo$site)

Kruskal-Wallis rank sum test

data:  phyloseq_class_desulfo$Abundance by phyloseq_class_desulfo$site
Kruskal-Wallis chi-squared = 3.0976, df = 2, p-value = 0.2125

##Gammaproteobacteria
phyloseq_class_gamma<- phyloseq_class %>% filter(Class==" c__Gammaproteobacteria")
View(phyloseq_class_gamma)

library(ggh4x)
phyloseq_class_gamma$transect<- factor(phyloseq_class_gamma$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)
plot_class<- ggplot(phyloseq_class_gamma, aes(x= site, y = Abundance, fill =Class)) + 
  facet_nested(.~transect+horizon, scales="free", independent = "x", labeller = label_wrap_gen(width=20)) + theme(strip.text.x = element_text(size = 14, face="bold"))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_boxplot() +
  scale_fill_manual(values=c("#39568CFF", "#55C667FF", "orange"))+
  theme(axis.title.x = element_text(size=17)) + theme(axis.text.x = element_text(size=15))+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=20, face="bold"),
                                     axis.title.y=element_text(size=20, face="bold")) + 
  ggtitle("Class Gammaproteobacteria")+ ylim(0,0.15)+
  theme(legend.title = element_blank())+ theme(legend.text = element_blank()) #guides(fill=guide_legend(ncol=1))
plot_class

ggsave(filename = "boxplot_gamma_WLE_filtered0.01.tiff", plot = plot_class,
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)



##gamma by transect
kruskal.test(phyloseq_class_gamma$Abundance ~ phyloseq_class_gamma$transect)

Kruskal-Wallis rank sum test

data:  phyloseq_class_gamma$Abundance by phyloseq_class_gamma$transect
Kruskal-Wallis chi-squared = 0.242, df = 2, p-value = 0.886

##gamma by site
kruskal.test(phyloseq_class_gamma$Abundance ~ phyloseq_class_gamma$site)

Kruskal-Wallis rank sum test

data:  phyloseq_class_gamma$Abundance by phyloseq_class_gamma$site
Kruskal-Wallis chi-squared = 2.9663, df = 2, p-value = 0.2269







##barplot ALL class plot

devtools::install_github("doehm/evoPalette")
library(evoPalette)
launch_evo_palette()

palette_new56 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
                  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
                  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "firebrick2",
                  "darkslategrey", "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4", "lightgreen", "magenta",
                  "lightgoldenrod1", "lightgoldenrod4", "lightpink3")

install.packages("ggh4x")
library(ggh4x)
plot_class<- ggplot(phyloseq_class, aes(x= tree_number, y = Abundance, fill =Class)) + 
  facet_grid2(transect~site+horizon, scales="free", independent = "x", labeller = label_wrap_gen(width=10)) + theme(strip.text.x = element_text(size = 18))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(palette_new56))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=17),
                                     axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))+ guides(fill=guide_legend(ncol=2))
plot_class



ggsave(filename = "barplot_class.tiff", plot = plot_class,
       width = 35,
       height = 25, units = c("cm"),
       dpi = 300)

####new class level plot for WLE 

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
#phyloseq_class_new<- phyloseq_class %>% mutate(transect1= ifelse(transect=="WC", "wetland", transect))

phyloseq_class$transect <- factor(phyloseq_class$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)
View(phyloseq_class)
plot_class<- ggplot(phyloseq_class, aes(x= tree_number, y = Abundance, fill =Class)) + 
  facet_nested(.~transect+horizon+site, scales="free", independent = "x", labeller = label_wrap_gen(width=10)) + theme(strip.text.x = element_text(size = 20, face="bold"))+ 
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

ggsave(filename = "barplot_class_WLE_new.jpg", plot = plot_class,
       width = 80,
       height = 45, units = c("cm"),
       dpi = 300)

ggsave(filename = "barplot_class_WLE_new_samecolorsasCB.jpg", plot = plot_class,
       width = 80,
       height = 45, units = c("cm"),
       dpi = 300)

##this gives a plot that does not add up to 1, did both: filtered to taxa which are greater than 1% and 2% of total
ggsave(filename = "barplot_class_WLE_new_filtered1percent.jpg", plot = plot_class,
       width = 80,
       height = 45, units = c("cm"),
       dpi = 300)

##barplot phylum level
phyloseq_phylum <- phyloseq_merged %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum) # Sort data frame alphabetically by phylum/family etc

phyloseq_phylum$Phylum[phyloseq_phylum$Abundance < 0.01] <- "< 1% abund."

View(phyloseq_phylum)
palette_new56 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#60ffaf","#b7ffdb","#825121","#117744",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#ae09ea","#521899","#1e0047","#195637",
                  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD",
                  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "firebrick2",
                  "darkslategrey", "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4", "lightgreen", "magenta",
                  "lightgoldenrod1", "lightpink3", "lavenderblush1", "lightgoldenrod4", "mediumorchid1","#fa7efc","#ffb7ef", "#fcfc00","#ffff9e")

install.packages("ggh4x")


plot_phylum<- ggplot(phyloseq_phylum, aes(x= tree_number, y = Abundance, fill =Phylum)) + 
  facet_grid2(transect~site+horizon, scales="free", independent = "x", labeller = label_wrap_gen(width=10)) + theme(strip.text.x = element_text(size = 18))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(palette_new56))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=17),
                                     axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))+ guides(fill=guide_legend(ncol=1))
plot_phylum



ggsave(filename = "barplot_phylum.tiff", plot = plot_phylum,
       width = 35,
       height = 25, units = c("cm"),
       dpi = 300)


##new phylum level for WLE
palette_new56 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#60ffaf","#b7ffdb","#825121","#117744",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#1e0047","lightpink3","#195637",
                  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#ae09ea","#521899",
                  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "firebrick2",
                  "darkslategrey", "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4", "lightgreen", "magenta",
                  "lightgoldenrod1", "lavenderblush1", "lightgoldenrod4", "mediumorchid1","#fa7efc","#ffb7ef", "#fcfc00","#ffff9e")

phyloseq_phylum$transect <- factor(phyloseq_phylum$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)
plot_phylum<- ggplot(phyloseq_phylum, aes(x= tree_number, y = Abundance, fill =Phylum)) + 
  facet_nested(.~transect+horizon+site, scales="free", independent = "x", labeller = label_wrap_gen(width=10)) + theme(strip.text.x = element_text(size = 20, face="bold"))+ 
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
plot_phylum


ggsave(filename = "barplot_phylum_WLE.tiff", plot = plot_phylum,
       width = 70,
       height = 25, units = c("cm"),
       dpi = 300)



######### Principal Coordinate analysis #########


## PCoaA analysis all samples 

#metadata<- read.csv("sample_metadata_dna_cdna.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones
#keep.samples <- as.vector(metadata$SampleID)
#keep.samples

#phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
#phyloseq_merged_final


##pcoa
phyloseq_pcoa <- ordinate(
  physeq = phyloseq_merged, ##considering using phyloseq_rarefy directly
  method = "PCoA", 
  distance = "bray"
)
##plot by site
pcoa_WLE<-plot_ordination(
  physeq = phyloseq_merged,
  ordination = phyloseq_pcoa,
  color = "transect",
  shape = "site",
  title = "PCoA WLE") + 
  scale_color_manual(values = c("#CC79A7", "#0072B2", "#ea7f17"))+
  geom_point(aes(color = transect), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))

pcoa_WLE


ggsave(filename="pcoa-site-shape.TIFF", plot=pcoa_WLE, width=8, height=6, units="in", dpi=300)

###new pcoa

sample_data(phyloseq_merged)$transect <- factor(sample_data(phyloseq_merged)$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)

##pcoa
phyloseq_pcoa <- ordinate(
  physeq = phyloseq_merged, ##considering using phyloseq_rarefy directly
  method = "PCoA", 
  distance = "bray"
)
##plot by site
pcoa_WLE<-plot_ordination(
  physeq = phyloseq_merged,
  ordination = phyloseq_pcoa,
  color = "site",
  shape = "transect",
  title = "PCoA WLE") + 
  scale_color_manual(values = c("#CC79A7", "#E69F00", "#0072B2", "black"))+
  geom_point(aes(color = site), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))

pcoa_WLE



ggsave(filename="pcoa-transect-shape-WLE.TIFF", plot=pcoa_WLE, width=8, height=6, units="in", dpi=300)






##basic stats

set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged, method = "bray") 
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged)) 
View(sample_df)
set.seed(1)
library(vegan)
adonis2(phyloseq_bray ~ transect*site, by="terms", data = sample_df)

#using merged phyloseq object which uses OTU_rarefied table instead of rarefy object
##gives exact same results as rarefy object below

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ transect * site, data = sample_df, by = "terms")
Df SumOfSqs      R2      F Pr(>F)    
transect       2   3.8931 0.16026 7.0940  0.001 ***
  site           2   2.0686 0.08516 3.7695  0.001 ***
  transect:site  4   2.6902 0.11074 2.4511  0.001 ***
  Residual      57  15.6402 0.64384                  
Total         65  24.2921 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1





##using rarefy instead of merged phyloseq object
adonis2(formula = phyloseq_bray ~ transect * site, data = sample_df, by = "terms")
Df SumOfSqs      R2      F Pr(>F)    
transect       2   3.8931 0.16026 7.0940  0.001 ***
  site           2   2.0686 0.08516 3.7695  0.001 ***
  transect:site  4   2.6902 0.11074 2.4511  0.001 ***
  Residual      57  15.6402 0.64384                  
Total         65  24.2921 1.00000                  
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#########################################################
################Alpha Diversity WLE #########################
#########################################################


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
  r <- rarefy_even_depth(phyloseq_merged, sample.size = 1000, verbose = FALSE, replace = TRUE) #no rarefaction a second time.
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed"))) ##changed r to phyloseq_merged, no second rarefaction
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson"))) ##changed r to phyloseq_merged, no second rarefaction
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

s <- read.csv("WLE_metadata.csv")
alphadiv <- merge(alpha, s, by = "SampleID") 
write.csv(alphadiv,file = "alphadiv_WLE.csv")


##subset if needed
alphadiv_wle<-alphadiv[alphadiv$region=="WLE",]
dim(alphadiv_wle)
alphadiv_cb<-alphadiv[alphadiv$region=="CB",]
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

### Inverse Simpson WLE
alphadiv_WLE <- read.csv("alphadiv_WLE.csv")
pd<-position_dodge(0.7)
alphadiv_WLE_InvS =  alphadiv_WLE %>% filter(measure=="Inverse Simpson")

View(alphadiv_WLE_InvS)


#### Richness WLE
alphadiv_WLE <- read.csv("alphadiv_WLE.csv")
pd<-position_dodge(0.7)
alphadiv_WLE_rich =  alphadiv_WLE %>% filter(measure=="Richness")

View(alphadiv_WLE_rich)



##combine richness and evenness together in panel


View(alphadiv_WLE_InvS)


View(alphadiv_WLE_rich)

alphadiv_combined= rbind(alphadiv_WLE_rich, alphadiv_WLE_InvS)
View(alphadiv_combined)



write.csv(alphadiv_combined, "alphadiv_combined.csv")
alphadiv_combined<- read.csv("alphadiv_combined.csv")

#alphadiv_new<- alphadiv %>% mutate(transect1=ifelse(transect=="WC", "wetland", transect))

alphadiv$transect <- factor(alphadiv$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)

plot_WLE<- ggplot(alphadiv, aes(x= transect,y=mean,color=site)) + theme_classic()+ scale_color_manual(values=c("#39568CFF", "#55C667FF", "orange")) +
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
  scale_x_discrete(name = 'transect', 
                     breaks = c('upland', 'transition', 'wetland/wetland-transition edge'), 
                     labels = c('upland', 'transition', 'wetland/wetland-\ntransition edge'))+
  labs(y = "measure")+ xlab("transect")+theme(axis.text.x = element_text(colour = "black", size = 17,angle = 90, vjust = 0.5, hjust=1))+
  theme( axis.text.y = element_text(size = 16, colour = "black"), 
         axis.title.y= element_text(size=16, color="black"),
         axis.title.x= element_text(size=16, color="black"),
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =14, colour = "black"),
         legend.text = element_text(size = 12, colour = "black")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")), color= guide_legend(order = 1)) 


plot_WLE


##Figure
ggsave(filename="alphadiv_rich_even_WLE.TIFF", plot=plot_WLE, width=25, height=15, units="cm", dpi=300)

##alpha diversity ANOVA stats
library(ggplot2)
library(grid)
library(lattice)
library(multcompView)
library(tidyverse)
library(dplyr)
alphadiv_rich <- read.csv("alphadiv_WLE.csv") ## combined richness and evenness estimates of all samples


pd<-position_dodge(0.7)


alphadiv_rich =  alphadiv_rich %>% filter(measure=="Inverse Simpson")
alphadiv_rich =  alphadiv_rich %>% filter(measure=="Richness")

##Conduct the following code for each of the datasets selected above
View(alphadiv_rich)
options(contrasts = c("contr.sum", "contr.poly")) ### must set this before running the model for type III tests
model = lm(mean~ transect+site + transect:site,
           data=alphadiv_rich)

shapiro.test(resid(model)) #W statistic > 0.9 means normality okay #check normality to flag outliers 


##inverse simpson
Shapiro-Wilk normality test

data:  resid(model)
W = 0.86948, p-value = 5.104e-06


##richness
Shapiro-Wilk normality test

data:  resid(model)
W = 0.91803, p-value = 0.0003315


ee= as.matrix(resid(model))
qqnorm (ee) #look for outliers (far away from curve) and test for normality 
ee

library(car)
#options(contrasts = c("contr.sum", "contr.poly"))
Anova(model, type=c("III"))

####Results inverse simpson

Anova Table (Type III tests)

Response: mean
Sum Sq Df  F value    Pr(>F)    
(Intercept)   213617  1 109.8522 6.533e-15 ***
  transect       19271  2   4.9550   0.01037 *  
  site            9253  2   2.3791   0.10177    
transect:site   4716  4   0.6063   0.65969    
Residuals     110841 57                       
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##Richness
Anova Table (Type III tests)

Response: mean
Sum Sq Df  F value    Pr(>F)    
(Intercept)   784464  1 122.0754 8.567e-16 ***
  transect       54333  2   4.2276   0.01941 *  
  site           34914  2   2.7166   0.07466 .  
transect:site  21689  4   0.8438   0.50328    
Residuals     366285 57                       
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


##add letters on graph 
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(multcompView)
library(multcomp)


#set up model
model <- lm(mean~transect, data = alphadiv_rich)

# get (adjusted) weight means per group
library(emmeans)
model_means <- emmeans(object = model,
                       specs = "transect")

cld <- cld(object = model_means,
           adjust = "Tukey",
           Letters = letters,
           alpha = 0.05)

View(cld)

combined<- full_join(alphadiv_rich,cld, by=c("transect"))
View(combined)
###

#plot<-combined%>%
#group_by("Soil")%>% unique(combined$.group)>1

library(dplyr)
combined<- as.tibble(combined)
class(combined)

data<- combined %>% mutate(group=.group)
View(data)
data$group <-str_remove(data$group, " ") ##run twice/thrice

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
summarized <- df_new %>% dplyr::select(c('transect', 'mean', 'group')) %>% group_by(transect)%>% summarise(max_mean=max(mean), across(group))%>% distinct(.keep_all=TRUE)
View(summarized)

summarized$transect <- factor(summarized$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)



#alphadiv_rich$site<-factor(alphadiv_rich$site,levels = c("CC", "PR", "OWC"), ordered=TRUE)      

#alphadiv_rich$transect<-factor(alphadiv_rich$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)         


df_new$transect<- factor(df_new$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)         
  
p<- ggplot(df_new, aes(transect, 
                      mean, 
                      colour = site)) +
  theme_classic()+ scale_color_manual(values=c("#39568CFF", "#55C667FF", "orange"))+
  geom_boxplot(outlier.shape = NA,
              fill = "white", 
              outlier.colour = NA, 
              position = position_dodge(width = 0.9)) + 
  geom_point(position = position_jitterdodge()) +
  scale_x_discrete(name = 'transect', 
                   breaks = c('upland', 'transition', 'wetland/wetland-transition edge'), 
                   labels = c('upland', 'transition', 'wetland/wetland-\ntransition edge'))+
  geom_text(data=summarized, aes(x=summarized$transect, y = 400, 
                                   group=summarized$transect, label = summarized$group, vjust=-0.3), 
            check_overlap = T, 
            position = position_dodge(width = 0.9), col = "black") + 
  ylab("Richness (Number of Observed ASVs)") +
  xlab("Transect") +
  theme(axis.text.x = element_text(colour = "black", size = 17,angle = 90, vjust = 0.5, hjust=1))+
  theme( axis.text.y = element_text(size = 16, colour = "black"), 
         axis.title.y= element_text(size=16, color="black"),
         axis.title.x= element_text(size=16, color="black"),
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =14, colour = "black"),
         legend.text = element_text(size = 12, colour = "black")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")), color= guide_legend(order = 1)) 

  
p

ggsave(filename = "alphadiv_letters_WLE_InvSimp.tiff", plot = p,
       width = 15,
       height = 15, units = c("cm"),
       dpi = 300)

ggsave(filename = "alphadiv_letters_WLE_Richness_ASV.tiff", plot = p,
       width = 15,
       height = 15, units = c("cm"),
       dpi = 300)


###############################
#########indicator taxa for WLE
###############################

###indicator for transect site in WLE


##load indicator species analysis package

install.packages("indicspecies")
library(indicspecies)

###########

##phyloseq object
library(phyloseq)
library(tidyverse)
otu = read.csv("OTU_rarefied1k.csv", sep=",", row.names=1) ##USING RAREFIED DATA here unlike unrarefied data used for site specific indicators  

otu_final = otu[rowSums(otu[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final)



#otu<- otu %>% select(starts_with("s"))
#otu<- otu %>% select(!c(s20230217096,s20230217120,s20230217121,s20230217122, s20230217123, s20230217124))
View(otu_final)
tax = read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
tax = as.matrix(tax)

metadata = read.csv("WLE_metadata.csv", sep=",", row.names=1) 


OTU = otu_table(otu_final, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 3257 taxa and 66 samples ]
sample_data() Sample Data:       [ 66 samples by 5 sample variables ]
tax_table()   Taxonomy Table:    [ 3257 taxa by 8 taxonomic ranks ]

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

#sample_data(phyloseq_merged)<- data.frame(sample_data(phyloseq_merged)) %>% mutate(transect1=ifelse(transect=="WC", "wetland", transect))

View(sample_data(phyloseq_merged))

#sample_data(phyloseq_merged)$transect1 <- factor(sample_data(phyloseq_merged)$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)


ind_transectWLE <-
  GetIndicators(dataframe=phyloseq_merged, "transect")


##checking

View(ind_transectWLE) ##WLE shows no indicators across transects

ind_transectWLE_subset<- ind_transectWLE %>% subset(s.upland==0 & s.transition==0 & s.wetland==1)
View(ind_transectWLE_subset) 

write.csv(ind_transectWLE_subset, "ind_WLE_wetland_unique.csv")



library(tidyverse)

write.csv(ind_transectWLE, "indicator_transect_WLE.csv")

##creating a dataframe for CB only indicators and passing on as a join to the WLE count table. 

##need to use unrarefied count tables for WLE here for the join

ind_transect_WLE<- read.csv("indicator_transect_WLE.csv", row.names=1)
View(ind_transect_WLE)
ind_transect_WLE_wetland<- ind_transect_WLE %>% subset(s.wetland==1)
View(ind_transect_WLE_wetland) ##2803 ASVs

merge_WLEind_wetland<- merge(ind_transect_WLE_wetland, otu_final, by=0)
View(merge_WLEind_wetland)

write.csv(merge_WLEind_wetland, "WLE_wetlandindAll_count.csv")

##remove extra columns 

merge_WLEind_wetland <- as.data.frame(merge_WLEind_wetland)
rnames <-merge_WLEind_wetland[,1] 
rownames(merge_WLEind_wetland) <- rnames  # assign row names
View(merge_WLEind_wetland)
merge_WLEind_wetland_new <- merge_WLEind_wetland %>% select(!c(Row.names, OTU, s.upland, s.transition, s.wetland,index, stat, sequencefeature, p.value, Domain, Phylum, Class, Order, Family, Genus, Species)) ##add other columns
View(merge_WLEind_wetland_new)

write.csv(merge_WLEind_wetland_new, "indicatorWLE_wetland_count_final.csv")

##plot the most abundant ones only 

merge_WLEind_wetland_new<- read.csv("indicatorWLE_wetland_count_final.csv", row.names=1)
View(merge_WLEind_wetland_new)
merge_WLEind_wetland_new<-merge_WLEind_wetland_new[order(rowSums(merge_WLEind_wetland_new), decreasing = TRUE),]
View(merge_WLEind_wetland_new)
merge_WLEind_wetland_new_20 <- merge_WLEind_wetland_new[c(1:20),] ##select top 20 most abundant Classes
View(merge_WLEind_wetland_new_20)

write.csv(merge_WLEind_wetland_new_20, "WLEind_wetland_top20.csv")

##max_standardize 
max_standardize_WLEwetland <-decostand(merge_WLEind_wetland_new_20, method = "max", MARGIN = 1) ##from Jackson's phil trans paper
View(max_standardize_WLEwetland)

tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("metadata_WLE.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(max_standardize_WLEwetland, taxa_are_rows = TRUE)
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
#phyloseq_otu_new<- phyloseq_otu %>% mutate(transect1= ifelse(transect=="WC", "wetland", transect))
View(phyloseq_otu)

phyloseq_otu$transect <- factor(phyloseq_otu$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)
plot<- ggplot(phyloseq_otu, aes(x= transect, y = Abundance, fill =transect)) + 
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



ggsave(filename = "ASV_WLEind_wetland_top20.jpg", plot = plot,
       width =48,
       height = 10, units = c("in"),
       dpi = 300)



##Since no indicators for transects are seen in WLE, next checking if CB wetland indicators can be detected in WLE

##creating a dataframe for CB only indicators and passing on as a join to the WLE count table. 

##need to use unrarefied count tables for WLE here for the join

ind_transect<- read.csv("indicator_transect_CB.csv", row.names=1)
View(ind_transect)
ind_transect_wetland<- ind_transect %>% subset(s.wetland==1)
View(ind_transect_wetland)


otu_all<- read.csv("COMPASS_CB_WLE_count.csv", row.names=1)
View(otu_all)
otu_WLE<- otu_all %>% dplyr::select(!starts_with("s"))
View(otu_WLE)
otu_WLE<- otu_WLE %>% dplyr::select(!starts_with("Zymo"))

otu_WLE_final = otu_WLE[rowSums(otu_WLE[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_WLE_final)## slightly more OTUs here than in rarefied WLE count table

merge_cbind_wle<- merge(ind_transect_wetland, otu_WLE_final, by=0)
View(merge_cbind_wle)

write.csv(merge_cbind_wle, "CBind_wetland_inWLEunrarefiedCount.csv")

##remove extra columns 

merge_cbind_wle <- as.data.frame(merge_cbind_wle)
rnames <-merge_cbind_wle[,1] 
rownames(merge_cbind_wle) <- rnames  # assign row names
View(merge_cbind_wle)
merge_cbind_wle_new <- merge_cbind_wle %>% dplyr::select(!c(Row.names, OTU, s.upland, s.transition, s.wetland, index, stat, sequencefeature, p.value, Domain, Phylum, Class, Order, Family, Genus, Species)) ##add other columns
View(merge_cbind_wle_new) ##101 OTUs

write.csv(merge_cbind_wle_new, "indicatorCB_wetland_inWLE_count_final.csv")

##plot the most abundant ones only 

merge_cbind_wle_new<- read.csv("indicatorCB_wetland_inWLE_count_final.csv", row.names=1)
merge_cbind_wle_new<-merge_cbind_wle_new[order(rowSums(merge_cbind_wle_new), decreasing = TRUE),]
View(merge_cbind_wle_new)
merge_cbind_wle_new_20 <- merge_cbind_wle_new[c(1:20),] ##select top 20 most abundant Classes
View(merge_cbind_wle_new_20)

write.csv(merge_cbind_wle_new_20, "CBind_wetland_inWLE_top20.csv")

##max_standardize 
max_standardize_WLE <-decostand(merge_cbind_wle_new_20, method = "max", MARGIN = 1) ##from Jackson's phil trans paper
View(max_standardize_WLE)

tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("WLE_metadata.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(max_standardize_WLE, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 20 taxa and 85 samples ]
sample_data() Sample Data:       [ 85 samples by 5 sample variables ]
tax_table()   Taxonomy Table:    [ 20 taxa by 8 taxonomic ranks ]

phyloseq_otu <- phyloseq_merged %>%
  #tax_glom(taxrank = "Class")   %>% 
  #transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()
View(phyloseq_otu)

##combining class and otu information

phyloseq_otu$ASVclass <- paste(phyloseq_otu$Class, phyloseq_otu$OTU)
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
#phyloseq_otu_new<- phyloseq_otu %>% mutate(transect1= ifelse(transect=="WC", "wetland", transect))
#View(phyloseq_otu_new)

phyloseq_otu$transect <- factor(phyloseq_otu$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)
plot<- ggplot(phyloseq_otu, aes(x= transect, y = Abundance, fill =transect)) + 
  facet_grid2(site~ASVclass, scales="free", independent = "x", labeller = label_wrap_gen(width=25)) + theme(strip.text.x = element_text(size = 18, face="bold"))+ 
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



ggsave(filename = "ASV_CBind_wetland_inWLE1.jpg", plot = plot,
       width =49,
       height = 10, units = c("in"),
       dpi = 300)



######################
###CB data checking taxa in all samples at class level and those in samples 96, 120, 121, and 122. This
##is before knowing what the sample annotations were. Still figuring them out at this stage.
######################

library(tidyverse)
library(reshape2)
data<- read.csv("COMPASS_CB_WLE_count.csv", sep=",",row.names=1)


View(data)
data_CB<- data %>% select(starts_with("s"))
View(data_CB)

any(rowSums(data_CB[])<1)


data_CB = data_CB[rowSums(data_CB[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(data_CB)


library(phyloseq)


tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("COMPASS_metadata_CB.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(data_CB, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged


phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 47110 taxa and 124 samples ]
sample_data() Sample Data:       [ 124 samples by 8 sample variables ]
tax_table()   Taxonomy Table:    [ 47110 taxa by 8 taxonomic ranks ]

##selecting 4 samples before rel abund calculations, suspecting that relabund calculations are failing because of memory time out
data_CB<- data_CB %>% select(c(s20230217096,s20230217120,s20230217121,s20230217122))
View(data_CB)
data_CB = data_CB[rowSums(data_CB[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(data_CB)


tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("COMPASS_metadata_CB.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(data_CB, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged


phyloseq_class <- phyloseq_merged %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."

View(phyloseq_class)

phyloseq_genus <- phyloseq_merged %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Genus) # Sort data frame alphabetically by phylum/family etc

write.csv(phyloseq_genus, "phyloseq_genus_CB.csv")

phyloseq_genus<- read.csv("phyloseq_genus_CB.csv")
View(phyloseq_genus)
phyloseq_genus$Genus[phyloseq_genus$Abundance < 0.005] <- "< 0.5% abund."

View(phyloseq_genus)



devtools::install_github("doehm/evoPalette")
library(evoPalette)
launch_evo_palette()

palette_new56 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
                  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
                  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "firebrick2",
                  "darkslategrey", "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4", "lightgreen", "magenta",
                  "lightgoldenrod1", "lightgoldenrod4", "lightpink3")

install.packages("ggh4x")
library(ggh4x)

##genus level

phyloseq_genus_4<- phyloseq_genus %>% filter(Sample== c("s20230217096","s20230217120","s20230217121","s20230217122"))
View(phyloseq_genus_4)

sum(phyloseq_genus$Abundance)



plot_genus<- ggplot(phyloseq_genus, aes(x= Sample, y = Abundance, fill =Genus)) + 
  #facet_grid2(~PCR.Plate, scales="free", independent = "x", labeller = label_wrap_gen(width=10)) 
  #theme(strip.text.x = element_text(size = 18))+ 
  #theme(strip.text.y = element_text(size = 18))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(palette_new56))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + #theme(axis.text.x = element_blank())+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=20))+
  #
  #guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=17),
                                     axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))+ guides(fill=guide_legend(ncol=1))
plot_genus



ggsave(filename = "barplot_genus_CB_4samples_0.5.tiff", plot = plot_genus,
       width = 35,
       height = 35, units = c("cm"),
       dpi = 300)

##class level
plot_class<- ggplot(phyloseq_class, aes(x= Sample, y = Abundance, fill =Class)) + 
  facet_grid2(~PCR.Plate, scales="free", independent = "x", labeller = label_wrap_gen(width=10)) + theme(strip.text.x = element_text(size = 18))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(palette_new56))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + #theme(axis.text.x = element_blank())+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=17),
                                     axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))+ guides(fill=guide_legend(ncol=2))
plot_class



ggsave(filename = "barplot_class_CB.tiff", plot = plot_class,
       width = 80,
       height = 30, units = c("cm"),
       dpi = 300)


##new metadata after accounting for samples in google spreadsheet
##after figuring out sample metadata

metadata1<- read.csv("CB_metadata.csv")
metadata2<- read.csv("COMPASS_metadata_CB.csv")

View(metadata1)

View(metadata2)

merge_metadata<- full_join(metadata1, metadata2, by="SampleID")

View(merge_metadata)

write.csv(merge_metadata, "combined_metadataCB.csv")

data<- read.csv("COMPASS_CB_WLE_count.csv", sep=",",row.names=1)
View(data)
data_CB<- data %>% select(starts_with("s"))
View(data_CB)
data_CB<- data_CB %>% select(!c(s20230217096,s20230217120,s20230217121,s20230217122, s20230217123, s20230217124))

library(phyloseq)


tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("combined_metadataCB.csv", sep=",", row.names=1) ##manually changed this file to remove the samples which were missing from google sheet and were also not sequenced,
#namely 129, 130, 137, 167, 205 and make the first column SampleID starting with "s---"

tax = as.matrix(tax) 

otu<- otu_table(data_CB, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 50421 taxa and 118 samples ]
sample_data() Sample Data:       [ 118 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 50421 taxa by 8 taxonomic ranks ]


sample_sum_df <- data.frame(sum = sample_sums(phyloseq_merged))
sum(sample_sum_df) ##426275 reads in WLE, 8845545 reads in CB -- lots more in CB, will lose much information during rarefaction for this dataset when doing combined analysis

# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

smin <- min(sample_sums(phyloseq_merged))
smin # 11292
smean <- mean(sample_sums(phyloseq_merged))
smean # 74962.25
smax <- max(sample_sums(phyloseq_merged))
smax # 136223


taxa_are_rows(phyloseq_merged)
#### Rarefaction curves

library(scales)
library(vegan)
library(dplyr)

mat<- as(t(otu_table(phyloseq_merged)), "matrix")
system.time(rarecurve(mat, step=1000, label = FALSE, col="cyan4")) ##looks good, very even reads across samples

##rarefaction

sample_df <- data.frame(sample_data(phyloseq_merged))
View(sample_df)
write.csv(sample_df, "samples_pre_rarefied_CB.csv")

library(vegan)
#rarecurve(t(otu_table(phyloseq_merged)), step=50)
phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged, sample.size = 40000, rngseed = TRUE, trimOTUs=FALSE)

#`set.seed(TRUE)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(TRUE); .Random.seed` for the full vector
...
#2 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 
  
  #s20230217052s20230217078




sample_df <- data.frame(sample_data(phyloseq_rarefy))

View(sample_df)
write.csv(sample_df, "samples_rarefied_CB.csv")

##use anti-join to select non-present samples in rarefied object
#This join is like eg: for df1-df2, it selects all rows from df1 that are not present in df2.
df1<- read.csv("samples_pre_rarefied_CB.csv")
df2 <- read.csv("samples_rarefied_CB.csv")
df= df1 %>% anti_join(df2,by="X")

View(df)

write.csv(df, "samples-lost-rarefaction-CB.csv")

OTU2 = as(otu_table(phyloseq_rarefy), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, "OTU_rarefied_CB.csv")

otu_final<- read.csv("OTU_rarefied_CB.csv", row.names=1) ##change name 1k--- depending on rarefaction curve
View(otu_final)
any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #5590 otus with no counts
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) #44831 otus 


##continuing with analysis using prior rarefied data (repetition of the above steps)

otu_final<- read.csv("OTU_rarefied1k.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #47164 otus with no counts
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final)

library(phyloseq)


tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("combined_metadataCB.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged


phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 44831 taxa and 116 samples ]
sample_data() Sample Data:       [ 116 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 44831 taxa by 8 taxonomic ranks ]


#metadata<- read.csv("sample_metadata_dna_cdna.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones
#keep.samples <- as.vector(metadata$SampleID)
#keep.samples

#phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
#phyloseq_merged_final

##barplot class level
phyloseq_class <- phyloseq_merged %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."

View(phyloseq_class)



devtools::install_github("doehm/evoPalette")
library(evoPalette)
launch_evo_palette()

palette_new56 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
                  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
                  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "firebrick2",
                  "lightpink3", "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4", "lightgreen", "magenta",
                  "lightgoldenrod1", "lightgoldenrod4","darkslategrey" )

install.packages("ggh4x")
library(ggh4x)
phyloseq_class_new<- phyloseq_class %>% mutate(transect1= ifelse(transect=="WC", "wetland", transect))

phyloseq_class_new$transect1 <- factor(phyloseq_class_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)
View(phyloseq_class_new)
plot_class<- ggplot(phyloseq_class_new, aes(x= tree_number, y = Abundance, fill =Class)) + 
  facet_nested(.~transect1+horizon+site, scales="free", independent = "x", labeller = label_wrap_gen(width=10)) + theme(strip.text.x = element_text(size = 20, face="bold"))+ 
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



ggsave(filename = "barplot_class_CB_new6.jpg", plot = plot_class,
       width = 80,
       height = 36, units = c("cm"),
       dpi = 300)




##Desulfobacteria class plot only
phyloseq_class_new<- phyloseq_class %>% mutate(transect1= ifelse(transect=="WC", "wetland", transect))
phyloseq_class_new$transect1 <- factor(phyloseq_class_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)

phyloseq_class_desulfo_cb_all<- phyloseq_class_new %>% filter(Class==" c__Desulfobacteria")
phyloseq_class_desulfo_cb_filtered0.01<- phyloseq_class_new %>% filter(Class==" c__Desulfobacteria")
View(phyloseq_class_desulfo_cb_all)
View(phyloseq_class_desulfo_cb_filtered0.01)
library(ggh4x)
plot_class<- ggplot(phyloseq_class_desulfo_cb_filtered0.01, aes(x= site, y = Abundance, fill =Class)) + 
  facet_nested(.~transect1+horizon, scales="free", independent = "x", labeller = label_wrap_gen(width=20)) + theme(strip.text.x = element_text(size = 14, face="bold"))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_boxplot() +
  scale_fill_manual(values=c("#39568CFF", "#55C667FF", "orange"))+
  theme(axis.title.x = element_text(size=17)) + theme(axis.text.x = element_text(size=15, angle=90))+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=20, face="bold"),
                                     axis.title.y=element_text(size=20, face="bold")) + 
  ggtitle("Class Desulfobacteria")+
  theme(legend.title = element_blank())+ theme(legend.text = element_blank()) #guides(fill=guide_legend(ncol=1))
plot_class

ggsave(filename = "boxplot_desulfo_CB_filtered0.01.tiff", plot = plot_class,
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)

##kruskal test for wetland desulfo between sites

phyloseq_class_desulfo_wetland_ohorizon<- phyloseq_class_desulfo_cb_filtered0.01 %>% filter(transect1=="wetland"& horizon=="O")

kruskal.test(phyloseq_class_desulfo_wetland_ohorizon$Abundance ~ phyloseq_class_desulfo_wetland_ohorizon$site)

##not significant
Kruskal-Wallis rank sum test

data:  phyloseq_class_desulfo_wetland_ohorizon$Abundance by phyloseq_class_desulfo_wetland_ohorizon$site
Kruskal-Wallis chi-squared = 0.038961, df = 1, p-value = 0.8435

##by transect
kruskal.test(phyloseq_class_desulfo_cb_filtered0.01$Abundance ~ phyloseq_class_desulfo_cb_filtered0.01$transect1)


Kruskal-Wallis rank sum test

data:  phyloseq_class_desulfo_cb_filtered0.01$Abundance by phyloseq_class_desulfo_cb_filtered0.01$transect1
Kruskal-Wallis chi-squared = 1.7696, df = 1, p-value = 0.1834


##by site 

kruskal.test(phyloseq_class_desulfo_cb_filtered0.01$Abundance ~ phyloseq_class_desulfo_cb_filtered0.01$site)

Kruskal-Wallis rank sum test

data:  phyloseq_class_desulfo_cb_filtered0.01$Abundance by phyloseq_class_desulfo_cb_filtered0.01$site
Kruskal-Wallis chi-squared = 7.5285, df = 2, p-value = 0.02319

##by horizon
kruskal.test(phyloseq_class_desulfo_cb_filtered0.01$Abundance ~ phyloseq_class_desulfo_cb_filtered0.01$horizon)
Kruskal-Wallis rank sum test

data:  phyloseq_class_desulfo_cb_filtered0.01$Abundance by phyloseq_class_desulfo_cb_filtered0.01$horizon
Kruskal-Wallis chi-squared = 7.5542, df = 2, p-value = 0.02289


##Acidobacteriae


phyloseq_class_acido_cb_filtered0.01<- phyloseq_class_new %>% filter(Class==" c__Acidobacteriae")

View(phyloseq_class_acido_cb_filtered0.01)
library(ggh4x)
plot_class<- ggplot(phyloseq_class_acido_cb_filtered0.01, aes(x= site, y = Abundance, fill =Class)) + 
  facet_nested(.~transect1+horizon, scales="free", independent = "x", labeller = label_wrap_gen(width=20)) + theme(strip.text.x = element_text(size = 14, face="bold"))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_boxplot() +
  scale_fill_manual(values=c("#39568CFF", "#55C667FF", "orange"))+
  theme(axis.title.x = element_text(size=17)) + theme(axis.text.x = element_text(size=15, angle=90))+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=20, face="bold"),
                                     axis.title.y=element_text(size=20, face="bold")) + 
  ggtitle("Class Acidobacteriae")+
  theme(legend.title = element_blank())+ theme(legend.text = element_blank()) #guides(fill=guide_legend(ncol=1))
plot_class

ggsave(filename = "boxplot_acido_CB_filtered0.01.tiff", plot = plot_class,
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)

##kruskal test for wetland acido between sites

#phyloseq_class_acido_wetland_ohorizon<- phyloseq_class_acido_cb_filtered0.01 %>% filter(transect1=="wetland"& horizon=="O")

#kruskal.test(phyloseq_class_acido_wetland_ohorizon$Abundance ~ phyloseq_class_acido_wetland_ohorizon$site)

##results


##by transect
kruskal.test(phyloseq_class_acido_cb_filtered0.01$Abundance ~ phyloseq_class_acido_cb_filtered0.01$transect1)

Kruskal-Wallis rank sum test

data:  phyloseq_class_acido_cb_filtered0.01$Abundance by phyloseq_class_acido_cb_filtered0.01$transect1
Kruskal-Wallis chi-squared = 4.7377, df = 2, p-value = 0.09359


##by site 

kruskal.test(phyloseq_class_acido_cb_filtered0.01$Abundance ~ phyloseq_class_acido_cb_filtered0.01$site)


Kruskal-Wallis rank sum test

data:  phyloseq_class_acido_cb_filtered0.01$Abundance by phyloseq_class_acido_cb_filtered0.01$site
Kruskal-Wallis chi-squared = 3.0042, df = 2, p-value = 0.2227

##by horizon
kruskal.test(phyloseq_class_acido_cb_filtered0.01$Abundance ~ phyloseq_class_acido_cb_filtered0.01$horizon)

Kruskal-Wallis rank sum test

data:  phyloseq_class_acido_cb_filtered0.01$Abundance by phyloseq_class_acido_cb_filtered0.01$horizon
Kruskal-Wallis chi-squared = 0.81538, df = 2, p-value = 0.6652



##Actinobacteria

phyloseq_class_actino_cb_filtered0.01<- phyloseq_class_new %>% filter(Class==" c__Actinobacteria")

View(phyloseq_class_actino_cb_filtered0.01)
library(ggh4x)
plot_class<- ggplot(phyloseq_class_actino_cb_filtered0.01, aes(x= site, y = Abundance, fill =Class)) + 
  facet_nested(.~transect1+horizon, scales="free", independent = "x", labeller = label_wrap_gen(width=20)) + theme(strip.text.x = element_text(size = 14, face="bold"))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_boxplot() +
  scale_fill_manual(values=c("#39568CFF", "#55C667FF", "orange"))+
  theme(axis.title.x = element_text(size=17)) + theme(axis.text.x = element_text(size=15, angle=90))+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=20, face="bold"),
                                     axis.title.y=element_text(size=20, face="bold")) + 
  ggtitle("Class Actinobacteria")+
  theme(legend.title = element_blank())+ theme(legend.text = element_blank()) #guides(fill=guide_legend(ncol=1))
plot_class

ggsave(filename = "boxplot_actino_CB_filtered0.01.tiff", plot = plot_class,
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)

##kruskal test for wetland acido between sites

#phyloseq_class_actino_wetland_ohorizon<- phyloseq_class_actino_cb_filtered0.01 %>% filter(transect1=="wetland"& horizon=="O")

#kruskal.test(phyloseq_class_actino_wetland_ohorizon$Abundance ~ phyloseq_class_actino_wetland_ohorizon$site)

##results


##by transect
kruskal.test(phyloseq_class_actino_cb_filtered0.01$Abundance ~ phyloseq_class_actino_cb_filtered0.01$transect1)

Kruskal-Wallis rank sum test

data:  phyloseq_class_actino_cb_filtered0.01$Abundance by phyloseq_class_actino_cb_filtered0.01$transect1
Kruskal-Wallis chi-squared = 0.91918, df = 1, p-value = 0.3377



##by site 

kruskal.test(phyloseq_class_actino_cb_filtered0.01$Abundance ~ phyloseq_class_actino_cb_filtered0.01$site)

Kruskal-Wallis rank sum test

data:  phyloseq_class_actino_cb_filtered0.01$Abundance by phyloseq_class_actino_cb_filtered0.01$site
Kruskal-Wallis chi-squared = 0.13869, df = 2, p-value = 0.933

##by horizon
kruskal.test(phyloseq_class_actino_cb_filtered0.01$Abundance ~ phyloseq_class_actino_cb_filtered0.01$horizon)


Kruskal-Wallis rank sum test

data:  phyloseq_class_actino_cb_filtered0.01$Abundance by phyloseq_class_actino_cb_filtered0.01$horizon
Kruskal-Wallis chi-squared = 0.29536, df = 2, p-value = 0.8627





##Verrucomicrobiae

phyloseq_class_verru_cb_filtered0.01<- phyloseq_class_new %>% filter(Class==" c__Verrucomicrobiae")

View(phyloseq_class_verru_cb_filtered0.01)
library(ggh4x)
plot_class<- ggplot(phyloseq_class_verru_cb_filtered0.01, aes(x= site, y = Abundance, fill =Class)) + 
  facet_nested(.~transect1+horizon, scales="free", independent = "x", labeller = label_wrap_gen(width=20)) + theme(strip.text.x = element_text(size = 14, face="bold"))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_boxplot() +
  scale_fill_manual(values=c("#39568CFF", "#55C667FF", "orange"))+
  theme(axis.title.x = element_text(size=17)) + theme(axis.text.x = element_text(size=15, angle=90))+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=20, face="bold"),
                                     axis.title.y=element_text(size=20, face="bold")) + 
  ggtitle("Class Verrucomicrobiae")+
  theme(legend.title = element_blank())+ theme(legend.text = element_blank()) #guides(fill=guide_legend(ncol=1))
plot_class

ggsave(filename = "boxplot_verru_CB_filtered0.01.tiff", plot = plot_class,
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)

##kruskal test for wetland acido between sites

#phyloseq_class_verru_wetland_ohorizon<- phyloseq_class_verru_cb_filtered0.01 %>% filter(transect1=="wetland"& horizon=="O")

#kruskal.test(phyloseq_class_verru_wetland_ohorizon$Abundance ~ phyloseq_class_verru_wetland_ohorizon$site)

##results


##by transect
kruskal.test(phyloseq_class_verru_cb_filtered0.01$Abundance ~ phyloseq_class_verru_cb_filtered0.01$transect1)

Kruskal-Wallis rank sum test

data:  phyloseq_class_verru_cb_filtered0.01$Abundance by phyloseq_class_verru_cb_filtered0.01$transect1
Kruskal-Wallis chi-squared = 7.3486, df = 2, p-value = 0.02537



##by site 

kruskal.test(phyloseq_class_verru_cb_filtered0.01$Abundance ~ phyloseq_class_verru_cb_filtered0.01$site)

Kruskal-Wallis rank sum test

data:  phyloseq_class_verru_cb_filtered0.01$Abundance by phyloseq_class_verru_cb_filtered0.01$site
Kruskal-Wallis chi-squared = 0.29195, df = 2, p-value = 0.8642


##by horizon
kruskal.test(phyloseq_class_verru_cb_filtered0.01$Abundance ~ phyloseq_class_verru_cb_filtered0.01$horizon)


Kruskal-Wallis rank sum test

data:  phyloseq_class_verru_cb_filtered0.01$Abundance by phyloseq_class_verru_cb_filtered0.01$horizon
Kruskal-Wallis chi-squared = 0.28917, df = 2, p-value = 0.8654


##barplot phylum level
phyloseq_phylum <- phyloseq_merged %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum) # Sort data frame alphabetically by phylum/family etc

phyloseq_phylum$Phylum[phyloseq_phylum$Abundance < 0.01] <- "< 1% abund."

View(phyloseq_phylum)
palette_new56 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#60ffaf","#b7ffdb","#825121","#117744",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#1e0047","lightpink3","#195637",
                  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#ae09ea","#521899",
                  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "firebrick2",
                  "darkslategrey", "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4", "lightgreen", "magenta",
                  "lightgoldenrod1", "lavenderblush1", "lightgoldenrod4", "mediumorchid1","#fa7efc","#ffb7ef", "#fcfc00","#ffff9e")

install.packages("ggh4x")

phyloseq_phylum_new<- phyloseq_phylum %>% mutate(transect1= ifelse(transect=="WC", "wetland", transect))
phyloseq_phylum_new$transect1 <- factor(phyloseq_phylum_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)
plot_phylum<- ggplot(phyloseq_phylum_new, aes(x= tree_number, y = Abundance, fill =Phylum)) + 
  facet_grid2(~transect1+site+horizon, scales="free", independent = "x", labeller = label_wrap_gen(width=10)) + theme(strip.text.x = element_text(size = 18))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(palette_new56))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=17),
                                     axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))+ guides(fill=guide_legend(ncol=1))
plot_phylum

##newplot

plot_phylum<- ggplot(phyloseq_phylum_new, aes(x= tree_number, y = Abundance, fill =Phylum)) + 
  facet_nested(.~transect1+horizon+site, scales="free", independent = "x", labeller = label_wrap_gen(width=10)) + theme(strip.text.x = element_text(size = 20, face="bold"))+ 
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
plot_phylum


ggsave(filename = "barplot_phylum_CB1.tiff", plot = plot_phylum,
       width = 70,
       height = 25, units = c("cm"),
       dpi = 300)

######### Principal Coordinate analysis #########


## PCoaA analysis all samples 

#metadata<- read.csv("sample_metadata_dna_cdna.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones
#keep.samples <- as.vector(metadata$SampleID)
#keep.samples

#phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
#phyloseq_merged_final


sample_data(phyloseq_merged)<- data.frame(sample_data(phyloseq_merged)) %>% mutate(transect1=ifelse(transect=="WC", "wetland", transect))

View(sample_data(phyloseq_merged))

sample_data(phyloseq_merged)$transect1 <- factor(sample_data(phyloseq_merged)$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)

##pcoa
phyloseq_pcoa <- ordinate(
  physeq = phyloseq_merged, 
  method = "PCoA", 
  distance = "bray"
)
##plot by site
pcoa_CB<-plot_ordination(
  physeq = phyloseq_merged,
  ordination = phyloseq_pcoa,
  color = "site",
  shape = "transect1",
  title = "PCoA CB 1") + 
  scale_color_manual(values = c("#CC79A7", "#E69F00", "#0072B2", "black"))+
  geom_point(aes(color = site), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))

pcoa_CB


### For main Figure 4 A
ggsave(filename="pcoa-transect-shape-CB1.TIFF", plot=pcoa_CB, width=8, height=6, units="in", dpi=300)

##basic stats

set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged, method = "bray") ##using merged object is fine here, gives same results as rarefy object
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged)) 
View(sample_df)
set.seed(1)
adonis2(phyloseq_bray ~ transect1*site, by="terms", data = sample_df)

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ transect1 * site, data = sample_df, by = "terms")
Df SumOfSqs      R2       F Pr(>F)    
transect1        2    7.067 0.14349 13.1231  0.001 ***
  site             2    5.958 0.12097 11.0636  0.001 ***
  transect1:site   4    7.416 0.15058  6.8857  0.001 ***
  Residual       107   28.812 0.58497                   
Total          115   49.254 1.00000                   
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


###combined data

##combined wle and cb metadata

wle_map<- read.csv("WLE_metadata.csv")

cb_map<- read.csv("combined_metadataCB.csv")

wle_cb_map<- full_join(wle_map, cb_map)

View(wle_cb_map)

map_edit<- wle_cb_map %>% mutate(region_new= ifelse(is.na(wle_cb_map$region), "CB", region)) %>% mutate(Sample=ifelse(is.na(wle_cb_map$Sample.ID), SampleID, Sample.ID))
View(map_edit)

class(wle_cb_map)
as.tibble(wle_cb_map)

wle_cb_map2<- map_edit[,-19]
View(wle_cb_map2)
rownames(wle_cb_map2) <- map_edit[,19]

View(wle_cb_map2)


data<- read.csv("COMPASS_CB_WLE_count.csv", sep=",",row.names=1)
data<- data %>% select(!c(s20230217096,s20230217120,s20230217121,s20230217122, s20230217123, s20230217124))
tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)


tax = as.matrix(tax) 

otu<- otu_table(data, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(wle_cb_map2)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged


##WLE had 121 samples in the metadata but only 85 were sequenced so the extras samples in the metadata file here are all from WLE. 
##121-86 from wle metadata csv file==35+203 == 244-6 samples==238 samples
##118 samples from CB and 85 from WLE == 203 samples
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 50421 taxa and 203 samples ]
sample_data() Sample Data:       [ 203 samples by 18 sample variables ]
tax_table()   Taxonomy Table:    [ 50421 taxa by 8 taxonomic ranks ]


sample_sum_df <- data.frame(sum = sample_sums(phyloseq_merged))
sum(sample_sum_df) ##426275 reads in WLE, 8845545 reads in CB -- lots more in CB, will lose much information during rarefaction for this dataset when doing combined analysis

##TOTAL READS - CB + WLE ==9271820 (matches the sum up of the separate data)


# Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

smin <- min(sample_sums(phyloseq_merged))
smin # 79
smean <- mean(sample_sums(phyloseq_merged))
smean # 45673.99
smax <- max(sample_sums(phyloseq_merged))
smax # 136223


taxa_are_rows(phyloseq_merged)
#### Rarefaction curves

library(scales)
library(vegan)
library(dplyr)

mat<- as(t(otu_table(phyloseq_merged)), "matrix")
system.time(rarecurve(mat, step=1000, label = FALSE, col="cyan4")) ##COMBINED RAREFACTION IS VERY UNEVEN, RAREFYING WILL OMIT A LOT OF READS DUE TO LOW LIBRARY SIZE IN WLE

##rarefaction

sample_df <- data.frame(sample_data(phyloseq_merged))
View(sample_df)
write.csv(sample_df, "samples_pre_rarefied_CB_WLE.csv")

library(vegan)
#rarecurve(t(otu_table(phyloseq_merged)), step=50)
phyloseq_rarefy<- rarefy_even_depth(phyloseq_merged, sample.size = 1000, rngseed = TRUE, trimOTUs=FALSE)

#`set.seed(TRUE)` was used to initialize repeatable random subsampling.
#Please record this for your records so others can reproduce.
#Try `set.seed(TRUE); .Random.seed` for the full vector
...
#2 samples removedbecause they contained fewer reads than `sample.size`.
#Up to first five removed samples are: 

#s20230217052s20230217078




sample_df <- data.frame(sample_data(phyloseq_rarefy))

View(sample_df)
write.csv(sample_df, "samples_rarefied_CB_WLE.csv")

##use anti-join to select non-present samples in rarefied object
#This join is like eg: for df1-df2, it selects all rows from df1 that are not present in df2.
df1<- read.csv("samples_pre_rarefied_CB_WLE.csv")
df2 <- read.csv("samples_rarefied_CB_WLE.csv")
df= df1 %>% anti_join(df2,by="X")

View(df)

write.csv(df, "samples-lost-rarefaction-CB-WLE.csv")

OTU2 = as(otu_table(phyloseq_rarefy), "matrix")

# Coerce to data.frame
OTUdf = as.data.frame(OTU2)

write.csv(OTUdf, "OTU_rarefied_CB_WLE.csv")

otu_final<- read.csv("OTU_rarefied_CB_WLE.csv", row.names=1) ##change name 1k--- depending on rarefaction curve
View(otu_final) ##184 samples remain, 203-19, 50421 otus , need to remove rowsums ==0
any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #5590 otus with no counts
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) #19639 otus remain after removing taxa not present in any samples, initially was 50421

write.csv(wle_cb_map2, "metadata_cb_wle.csv")

tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("metadata_cb_wle.csv", sep=",", row.names=1)
View(metadata)
tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged


phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19639 taxa and 184 samples ]
sample_data() Sample Data:       [ 184 samples by 18 sample variables ]
tax_table()   Taxonomy Table:    [ 19639 taxa by 8 taxonomic ranks ]




##barplot class level
phyloseq_class <- phyloseq_merged %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class) # Sort data frame alphabetically by phylum/family etc

phyloseq_class$Class[phyloseq_class$Abundance < 0.01] <- "< 1% abund."

View(phyloseq_class)



devtools::install_github("doehm/evoPalette")
library(evoPalette)
launch_evo_palette()

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
phyloseq_class_new<- phyloseq_class %>% mutate(transect1= ifelse(transect=="WC"|transect=="wetland", "wetland/wetland-transition edge", transect))
View(phyloseq_class_new)
plot_class<- ggplot(phyloseq_class_new, aes(x= tree_number, y = Abundance, fill =Class)) + 
  facet_grid2(region_new+site~transect1+horizon, scales="free", independent = "x", labeller = label_wrap_gen(width=10)) + theme(strip.text.x = element_text(size = 10))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(palette_new72))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=17),
                                     axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))+ guides(fill=guide_legend(ncol=2))
plot_class



ggsave(filename = "barplot_class_CB_WLE_new.jpg", plot = plot_class,
       width = 70,
       height = 35, units = c("cm"),
       dpi = 300)

##barplot phylum level
phyloseq_phylum <- phyloseq_merged %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum) # Sort data frame alphabetically by phylum/family etc

phyloseq_phylum$Phylum[phyloseq_phylum$Abundance < 0.01] <- "< 1% abund."

View(phyloseq_phylum)
palette_new56 = c("#560d0d","#a35151", "#dba4a4", "#cc1c1c","#111b77","#283dff","#636bb7",
                  "#bfc5ff","#014443","#60ffaf","#b7ffdb","#825121","#117744",
                  "#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
                  "#5b5b19","#ae09ea","#521899","#1e0047","#195637",
                  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
                  "#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD",
                  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "firebrick2",
                  "darkslategrey", "dimgrey","hotpink3", "indianred1", "indianred", "lightcyan4", "lightgreen", "magenta",
                  "lightgoldenrod1", "lightpink3", "lavenderblush1", "lightgoldenrod4", "mediumorchid1","#fa7efc","#ffb7ef", "#fcfc00","#ffff9e")

install.packages("ggh4x")

phyloseq_phylum_new<- phyloseq_phylum %>% mutate(transect1= ifelse(transect=="WC"|transect=="wetland", "wetland/wetland-transition edge", transect))
plot_phylum<- ggplot(phyloseq_phylum_new, aes(x= tree_number, y = Abundance, fill =Phylum)) + 
  facet_grid2(region_new+site~transect1+horizon, scales="free", independent = "x", labeller = label_wrap_gen(width=15)) + theme(strip.text.x = element_text(size = 10))+ 
  theme(strip.text.y = element_text(size = 18))+
  geom_bar(stat = "identity") +
  scale_fill_manual(values=as.vector(palette_new56))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_blank())+
  #
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(axis.text.y=element_text(size=17),
                                     axis.title.y=element_text(size=20)) + 
  theme(legend.title = element_text(size=20))+ theme(legend.text = element_text(size=17))+ guides(fill=guide_legend(ncol=1))
plot_phylum



ggsave(filename = "barplot_phylum_CB_WLE.jpg", plot = plot_phylum,
       width = 50,
       height = 25, units = c("cm"),
       dpi = 300)

######### Principal Coordinate analysis #########


## PCoaA analysis all samples CB and WLE

#metadata<- read.csv("sample_metadata_dna_cdna.csv") ##change metadata file as needed to analyze for plots, keep all or selected ones
#keep.samples <- as.vector(metadata$SampleID)
#keep.samples

#phyloseq_merged_final <- prune_samples(keep.samples, phyloseq_merged) #changed to merged phyloseq instead of rarefy
#phyloseq_merged_final


sample_data(phyloseq_merged)<- data.frame(sample_data(phyloseq_merged)) %>% mutate(transect1=ifelse(transect=="WC"|transect=="wetland", "wetland/wetland-transition edge", transect))

View(sample_data(phyloseq_merged))


##pcoa
phyloseq_pcoa <- ordinate(
  physeq = phyloseq_merged, 
  method = "PCoA", 
  distance = "bray"
)
##plot by site
pcoa_CB_WLE<-plot_ordination(
  physeq = phyloseq_merged,
  ordination = phyloseq_pcoa,
  color = "transect1",
  shape = "site")+
  #title = "PCoA CB_WLE 2") + 
  scale_color_manual(values = c("#CC79A7", "#0072B2", "#ea7f17", "#014443", "goldenrod","firebrick"))+
  scale_shape_manual(values=c(15, 16, 17, 18, 20, 25))+
  geom_point(aes(color = transect1), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))

pcoa_CB_WLE



ggsave(filename="pcoa-site-shape-CB-WLE.TIFF", plot=pcoa_CB_WLE, width=8, height=6, units="in", dpi=300)

##basic stats

set.seed(1)
# Calculate bray curtis distance matrix
phyloseq_bray <- phyloseq::distance(phyloseq_merged, method = "bray") ##using merged object is fine here, gives same results as rarefy object
phyloseq_bray

sample_df <- data.frame(sample_data(phyloseq_merged)) 
View(sample_df)
set.seed(1)
library(vegan)
adonis(phyloseq_bray ~ region_new*transect1*site, by="margin", data = sample_df)##gives the exact result as adonis2 by "terms"

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = phyloseq_bray ~ region_new * transect1 * site, data = sample_df, by = "terms")
Df SumOfSqs      R2       F Pr(>F)    
region_new             1    7.286 0.08770 24.4631  0.001 ***
  transect1              2    5.582 0.06719  9.3712  0.001 ***
  site                   4    7.269 0.08750  6.1018  0.001 ***
  region_new:transect1   2    4.231 0.05093  7.1036  0.001 ***
  transect1:site         8    9.266 0.11154  3.8888  0.001 ***
  Residual             166   49.441 0.59513                   
Total                183   83.076 1.00000                   
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###################################################################
###deseq for taxa that are fold enriched in one region versus the other
####################################################################

##SAME PHYLOSEQ OBJECT AS ABOVE 

phyloseq_merged
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 19639 taxa and 184 samples ]
sample_data() Sample Data:       [ 184 samples by 18 sample variables ]
tax_table()   Taxonomy Table:    [ 19639 taxa by 8 taxonomic ranks ]

##CONTINUE WITH DESEQ

head(sample_data(phyloseq_merged)$region_new, 25)

if (!requireNamespace("BiocManager", quietly = TRUE))
  
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
packageVersion("DESeq2")

##convert phyloseq format to DESEq dataset with dispersions estimated using the experimental design formula
diagdds = phyloseq_to_deseq2(phyloseq_merged, ~ region_new)
diagdds$region_new
diagdds

diagdds$region_new <-factor(diagdds$region_new, levels = c("CB", "WLE"))
diagdds$region_new #make sure that Control is the first level in the treatment factor, so that the
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


write.csv(sigtab, "sigtab_CB_WLE_enriched.depleted.CB.csv")

## Figure

plot<- ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + 
  #scale_x_reordered() +
  geom_point(size=3) + scale_color_manual(values=as.vector(P52))+ 
  ggtitle("DeSeq analysis")+
  theme(plot.title = element_text(hjust = 0.5, size=18))+
  theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1), axis.title.x = element_text(size=18))+
  theme(axis.text.y=element_text(size=13), axis.title.y=element_text(size=18)) + 
  theme(legend.title = element_text(size=18))+ theme(legend.text = element_text(size=13)) +
  guides(col=guide_legend(ncol=1)) +theme(strip.text.x = element_text(size = 30)) +geom_hline(yintercept=0, linetype='dashed', col = 'black')+
  ylab("Log2FoldChange")+coord_flip() 

ggsave(filename = "Deseq2.tiff", plot = plot,
       width = 90,
       height = 70, units = c("cm"),
       dpi = 300)

###################################
##indicator analysis may give a sound answer to which OTUs are detected as indicators in both CB and WLE because 
##it gives 1,0 ; 1,1 ; 0,1 combinations
###################################

##load indicator species analysis package

install.packages("indicspecies")
library(indicspecies)

###########

##phyloseq object
library(phyloseq)
library(tidyverse)
otu = read.csv("COMPASS_CB_WLE_count.csv", sep=",", row.names=1) ##USING UNRAREFIED DATA, CAN ALSO CHECK WITH RAREFIED 

otu<- otu %>% select(!c(s20230217096,s20230217120,s20230217121,s20230217122, s20230217123, s20230217124))

tax = read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
tax = as.matrix(tax)

metadata = read.csv("metadata_cb_wle.csv", sep=",", row.names=1) # use for crop indicators 


OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 50421 taxa and 203 samples ]
sample_data() Sample Data:       [ 203 samples by 18 sample variables ]
tax_table()   Taxonomy Table:    [ 50421 taxa by 8 taxonomic ranks ]

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
  GetIndicators(dataframe=phyloseq_merged, "region_new")


##checking

View(ind_region)

ind_region_subset<- ind_region %>% subset(s.WLE==1 & s.CB==1)
View(ind_region_subset) ##no common indicators across two regions


library(tidyverse)

write.csv(ind_region, "indicator_region_CB_WLE.csv")

##creating a dataframe for CB only indicators and passing on as a join to the WLE count table. 

##need to use unrarefied count tables for WLE here for the join

ind_region<- read.csv("indicator_region_CB_WLE.csv", row.names=1)
View(ind_region)
ind_region_CB<- ind_region %>% subset(s.WLE==0 & s.CB==1)
View(ind_region_CB)


otu_all<- read.csv("COMPASS_CB_WLE_count.csv", row.names=1)
View(otu_all)
otu_WLE<- otu_all %>% select(!starts_with("s"))
View(otu_WLE)
otu_WLE<- otu_WLE %>% select(!starts_with("Zymo"))

otu_WLE_final = otu_WLE[rowSums(otu_WLE[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_WLE_final)

merge_cbind_wle<- merge(ind_region_CB, otu_WLE_final, by=0)
View(merge_cbind_wle)

write.csv(merge_cbind_wle, "CBind_inWLEunrarefiedCount.csv")

##remove extra columns 

merge_cbind_wle <- as.data.frame(merge_cbind_wle)
rnames <-merge_cbind_wle[,1] 
rownames(merge_cbind_wle) <- rnames  # assign row names
View(merge_cbind_wle)
merge_cbind_wle_new <- merge_cbind_wle %>% select(!c(Row.names, OTU, s.CB, s.WLE, index, stat, sequencefeature, p.value, Domain, Phylum, Class, Order, Family, Genus, Species)) ##add other columns
View(merge_cbind_wle_new)

write.csv(merge_cbind_wle_new, "indicatorCB_inWLE_count_final.csv")

##get these otus from WLE to assess treatment variation 
##merge to phyloseq object



tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("WLE_metadata.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(merge_cbind_wle_new, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 106 taxa and 85 samples ]
sample_data() Sample Data:       [ 85 samples by 5 sample variables ]
tax_table()   Taxonomy Table:    [ 106 taxa by 8 taxonomic ranks ]



##max_standardize 
max_standardize_WLE <-decostand(merge_cbind_wle_new, method = "max", MARGIN = 1) 
View(max_standardize_WLE)

tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("WLE_metadata.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(max_standardize_WLE, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 106 taxa and 85 samples ]
sample_data() Sample Data:       [ 85 samples by 5 sample variables ]
tax_table()   Taxonomy Table:    [ 106 taxa by 8 taxonomic ranks ]
  
phyloseq_otu <- phyloseq_merged %>%
#tax_glom(taxrank = "Class")   %>% 
#transform_sample_counts(function(x) {x/sum(x)} ) %>%
psmelt()
View(phyloseq_otu)

##combining class and otu information

phyloseq_otu$ASVclass <- paste(phyloseq_otu$Class, phyloseq_otu$OTU)
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
#phyloseq_otu_new<- phyloseq_otu %>% mutate(transect1= ifelse(transect=="WC", "wetland", transect))
#View(phyloseq_otu_new)

phyloseq_otu$transect <- factor(phyloseq_otu$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)
plot<- ggplot(phyloseq_otu, aes(x= transect, y = Abundance, fill =transect)) + 
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



ggsave(filename = "ASV_CBind_inWLE.jpg", plot = plot,
       width = 120,
       height = 40, units = c("cm"),
       dpi = 300)


##plot the most abundant ones only 

merge_cbind_wle_new<- read.csv("indicatorCB_inWLE_count_final.csv", row.names=1)
merge_cbind_wle_new<-merge_cbind_wle_new[order(rowSums(merge_cbind_wle_new), decreasing = TRUE),]
View(merge_cbind_wle_new)
merge_cbind_wle_new_20 <- merge_cbind_wle_new[c(1:20),] ##select top 20 most abundant Classes
View(merge_cbind_wle_new_20)

write.csv(merge_cbind_wle_new_20, "CBind_inWLE_top20.csv")

##max_standardize 
max_standardize_WLE <-decostand(merge_cbind_wle_new_20, method = "max", MARGIN = 1) ##from Jackson's phil trans paper
View(max_standardize_WLE)

tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("WLE_metadata.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(max_standardize_WLE, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 20 taxa and 85 samples ]
sample_data() Sample Data:       [ 85 samples by 5 sample variables ]
tax_table()   Taxonomy Table:    [ 20 taxa by 8 taxonomic ranks ]

phyloseq_otu <- phyloseq_merged %>%
  #tax_glom(taxrank = "Class")   %>% 
  #transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()
View(phyloseq_otu)

##combining class and otu information

phyloseq_otu$ASVclass <- paste(phyloseq_otu$Class, phyloseq_otu$OTU)
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
#phyloseq_otu_new<- phyloseq_otu %>% mutate(transect1= ifelse(transect=="WC", "wetland", transect))
#View(phyloseq_otu_new)

phyloseq_otu$transect <- factor(phyloseq_otu$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)
plot<- ggplot(phyloseq_otu, aes(x= transect, y = Abundance, fill =transect)) + 
  facet_grid2(site~ASVclass, scales="free", independent = "x", labeller = label_wrap_gen(width=25)) + theme(strip.text.x = element_text(size = 18, face="bold"))+ 
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



ggsave(filename = "ASV_CBind_inWLE1.jpg", plot = plot,
       width =45,
       height = 10, units = c("in"),
       dpi = 300)


################################
##analyzing CB for specific taxa 
#################################

otu_final<- read.csv("OTU_rarefied_CB.csv", row.names=1) ##change name 1k--- depending on rarefaction curve
View(otu_final)
any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #5590 otus with no counts
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) #44831 otus 


library(phyloseq)


tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("combined_metadataCB.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged


phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 44831 taxa and 116 samples ]
sample_data() Sample Data:       [ 116 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 44831 taxa by 8 taxonomic ranks ]

###deseq analysis for one site may not be feasible, need to pick two variables of interest

##simper might be helpful. indicator could be better. 
##simper will give taxa contributing to the most variation across pairs of variables

#########################
##########simper CB #####
#########################

physeq_class_all <- phyloseq_merged %>%
  tax_glom(taxrank = "Class")   %>% 
  #transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()
# agglomerate at family/phylum level

View(physeq_class_all)

physeq_class<- physeq_class_all %>% select("Abundance", "Class", "Sample")

#write.csv(physeq_class, "physeq_class_simper.csv")

##Class"unclutured" has 2 values for each sample, from Desulfobacterota and Armatimonadota, so removed uncultured row for simper analysis 
##to do this I first had to rename uncultured in the csv file to _uncultured other wise it wasnt working ##i made an error with phyloseq object names so was not working, no need to rename to _uncultured, saves a step

#physeq_class<- read.csv("physeq_class_simper.csv")

physeq_class <- physeq_class %>% filter(Class!=" uncultured")

View(physeq_class)
physeq_wide_class<- pivot_wider(physeq_class, names_from=Sample, values_from=Abundance)

View(physeq_wide_class)

##identify duplicates (this code showed that the uncultured row was the problem)

physeq_wide_duplicates<- physeq_class %>%
  dplyr::group_by(Class, Sample) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 
View(physeq_wide_duplicates)

write.csv(physeq_wide_class, "physeq_wide_class_simper.csv") ##did sum calculations to see read counts per sample


physeq_wide_class_mat <- as.data.frame(physeq_wide_class)
rnames <-physeq_wide_class_mat[,1] 
rownames(physeq_wide_class_mat) <- rnames  # assign row names
View(physeq_wide_class_mat)
physeq_wide_class_new <- physeq_wide_class_mat %>% select(!Class)
View(physeq_wide_class_new)

physeq_wide_class_trans <- t(physeq_wide_class_new)
View(physeq_wide_class_trans)

##instead of class level lets do OTU level

View(otu_final)
physeq_wide_class_trans <- t(otu_final)
View(physeq_wide_class_trans)


#physeq_wide_class_trans <- cbind(SampleID = rownames(physeq_wide_class_trans), physeq_wide_class_trans, row.names=FALSE)
#View(physeq_wide_class_trans)

#rownames(physeq_wide_class_trans) <- NULL

#metadata<- read.csv("metadata_simper.csv", row.names=1) ##cannot use this as sample numbers dont match has to be 212 in metadata also
#View(metadata)
library(vegan)

sample_data = data.frame(sample_data(phyloseq_merged))

View(sample_data)

sample_data<- sample_data %>% mutate(transect1=ifelse(transect=="WC", "wetland", transect))

View(sample_data)

# running the simper analysis on the dataframe and the variable of interest 
transect_simper <- simper(physeq_wide_class_trans, sample_data$transect1, permutations = 100) ##or by site


sink(file="simper_transect_CB_output.txt")
# printing the top OTUs
print(transect_simper)
sink(file=NULL)


##plotting important taxa 
##first do relative abundance on all samples at OTU level. 
#left join taxa from simper excel file to the OTU table rarefied CB and select the top 10 taxa that contribute the most variation. 
##plot relative abundance as boxchart across upland, wetland, transition for each OTU

df1<- read.csv("transition_wetland_CB_simper.csv")

df2<- read.csv("upland_transition_CB_simper.csv")

df3<- read.csv("upland_wetland_CB_simper.csv")

View(df1)
View(df2)
View(df3)

df_merge<- inner_join(df1, df2, by="ASV")

View(df_merge)


df_merge_edited<- inner_join(df_merge, df3, by="ASV")

View(df_merge_edited)

write.csv(df_merge_edited, "simper_CB_bytransect.csv")


View(otu_final)


df_merge_edited<- read.csv("simper_CB_bytransect.csv", row.names=2)
otu_final_simper<- merge(otu_final, df_merge_edited, by=0)
View(otu_final_simper)


rnames <-otu_final_simper[,1] 
rownames(otu_final_simper) <- rnames  # assign row names
View(otu_final_simper)
otu_final_simper_new <- otu_final_simper %>% select(!c(X, perc.x, perc.y, perc, Row.names))
View(otu_final_simper_new)

##max_standardize 
max_standardize_simper_CB <-decostand(otu_final_simper_new, method = "max", MARGIN = 1) ##from Jackson's phil trans paper
View(max_standardize_simper_CB)

tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("combined_metadataCB.csv", sep=",", row.names=1)

tax = as.matrix(tax) 

otu<- otu_table(max_standardize_simper_CB, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 16 taxa and 116 samples ]
sample_data() Sample Data:       [ 116 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 16 taxa by 8 taxonomic ranks ]
> 

phyloseq_otu <- phyloseq_merged %>%
  #tax_glom(taxrank = "Class")   %>% 
  #transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()
View(phyloseq_otu)

##combining class and otu information

phyloseq_otu$ASVclass <- paste(phyloseq_otu$Class, phyloseq_otu$OTU)
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



ggsave(filename = "ASV_simper_CB_1.jpg", plot = plot,
       width = 120,
       height = 40, units = c("cm"),
       dpi = 300)






#########################################################
################Alpha Diversity #########################
#########################################################


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
  r <- rarefy_even_depth(phyloseq_merged, sample.size = 40000, verbose = FALSE, replace = TRUE) #no rarefaction a second time.
  
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed"))) ##changed r to phyloseq_merged, no second rarefaction
  richness[ ,i] <- rich
  
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson"))) ##changed r to phyloseq_merged, no second rarefaction
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

s <- read.csv("metadata_CB_alphadiv.csv")
alphadiv <- merge(alpha, s, by = "SampleID") 
write.csv(alphadiv,file = "alphadiv_CB.csv")


##subset if needed
alphadiv_wle<-alphadiv[alphadiv$region=="WLE",]
dim(alphadiv_wle)
alphadiv_cb<-alphadiv[alphadiv$region=="CB",]
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
alphadiv_CB <- read.csv("alphadiv_CB.csv")
pd<-position_dodge(0.7)
alphadiv_CB_InvS =  alphadiv_CB %>% filter(measure=="Inverse Simpson")

View(alphadiv_CB_InvS)


#### Richness CB
alphadiv_CB <- read.csv("alphadiv_CB.csv")
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
alphadiv_rich <- read.csv("alphadiv_CB.csv") ## combined richness and evenness estimates of all samples


pd<-position_dodge(0.7)


alphadiv_rich =  alphadiv_rich %>% filter(measure=="Inverse Simpson")
alphadiv_rich =  alphadiv_rich %>% filter(measure=="Richness")

##Conduct the following code for each of the datasets selected above
View(alphadiv_rich)

alphadiv_rich<- alphadiv_rich %>% mutate(transect1=ifelse(transect=="WC", "wetland", transect))
options(contrasts = c("contr.sum", "contr.poly")) ### must set this before running the model for type III tests
model = lm(mean~ transect1+site + transect1:site +transect1 + transect1:site,
           data=alphadiv_rich)

shapiro.test(resid(model)) #W statistic > 0.9 means normality okay #check normality to flag outliers 

##inverse simpson with horizon in model
Shapiro-Wilk normality test

data:  resid(model)
W = 0.92255, p-value = 4.761e-06

##inverse simpson without horizon in model
Shapiro-Wilk normality test

data:  resid(model)
W = 0.89548, p-value = 1.715e-07

##richness without horizon in model
Shapiro-Wilk normality test

data:  resid(model)
W = 0.97544, p-value = 0.03151





ee= as.matrix(resid(model))
qqnorm (ee) #look for outliers (far away from curve) and test for normality 
ee

library(car)


#options(contrasts = c("contr.sum", "contr.poly"))
Anova(model, type=c("III"))

####Results inverse simpson, works without horizon as horizon levels are unequal among transects and sites
Anova Table (Type III tests)

Response: mean
Sum Sq  Df  F value    Pr(>F)    
(Intercept)    3529996   1 470.4282 < 2.2e-16 ***
  transect1       447670   2  29.8296  5.06e-11 ***
  site             90257   2   6.0141 0.0033478 ** 
  transect1:site  157552   4   5.2491 0.0006716 ***
  Residuals       802906 107                       
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



##Richness
Anova Table (Type III tests)

Response: mean
Sum Sq  Df   F value    Pr(>F)    
(Intercept)    164360079   1 2162.0538 < 2.2e-16 ***
  transect1       13514325   2   88.8862 < 2.2e-16 ***
  site              259164   2    1.7046    0.1868    
transect1:site   3078004   4   10.1223 5.481e-07 ***
  Residuals        8134177 107                        
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
  geom_text(data=summarized, aes(x=summarized$transect1, y =2500, 
                                 group=summarized$transect1, label = summarized$group, vjust=-0.3), 
            check_overlap = T, 
            position = position_dodge(width = 0.9), col = "black") + 
  ylab("Richness (Number of Observed ASVs)") +
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

ggsave(filename = "alphadiv_letters_CB_InvSimp.tiff", plot = p,
       width = 15,
       height = 15, units = c("cm"),
       dpi = 300)

ggsave(filename = "alphadiv_letters_CB_Richness_ASV.tiff", plot = p,
       width = 15,
       height = 15, units = c("cm"),
       dpi = 300)







#######################################
###CAP analysis with soil chemistry data
#######################################


##sample names check with microbiome and soil chemistry data

data<- read.csv("allhorizon_soilchem.csv")
data_cb<- data %>% filter(region=="CB")

unique(data_cb$SampleID) ##125 SAMPLES

data<- read.csv("allhorizon_soilchem.csv")
data_wle<- data %>% filter(region=="WLE" & horizon=="A")

unique(data_wle$SampleID) ##95 SAMPLES


##reading in new phyloseq object with all samples in CB
otu_final<- read.csv("OTU_rarefied_CB.csv", row.names=1) ##change name to WLE or CB as intended
View(otu_final)
any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #5590 otus with no counts
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) #44831 otus for CB, 3257 for WLE

data_chem<- read.csv("allhorizon_soilchem.csv")
View(data_chem)
metadata_cb_wle<- read.csv("metadata_cb_wle.csv")
View(metadata_cb_wle)
metadata_cap<- merge(data_chem, metadata_cb_wle, by="SampleID")

View(metadata_cap)

unique(metadata_cap$SampleID) #236 samples, 2 are missing from soil chem

samples_lost<- anti_join(metadata_cb_wle, metadata_cap, by = "SampleID")
View(samples_lost) ##COMPASS.Dec2021.104 and COMPASS.Dec2021.105


##subset to CB samples

cap_cb<- metadata_cap %>% filter(region.x=="CB")

View(cap_cb)
##subset to WLE samples
cap_wle<- metadata_cap %>% filter(region.x=="WLE")


View(cap_wle)


##melt to wide format 

cap_cb_new<- cap_cb%>% dplyr::select(!analysis) %>% pivot_wider(names_from="name", values_from="value")  
View(cap_cb_new)
cap_wle_new<- cap_wle%>% dplyr::select(!analysis) %>% pivot_wider(names_from="name", values_from="value")  
View(cap_wle_new)

write.csv(cap_cb_new, "cap_cb_new.csv") ##did not save this a second time

write.csv(cap_wle_new, "cap_wle_new.csv") 

library(phyloseq)


tax<- read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
metadata <- read.csv("cap_wle_new.csv", sep=",", row.names=8)

##removing columns that are >50% NAs, adjust separately for CB and WLE
##WLE
metadata_new<- metadata %>% dplyr::select(!c(region.y, Description, #Nitrate_meq100g, #Sulfate_meq100g, Phosphate_meq100g, #Bromide_meq100g, Ammonia_meq100g
                                      Sample.ID, Sample.ID.1, project, date, DNA.Well, DNA.Plate, PCR.Plate, PCR.Well, Index.ID, Barcode, Potassium_meq100g, dic_ugg, Sodium_meq100g, Fluoride_meq100g, Chloride_meq100g, Magnesium_meq100g, Calcium_meq100g,
                                      percentOM)) 
#CB
metadata_new<- metadata %>% dplyr::select(!c(region.y, Description, #Nitrate_meq100g, #Sulfate_meq100g, #Bromide_meq100g, Ammonia_meq100g, Phosphate_meq100g,
                                              Potassium_meq100g, dic_ugg, Sodium_meq100g, Fluoride_meq100g, Chloride_meq100g, Magnesium_meq100g, Calcium_meq100g,
                                             percentOM)) 


View(metadata_new)

##removing rows that are NAs, need to do this because ordinate function does not work with NAs
metadata_new=na.omit(metadata_new)

#metadata_new<- na.pass(metadata)
View(metadata_new)
tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(metadata_new)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

##without na samples being removed CB
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 44831 taxa and 116 samples ]
sample_data() Sample Data:       [ 116 samples by 42 sample variables ]
tax_table()   Taxonomy Table:    [ 44831 taxa by 8 taxonomic ranks ]


##after NA samples are removed CB, rows and columns, first try where many columns were removed retaining maximum samples
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 44831 taxa and 85 samples ]
sample_data() Sample Data:       [ 85 samples by 40 sample variables ]
tax_table()   Taxonomy Table:    [ 44831 taxa by 8 taxonomic ranks ]

##keeping more columns in cb, BUT removing rows(samples, using na.omit) that dont have that data so sample number lowers
##two samples less than the chemistry data file. probably coz the OTU tables does not have those two samples after rarefaction, but they were present in the chemistry data

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 44831 taxa and 65 samples ]
sample_data() Sample Data:       [ 65 samples by 44 sample variables ]
tax_table()   Taxonomy Table:    [ 44831 taxa by 8 taxonomic ranks ]
> 

##after NA rows and columns are removed, for WLE

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 3257 taxa and 57 samples ]
sample_data() Sample Data:       [ 57 samples by 33 sample variables ]
tax_table()   Taxonomy Table:    [ 3257 taxa by 8 taxonomic ranks ]

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

##CB first trial
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 44831 taxa and 85 samples ]
sample_data() Sample Data:       [ 85 samples by 40 sample variables ]
tax_table()   Taxonomy Table:    [ 44831 taxa by 8 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 44831 tips and 44830 internal nodes ]



##CB second trial 
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 44831 taxa and 65 samples ]
sample_data() Sample Data:       [ 65 samples by 44 sample variables ]
tax_table()   Taxonomy Table:    [ 44831 taxa by 8 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 44831 tips and 44830 internal nodes ]




##WLE
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 3257 taxa and 57 samples ]
sample_data() Sample Data:       [ 57 samples by 33 sample variables ]
tax_table()   Taxonomy Table:    [ 3257 taxa by 8 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 3257 tips and 3256 internal nodes ]



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
cap.ord.df = data.frame(vegan::scores(cap_ord, display="sites")) %>%
  tibble::rownames_to_column(var="Sample.ID.1") %>%
  dplyr::select(Sample.ID.1, CAP1, CAP2) %>%
  left_join(sample_data(phyloseq_merged), by = "Sample.ID.1")

# CAP plot WLE
## For this I use my own code rather than phyloseq plotting for more control. To do this I have to pull out the x (CAP1) and y (CAP2) axes from the ordination and add in my metadata
cap.ord.df = data.frame(vegan::scores(cap_ord, display="sites")) %>%
  tibble::rownames_to_column(var="SampleID") %>%
  dplyr::select(SampleID, CAP1, CAP2) %>%
  left_join(sample_data(phyloseq_merged), by = "SampleID")

## Calculate eignevalues and fraction of variation explained by each CAP axis. I use this for the axis labels in the plot
eigvec = vegan::eigenvals(cap_ord)
fracvar = round(eigvec/sum(eigvec)*100, 2)

## Plot initial figure of points
library(ggplot2)
cap_plot<- ggplot(data=cap.ord.df, aes(x=CAP1, y=CAP2)) +
  geom_point(aes(fill=transect.y, shape=site.x)) +
 scale_shape_manual(values=c("PR" = 24, "CC" = 21, "OWC" = 22)) + ##change site names for CB
 # scale_fill_gradient(low="white", high="black") +
  labs(x=paste("CAP1 (", fracvar[1], "%)", sep=""),y=paste("CAP2 (", fracvar[2], "%)", sep=""))

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
  guides(shape = guide_legend(order = 1),
         fill = guide_legend(order = 2, override.aes=list(shape=22), ncol=2))

cap.plot
ggsave(filename = "cap_wle_new.jpg", plot = cap.plot,
       width = 17,
       height = 22, units = c("cm"),
       dpi = 300)

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
otu = read.csv("OTU_rarefied_CB.csv", sep=",", row.names=1) ##USING RAREFIED DATA here unlike unrarefied data used for site specific indicators  

otu_final = otu[rowSums(otu[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final)



#otu<- otu %>% select(starts_with("s"))
#otu<- otu %>% select(!c(s20230217096,s20230217120,s20230217121,s20230217122, s20230217123, s20230217124))
View(otu_final)
tax = read.csv("COMPASS_CB_WLE_taxonomy_edited.csv", sep=",", row.names=1)
tax = as.matrix(tax)

metadata = read.csv("combined_metadataCB.csv", sep=",", row.names=1) 


OTU = otu_table(otu_final, taxa_are_rows = TRUE)
TAX = tax_table(tax)
meta = sample_data(metadata)

##merge
phyloseq_merged = phyloseq(OTU, TAX, meta)
phyloseq_merged

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 44831 taxa and 116 samples ]
sample_data() Sample Data:       [ 116 samples by 15 sample variables ]
tax_table()   Taxonomy Table:    [ 44831 taxa by 8 taxonomic ranks ]

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

sample_data(phyloseq_merged)<- data.frame(sample_data(phyloseq_merged)) %>% mutate(transect1=ifelse(transect=="WC", "wetland", transect))

View(sample_data(phyloseq_merged))

sample_data(phyloseq_merged)$transect1 <- factor(sample_data(phyloseq_merged)$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)


ind_transectCB <-
  GetIndicators(dataframe=phyloseq_merged, "transect1")


##checking

View(ind_transectCB)

ind_transectCB_subset<- ind_transectCB %>% subset(s.upland==0 & s.transition==0 & s.wetland==1)
View(ind_transectCB_subset) ##2770 Taxa uniquely associated with wetland

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



##nannodrop data

data<- read.csv("nanodrop_wle.csv")

data1<- read.csv("WLE_metadata.csv")
merge<- left_join(data1, data, by="SampleID")
merge<- as.data.frame(merge)
class(merge)

library(tidyverse)
View(merge)
merge
library(ggplot2)

merge$transect <- factor(merge$transect,levels = c("upland", "transition", "wetland/wetland-transition edge"), ordered=TRUE)
plot<- ggplot(merge, aes(x=transect ,y=concentration, color=site)) + theme_classic()+ scale_color_manual(values=c("#39568CFF", "#55C667FF", "orange")) +
  #scale_shape_manual(values=c(1,16))+ #scale_size(range = c(1,10), # point size range
  #label = scales::dollar)+
  #expand_limits(shape = seq(2, 6, by = 1))+
  #scale_alpha_manual(values = c(2,3,4,5,6))+
  geom_boxplot()+
  geom_jitter() +  #scale_size_manual(values = c(6))+ 
  #geom_smooth(method = lm, se=FALSE, aes(group=summarise), inherit.aes=TRUE) + 
theme(
    strip.text.x = element_text(
      size = 16
    ))+
  theme(plot.title = element_text(size= 18, hjust = 0.5))+
  scale_x_discrete(name = 'transect', 
                   breaks = c('upland', 'transition', 'wetland/wetland-transition edge'), 
                   labels = c('upland', 'transition', 'wetland/wetland-\ntransition edge'))+
  labs(y=expression(paste(''*mu~'g DNA ml'^-1*'nucleic acid')), color='site')+ xlab("transect")+theme(axis.text.x = element_text(colour = "black", size = 17,angle = 90, vjust = 0.5, hjust=1))+
  theme( axis.text.y = element_text(size = 16, colour = "black"), 
         axis.title.y= element_text(size=16, color="black"),
         axis.title.x= element_text(size=16, color="black"),
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =14, colour = "black"),
         legend.text = element_text(size = 12, colour = "black")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")), color= guide_legend(order = 1)) 


plot

ggsave(filename = "DNA_con_WLE.tiff", plot = plot,
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)


##nanodrop cb
library(tidyverse)

library(ggplot2)

data<- read.csv("nanodrop_cb.csv", header=TRUE)

data1<- read.csv("combined_metadataCB.csv")
merge<- left_join(data, data1, by="SampleID")
merge<- as.data.frame(merge)
class(merge)
View(merge)



merge_new<- merge %>% mutate(transect1=ifelse(transect=="WC", "wetland", transect))

merge_new$transect1 <- factor(merge_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)
View(merge_new)
plot<- ggplot(merge_new, aes(x=transect1 ,y=concentration, color=site)) + theme_classic()+ scale_color_manual(values=c("#39568CFF", "#55C667FF", "orange")) +
  #scale_shape_manual(values=c(1,16))+ #scale_size(range = c(1,10), # point size range
  #label = scales::dollar)+
  #expand_limits(shape = seq(2, 6, by = 1))+
  #scale_alpha_manual(values = c(2,3,4,5,6))+
  geom_boxplot()+
  geom_jitter() +  #scale_size_manual(values = c(6))+ 
  #geom_smooth(method = lm, se=FALSE, aes(group=summarise), inherit.aes=TRUE) + 
  theme(
    strip.text.x = element_text(
      size = 16
    ))+
  theme(plot.title = element_text(size= 18, hjust = 0.5))+
  scale_x_discrete(name = 'transect', 
                   breaks = c('upland', 'transition', 'wetland'), 
                   labels = c('upland', 'transition', 'wetland'))+
  labs(y=expression(paste(''*mu~'g DNA ml'^-1*'nucleic acid')), color='site')+ xlab("transect")+theme(axis.text.x = element_text(colour = "black", size = 17,angle = 90, vjust = 0.5, hjust=1))+
  theme( axis.text.y = element_text(size = 16, colour = "black"), 
         axis.title.y= element_text(size=16, color="black"),
         axis.title.x= element_text(size=16, color="black"),
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"), 
         legend.title = element_text(size =14, colour = "black"),
         legend.text = element_text(size = 12, colour = "black")) +
  guides(linetype = guide_legend(override.aes= list(color = "black")), color= guide_legend(order = 1)) 


plot

ggsave(filename = "DNA_con_CB.tiff", plot = plot,
       width = 20,
       height = 15, units = c("cm"),
       dpi = 300)





