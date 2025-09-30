
###calculating relative abundances of FTICR classes

library(tidyverse)
icr_data_long<- read.csv("icr_long_all_samples.csv")
icr_meta<- read.csv("icr_meta.csv")
sample_key<- read.csv("sample-metadata-cb-wle.csv")

##run this function to calculate icr relative abundances
compute_icr_relabund = function(icr_data_long, icr_meta){
  
  icr_data_long %>% 
    # filter(horizon != "B") %>% 
    # add the Class column to the data
    left_join(dplyr::select(icr_meta, formula, Class), by = "formula") %>% 
    # calculate abundance of each Class as the sum of all counts
    group_by(sample_label, Class) %>%
    dplyr::summarise(abund = sum(presence)) %>%
    ungroup %>% 
    # create a new column for total counts per core assignment
    # and then calculate relative abundance  
    group_by(sample_label) %>% 
    dplyr::mutate(total = sum(abund),
                  relabund  = ((abund/total)*100))
}


##run the function on data
icr_relabund<- compute_icr_relabund(icr_data_long, icr_meta)
write.csv(icr_relabund, "icr_relabun.csv")
unique(icr_relabund$sample_label) #238

##merge icr rel abund with metadata

metadata<- read.csv("sample-metadata-cb-wle-cap-fticr.csv")
#metadata<- metadata %>% mutate(sample_label=Sample.description)
metadata_merge<- icr_relabund %>% left_join(metadata, by="sample_label") ##this metadata_merge has way more samples (around 238) than what matches with metadata (203) because it keeps all icr data for which there is no metadata

#TESTING
test<- merge(icr_relabund, metadata, by="sample_label")
unique(test$sample_label)#203 samples

##convert to wide format for easier plotting

metadata_merge_wide<- metadata_merge %>% dplyr::select(!c(total, abund)) %>% pivot_wider(names_from="Class", values_from="relabund")  


##sample names/numbers check with microbiome and soil chemistry data
metadata_merge_wide_cb<- metadata_merge_wide %>% filter(region=="CB") #116 samples

metadata_merge_wide_wle<- metadata_merge_wide %>% filter(region=="WLE") # 87 samples

write.csv(metadata_merge_wide_cb, "cap_cb_fticr.csv")
write.csv(metadata_merge_wide_wle, "cap_wle_fticr.csv")

#######################################
###CAP analysis with FTCIR data
#######################################

##reading in tables for WLE and CB
otu_final<- read.csv("OTU_CB_rarefied20k.csv", row.names=1) ##change name to WLE or CB as intended
View(otu_final)
any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #20954 otus with no counts (wle), 15252 otus with no counts (cb)
View(otu1)

max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(otu_final) #19262 otus for WLE, 24964 otus for CB

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
metadata <- read.csv("cap_cb_fticr.csv", sep=",", row.names=3) ##change as needed for cb and wle

#metadata_test<- read.csv("cap_cb_fticr.csv", sep=",")
#otu_test<- read.csv("otu_sample_cb.csv")
#anti<- metadata_test %>% anti_join(otu_test, by="sample.id")


##removing columns that are >50% NAs, adjust separately for CB and WLE
##WLE
#metadata_new<- metadata %>% dplyr::select(!c(region.y, Sample.description, #Nitrate_meq100g, #Sulfate_meq100g, Phosphate_meq100g, #Bromide_meq100g, Ammonia_meq100g
                                            # tree_number.y, transect.y, site.y, horizon.y, project, date, DNA.Well, DNA.Plate, PCR.Plate, PCR.Well, Index.ID, Barcode, Potassium_meq100g, dic_ugg, Sodium_meq100g, Fluoride_meq100g, Chloride_meq100g, Magnesium_meq100g, Calcium_meq100g,
                                             #percentOM, gwc_perc)) 
#CB
#metadata_new<- metadata %>% dplyr::select(!c(region.y, Sample.description, #Nitrate_meq100g, #Sulfate_meq100g, #Bromide_meq100g, Ammonia_meq100g, Phosphate_meq100g,
                                       #      tree_number.y, transect.y, site.y, horizon.y, project, date, DNA.Well, DNA.Plate, PCR.Plate, PCR.Well, Index.ID, Barcode,Potassium_meq100g, dic_ugg, Sodium_meq100g, Fluoride_meq100g, Chloride_meq100g, Magnesium_meq100g, Calcium_meq100g,
                                         #    percentOM, gwc_perc)) 


#View(metadata_new)

##removing rows that are NAs, need to do this because ordinate function does not work with NAs
metadata_new=na.omit(metadata)


##changing site names for Erie and Chesapeake as needed 
metadata_new<- metadata_new%>% mutate(site_new=recode(site, "PR"="PTR", "CC"="CRC"))
metadata_new<- metadata_new%>% mutate(site_new=recode(site, "GCREW"="GCW"))



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
otu_table()   OTU Table:         [ 19262 taxa and 81 samples ]
sample_data() Sample Data:       [ 81 samples by 14 sample variables ]
tax_table()   Taxonomy Table:    [ 19262 taxa by 8 taxonomic ranks ]

## after NA are removed for CB

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 24964 taxa and 109 samples ]
sample_data() Sample Data:       [ 109 samples by 14 sample variables ]
tax_table()   Taxonomy Table:    [ 24964 taxa by 8 taxonomic ranks ]


##using bray curtis instead of weighted unifrac to match the PCoA plots for microbial community data
set.seed(1)
phyloseq_bray <- phyloseq::distance(phyloseq_merged, method = "bray")
phyloseq_bray

library("ape")
set.seed(1)
#tree<- read.tree("tree_cb_wle_full_rooted.nwk")
#random_tree = rtree(ntaxa(phyloseq_merged), rooted=TRUE, tip.label=taxa_names(phyloseq_merged))
#plot(random_tree)
#phyloseq_tree = merge_phyloseq(phyloseq_merged, random_tree)
#phyloseq_tree = merge_phyloseq(phyloseq_merged, tree)
#phyloseq_tree

##CB 
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 24964 taxa and 66 samples ]
#sample_data() Sample Data:       [ 66 samples by 29 sample variables ]
#tax_table()   Taxonomy Table:    [ 24964 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 24964 tips and 24962 internal nodes ]


##WLE


#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 19262 taxa and 71 samples ]
#sample_data() Sample Data:       [ 71 samples by 29 sample variables ]
#tax_table()   Taxonomy Table:    [ 19262 taxa by 8 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 19262 tips and 19261 internal nodes ]


#calculate weighted unifrac distance matrix
#phyloseq_wunifrac <- phyloseq::distance(phyloseq_tree, method = "wunifrac" )


##CAP plot

# First select the variables that significantly explain variation in community composition.
## To do this I'm using the ordistep function which steps through the variables removing ones that don't add to the significance. 
## There are better, more thurough explanations of that online.
### First set new names for all variables that fit nicer on the final figure.
##CB only
var_goodnames = data.frame(labels = c("aliphatic", "aromatic",
                                      "condensed.aromatic", "unsaturated.lignin"),
                           goodnames = c("aliphatic", "aromatic",
                                         "condensed_aromatic", "unsaturated_lignin"))
#WLE only
var_goodnames = data.frame(labels = c("aliphatic", "aromatic",
                                      "condensed.aromatic", "unsaturated.lignin"),
                           goodnames = c("aliphatic", "aromatic",
                                         "condensed_aromatic", "unsaturated_lignin"))



### Full model with all variables
set.seed(4242)
cap_ord.full <- ordinate(physeq = phyloseq_merged, method = "CAP", distance=phyloseq_bray, #distance = phyloseq_wunifrac, 
                         formula = ~ aliphatic+ aromatic+
                         condensed.aromatic+unsaturated.lignin)
### Null model with no variables
set.seed(4242)
cap_ord.null <- ordinate(physeq = phyloseq_merged, method = "CAP", distance=phyloseq_bray, #distance = phyloseq_wunifrac, 
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
cap_ord <- ordinate(physeq = phyloseq_merged, method = "CAP", distance = phyloseq_bray, formula = goodform)

# CAP plot CB
## For this I use my own code rather than phyloseq plotting for more control. To do this I have to pull out the x (CAP1) and y (CAP2) axes from the ordination and add in my metadata
key_match<- read.csv("key_match_cb_fticr.csv")
View(key_match)
cap.ord.df1 <- data.frame(vegan::scores(cap_ord, display="sites")) %>%
  tibble::rownames_to_column(var="sample.id") %>%
  dplyr::select(sample.id,CAP1, CAP2) %>% 
  left_join(key_match, by="sample.id")
cap.ord.df = cap.ord.df1 %>%
  left_join(sample_data(phyloseq_merged), by = "sample_label")

# CAP plot WLE
## For this I use my own code rather than phyloseq plotting for more control. To do this I have to pull out the x (CAP1) and y (CAP2) axes from the ordination and add in my metadata

key_match<- read.csv("key_match_wle_fticr.csv")
View(key_match)
cap.ord.df1 <- data.frame(vegan::scores(cap_ord, display="sites")) %>%
  tibble::rownames_to_column(var="sample.id") %>%
  dplyr::select(sample.id,CAP1, CAP2) %>% 
  left_join(key_match, by="sample.id")
cap.ord.df = cap.ord.df1 %>%
 left_join(sample_data(phyloseq_merged), by = "sample_label")

View(cap.ord.df)
cap.ord.df$transect <- gsub(fixed("wte"), "wetland", cap.ord.df$transect)
cap.ord.df$transect <- gsub(fixed("wc"), "wetland", cap.ord.df$transect)

## Calculate eignevalues and fraction of variation explained by each CAP axis. I use this for the axis labels in the plot
eigvec = vegan::eigenvals(cap_ord)
fracvar = round(eigvec/sum(eigvec)*100, 2)

## Plot initial figure of points
library(ggplot2)
cap.ord.df$transect <- factor(cap.ord.df$transect,levels = c("upland", "transition", "wetland"), ordered=TRUE)



#cap.ord.df<- cap.ord.df %>% filter(!SampleID=="COMPASS.Dec2021.016")##filter out the outlier for wle COMPASS.Dec2021.016

cap_plot<- ggplot(data=cap.ord.df, aes(x=CAP1, y=CAP2)) +
  geom_point(aes(fill=transect, shape=site_new), size=6, color="black")+
  scale_fill_manual(values=c("#CC79A7","#E69F00","#0072B2", "#CC6677"))+
  scale_shape_manual(values=c(21, 22,24))+ ##change site names for CB
  ggtitle("Chesapeake")+
  labs(x=paste("CAP1 (", fracvar[1], "%)", sep=""),y=paste("CAP2 (", fracvar[2], "%)", sep=""))+
  labs(shape="site", fill="transect")+
  theme_classic()+
  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14))+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))+
  theme(legend.text = element_text(size=12),  legend.title = element_text(size=16))+
  labs(fill="transect", shape="site")+
  guides(fill=guide_legend(override.aes = list(shape=21), order=1),
         shape = guide_legend(override.aes = list(fill = "black"), order=2))+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16 ))


#cap_plot$labels$site.x <- "site"
#cap_plot$labels$transect.x <- "transect"

cap_plot


# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

## Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
  #mutate(labels = gsub("\\.", ":", labels))
colnames(arrowdf) = c("labels", "xend", "yend")
arrowdf = arrowdf %>%
  left_join(var_goodnames, by = "labels") %>%
  rename(old_labels = labels) %>%
  rename(labels = goodnames)

## Define the arrow aesthetic mapping
arrow_map <- aes(xend = xend, yend = yend, x = 0, y = 0, 
                 color = NULL)

label_map <- aes(x = xend + 0.02*xend/abs(xend), y = yend, ##0.02 IN CB, 0.05 IN WLE
                 color = NULL, label = labels)

arrowhead = arrow(length = unit(0.02, "npc"), type = "closed")

library(ggrepel)
# Make a new graphic including labels and axes for variables.
cap.plot = cap_plot + 
  geom_segment(mapping = arrow_map, linewidth = 1.2, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_segment(mapping = arrow_map, linewidth = 0.5, data = arrowdf, color = "white", arrow = arrowhead) + 
  geom_label_repel(mapping = label_map, data = filter(arrowdf, xend < 0), show.legend = FALSE, size=7*6/14, hjust=0.5, vjust=1, fill="white", color="black") + #HJUST = 1 IN CB, 0.5 in wle
  geom_label_repel(mapping = label_map, data = filter(arrowdf, xend > 0), show.legend = FALSE, size=7*6/14, hjust=0.6, fill="white", color="black") + #HJUST =0 IN CB, 0.6 in wle
  #geom_label_repel()+
  labs(title = paste("Variables explain ", round(100*RsquareAdj(cap_ord)$r.squared, 3), "% of the variation", sep=""), subtitle="in Chesapeake")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16 ),
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=16 ))
# theme_classic() + labs(shape="site", color="transect")+
#theme(legend.text=element_text(size=11),
#legend.title = element_text(size = 12),
#  axis.title.x = element_text(size=12),
#  axis.title.y=element_text(size=12),
# axis.text.x = element_text(size=11),
#  axis.text.y=element_text(size=11))
# theme(legend.position="bottom",
# legend.direction = "vertical") +
#guides(shape = guide_legend(order = 1, nrow=1),
#   fill = guide_legend(order = 2, override.aes=list(shape=22), nrow=1))

cap.plot


ggsave(filename = "cap_cb_bray_fticr.TIFF", plot = cap.plot, ##Final Figure
       width = 17,
       height = 16, units = c("cm"),
       dpi = 300)
