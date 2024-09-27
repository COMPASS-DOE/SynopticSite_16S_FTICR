### Thinking about metabolite tree generation

#options(digits = 10)

#filter.cal = F # Remove poorly calibrated samples (i.e., low QC)

require(phangorn) # For tree based functions
require(ggtree) # For tree visualization (on bioconductor)
require(vegan) # For vegdist


# ################## #
#### Load in data ####
# ################## #

# Set sample name
Sample_Name = "COMPASS_FTICR"

# Load in data
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/FTICR")
#mol = read.csv("icr_meta_100.csv", row.names=1) # Load in molecular data
mol <- read.table("icr_meta_final.csv", sep=",", header = TRUE, row.names = 1, dec = ",")

# Removing peaks that have no formula assignments
#mol = mol[!is.na(mol$formula),]

# Setting objects for useful parameters
Mol.Info = mol[,c("C", "H", "O", "N", "S", "P", "DBE"), drop = F]
Mol.Ratio = mol[,c("HC", "OC")]


# ##################### #
#### Generating tree ####
# ##################### #

# Pairwise distance between peaks


#Mol.Info = as.numeric(apply(Mol.Info, 2, scale), row.names = row.names(Mol.Info))
library(tidyverse)

Mol.Info = apply(Mol.Info, 2, scale)
#Mol.Info = as.numeric(Mol.Info)
Mol.Info = as.data.frame(cbind(row.names(mol), Mol.Info)) # Generating a distance matrix based upon the provided parameters
row.names(Mol.Info)<- Mol.Info[,1]
#Mol.Info<- Mol.Info %>% select(!V1)
Mol.Info<- Mol.Info[,-1]
Mol.Info
class(Mol.Info)
Mol.Info<- as_tibble(Mol.Info)
Mol.Info
#Mol.Info<- as.matrix(Mol.Info)

Mol.Info.numeric <- apply(Mol.Info, 2, as.numeric)
rownames(Mol.Info.numeric) <- rownames(mol)



library(usethis) 
usethis::edit_r_environ()

# Create tree
library(ape)
library(vegan)
library(phangorn)
library(ggtree)
library(fastcluster)
#rm(mol)
#tree = as.phylo(hclust(vegdist(Mol.Info, "euclidean"), "average")) # Converting the distance matrix into a tree
tree= vegan::vegdist(Mol.Info.numeric, "euclidean")
tree1= hclust(tree, "average")
tree2= as.phylo(tree1)
# Quick visualization (with Van Krevelen compound classes)
mol = cbind(row.names(mol), mol) # Peak names need to be the first column for ggtree to work

col = colorRampPalette(c("dodgerblue4", "goldenrod3", "firebrick3", "orange"))(length(unique(mol$Class)))
tree_plot<- ggtree(tree2, layout = "circular") %<+% mol + geom_tippoint(aes(color = Class), na.rm = T) +
  scale_color_manual(values = col)+
  theme(legend.position = "right")


ggsave("tree_final_100features_test.TIFF",tree_plot, width=8, height=8)

# Writing tre
write.tree(tree2, paste(Sample_Name, "_MCD_UPGMA.tre", sep = ""))

##generating tree for 16S data, but this is sample by sample, need this as asv to asv because this is feature NTI? same as feature to feature comparison in FTICR?

library(phyloseq)

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

library("ape")
set.seed(1)
random_tree = rtree(ntaxa(phyloseq_merged), rooted=TRUE, tip.label=taxa_names(phyloseq_merged))
plot(random_tree)
phyloseq_tree = merge_phyloseq(phyloseq_merged, random_tree)
phyloseq_tree
plot_tree(phyloseq_tree)

phyloseq_wunifrac <- phyloseq::distance(phyloseq_tree, method = "wunifrac" )



