
#############################################
#######Generating trees #####################
#############################################

################################################################################
#####FTCIR final tree was generated on the cluster, but also worked locally below
##################################################################################


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
mol <- read.table("icr_meta_final.csv", sep=",", header = TRUE, row.names = 1, dec = ",")# Load in molecular data

# Removing peaks that have no formula assignments
#mol = mol[!is.na(mol$formula),]

# Setting objects for useful parameters
Mol.Info = mol[,c("C", "H", "O", "N", "S", "P", "DBE"), drop = F]
Mol.Ratio = mol[,c("HC", "OC")]


# ##################### #
#### Generating tree ####
# ##################### #

# Pairwise distance between peaks

library(tidyverse)

Mol.Info = apply(Mol.Info, 2, scale)

Mol.Info = as.data.frame(cbind(row.names(mol), Mol.Info)) # Generating a distance matrix based upon the provided parameters
row.names(Mol.Info)<- Mol.Info[,1]

Mol.Info<- Mol.Info[,-1]
Mol.Info
class(Mol.Info)
Mol.Info<- as_tibble(Mol.Info)
Mol.Info


Mol.Info.numeric <- apply(Mol.Info, 2, as.numeric)
rownames(Mol.Info.numeric) <- rownames(mol)

# Create tree
library(ape)
library(vegan)
library(phangorn)
library(ggtree)
library(fastcluster)

# Converting the distance matrix into a tree
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

####################################################################
##generating tree for 16S data, but this is sample by sample, need this as asv to asv because this is feature NTI? same as feature to feature comparison in FTICR?
##tree for 16S was generated on the cluster, did not use random tree below 
####################################################################

######################################################################
###code below not needed,16S tree generated from sequence file in the HPC
#####################################################################
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

###############################################################################
##### BNTI create nulls - submitted via jobs to the HPC - for both FTICR and 16S
##############################################################################

################################################################################################################
#### bnti merge nulls for FTICR and 16S - code is saved as separate R files within "merge nulls - bnti" folders
##############################################################################################################

#####################
###### WGCNA ##########
#####################

### Analyzing single sample bNTI feature

# Load in libraries
library(ggplot2); library(reshape2); library(gplots); library(ggpubr); library(ggtree)
library(plyr); library(dplyr); library(stringr)
library(vegan)
library(picante)
library(Hmisc)
library(igraph); library(qgraph); library(WGCNA)

# Focus sample
#focus = "ECA_0Cyc_R2"

# Flags
plot.results = F
spear.cor = F
WGCNA = T


# ############################# #
#### Load and clean the data ####
# ############################# #

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Code from first meeting with Bob")

# Sample-resolved bNTI feature
dna.bnti = read.csv("bNTI_Feat_16s_withoutConsp_bNTI_OTU_by_samp_10rep.csv", row.names = 1)
#rna.bnti = read.csv("", row.names = 1)
icr.bnti = read.csv("bNTI_Feat_ICR_withoutConsp_CB_bNTI_feature_by_samp_10rep.csv", row.names = 1)

##########################################
##loading in bnti files for cb upland only
##########################################

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - 16S - CB/upland CB/100 perm")
dna.bnti=read.csv("bNTI_16S_OTU_Synoptic_CB_upland_allsites_16S_100_bNTI_feature_by_samp_transect-upland_100rep.csv", row.names=1)#change file names as needed

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - FTICR - CB/upland CB/100 perm")
icr.bnti<- read.csv("bNTI_Feat_OTU_Synoptic_CB_upland_allsites_FTICR_100_bNTI_feature_by_samp_transect-upland_100rep.csv", row.names=1)#change file names as needed

# Load in trees/dendrogram
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - 16S - CB/upland CB/100 perm")
asv.tree = read.tree("tree_cb_full_rooted_16s.nwk")
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - FTICR - CB/upland CB/100 perm")
icr.tree = read.tree("COMPASS_FTICR_Synoptic_Final.tre")


##########################################
##loading in bnti files for wle upland only
##########################################

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - 16S WLE/upland WLE/100 perm")
dna.bnti=read.csv("bNTI_16S_OTU_Synoptic_WLE_upland_allsites_16S_100_bNTI_feature_by_samp_transect-upland_100rep.csv", row.names=1)#change file names as needed

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - FTICR WLE/upland WLE/100 perm")
icr.bnti<- read.csv("bNTI_Feat_OTU_Synoptic_WLE_upland_allsites_FTICR_100_bNTI_feature_by_samp_transect-upland_100rep.csv", row.names=1)#change file names as needed

# Load in trees/dendrogram
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - 16S WLE/upland WLE/100 perm")
asv.tree = read.tree("tree_wle_full_rooted_16s.nwk")
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - FTICR WLE/upland WLE/100 perm")
icr.tree = read.tree("COMPASS_FTICR_Synoptic_Final.tre")



##########################################
##loading in bnti files for wle transition only
##########################################

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - 16S WLE/transition WLE")
dna.bnti=read.csv("bNTI_16S_OTU_Synoptic_WLE_transition_allsites_16S_100_bNTI_feature_by_samp_transect-transition_100rep.csv", row.names=1)#change file names as needed
asv.tree = read.tree("tree_wle_full_rooted_16s.nwk")

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - FTICR WLE/transition WLE")
icr.bnti<- read.csv("bNTI_Feat_OTU_Synoptic_WLE_transition_allsites_FTICR_100_bNTI_feature_by_samp_transect-transition_100rep.csv", row.names=1)#change file names as needed
icr.tree = read.tree("COMPASS_FTICR_Synoptic_Final.tre")

##########################################
##loading in bnti files for wle wetland only
##########################################

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - 16S WLE/wetland WLE")
dna.bnti=read.csv("bNTI_16S_OTU_Synoptic_WLE_wetland_allsites_16S_100_bNTI_feature_by_samp_transect-wetland_100rep.csv", row.names=1)#change file names as needed
asv.tree = read.tree("tree_wle_full_rooted_16s.nwk")

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - FTICR WLE/wetland WLE")
icr.bnti<- read.csv("bNTI_Feat_OTU_Synoptic_WLE_wetland_allsites_FTICR_100_bNTI_feature_by_samp_transect-wetland_100rep.csv", row.names=1)#change file names as needed
icr.tree = read.tree("COMPASS_FTICR_Synoptic_Final.tre")

##########################################
##loading in bnti files for cb transition only
##########################################

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - 16S - CB/transition CB")
dna.bnti=read.csv("bNTI_16S_OTU_Synoptic_CB_transition_allsites_16S_100_bNTI_feature_by_samp_transect-transition_100rep.csv", row.names=1)#change file names as needed
asv.tree = read.tree("tree_cb_full_rooted_16s.nwk")

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - FTICR - CB/transition CB")
icr.bnti<- read.csv("bNTI_Feat_OTU_Synoptic_CB_transition_allsites_FTICR_100_bNTI_feature_by_samp_transect-transition_100rep.csv", row.names=1)#change file names as needed
icr.tree = read.tree("COMPASS_FTICR_Synoptic_Final.tre")

##########################################
##loading in bnti files for cb wetland only
##########################################

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - 16S - CB/wetland CB")
dna.bnti=read.csv("bNTI_16S_OTU_Synoptic_CB_wetland_allsites_16S_100_bNTI_feature_by_samp_transect-wetland_100rep.csv", row.names=1)#change file names as needed
asv.tree = read.tree("tree_cb_full_rooted_16s.nwk")

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - FTICR - CB/wetland CB")
icr.bnti<- read.csv("bNTI_Feat_OTU_Synoptic_CB_wetland_allsites_FTICR_100_bNTI_feature_by_samp_transect-wetland_100rep.csv", row.names=1)#change file names as needed
icr.tree = read.tree("COMPASS_FTICR_Synoptic_Final.tre")


# Taxonomic information
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting")
data = read.table("taxonomy-dn-99.csv", sep = ",", row.names = 1, header = T)
row.tax = row.names(data)
#tax = data[,"taxonomy"]

#tax = strsplit(tax, "; ")
#tax = ldply(tax, rbind)
#colnames(tax) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

for(i in 2:ncol(data)){
  data[is.na(data[,i]), i] = data[is.na(data[,i]), (i-1)]
}

row.names(data) = row.tax

rm(i,row.tax)

# Molecular information

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Generating trees")
mol = read.csv("icr_meta_final.csv")



# ########### #
#### WGCNA ####
# ########### #
dna.bnti<- dna.bnti %>% select(!s20230217034) ##remove as needed if working with full CB dataset
icr.bnti.final <- icr.bnti[, order((colnames(icr.bnti)))]


dna.bnti.transpose<- t(dna.bnti)
icr.bnti.transpose<- t(icr.bnti.final)


###match the two data frames to make sure that the same sample names are present in both
##just ordering is not enough. if different samples are present, they will still be ordered but the final merge of sample names may be messed up 

##when working with cb upland only (remove samples as needed)
icr.bnti <- icr.bnti[, order((colnames(icr.bnti)))]
dna.bnti <- dna.bnti[, order((colnames(dna.bnti)))]
icr.bnti.final<- icr.bnti %>% select(!s20230217078)
dna.bnti.transpose<- t(dna.bnti)
icr.bnti.transpose<- t(icr.bnti.final)

##when working with wle upland only (remove samples as needed)
icr.bnti <- icr.bnti[, order((colnames(icr.bnti)))]
dna.bnti <- dna.bnti[, order((colnames(dna.bnti)))]
dna.bnti.transpose<- t(dna.bnti)
icr.bnti.transpose<- t(icr.bnti)

##when working with wle transition only (remove samples as needed)
icr.bnti <- icr.bnti[, order((colnames(icr.bnti)))]
dna.bnti <- dna.bnti[, order((colnames(dna.bnti)))]
icr.bnti.final<- icr.bnti %>% select(!c(s20230901075,s20230901078))
dna.bnti.transpose<- t(dna.bnti)
icr.bnti.transpose<- t(icr.bnti.final)

##when working with wle wetland only (remove samples as needed)
icr.bnti <- icr.bnti[, order((colnames(icr.bnti)))]
dna.bnti <- dna.bnti[, order((colnames(dna.bnti)))]
icr.bnti.final<- icr.bnti %>% select(!c(s20230901081,s20230901082,s20230901083))
dna.bnti.transpose<- t(dna.bnti)
icr.bnti.transpose<- t(icr.bnti.final)

##when working with cb transition only (remove samples as needed)
##blockwise modules failed with error 
##Error in goodGenes(datExpr, weights, goodSamples, goodGenes, minFraction = minFraction,  : 
##Too few genes with valid expression levels in the required number of samples.
icr.bnti <- icr.bnti[, order((colnames(icr.bnti)))]
dna.bnti <- dna.bnti[, order((colnames(dna.bnti)))]
dna.bnti.final<- dna.bnti %>% select(!c(s20230217034,s20230217117))
dna.bnti.transpose<- t(dna.bnti.final)
icr.bnti.transpose<- t(icr.bnti)

##when working with cb wetland only (remove samples as needed)
icr.bnti <- icr.bnti[, order((colnames(icr.bnti)))]
dna.bnti <- dna.bnti[, order((colnames(dna.bnti)))]
dna.bnti.transpose<- t(dna.bnti)
icr.bnti.transpose<- t(icr.bnti)

##R square cut offs for each transect

##cb upland 
##R square cut off used == 0.8 corresponding to power 6

##wle upland 

##R square cut off used == 0.8 corresponding to power 6

##wle transition

##R square cut off used == 0.696 corresponding power 4

##wle wetland 
##R square cut off used == 0.8 corresponding to power 1

##cb transition

##cb wetland
##R square cut off used: 0.8 corresponding to power 6

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting")

if(WGCNA){ # Switch controlling whether WCGNA will be run or not
  input = cbind(dna.bnti.transpose, icr.bnti.transpose) #bnti.icr.wide.final)
  thresh = 0 ##changing from 0.4
  sign = "signed" #changing unsigned to signed to see if it improves the module colors/assignments
  allowWGCNAThreads()
  #powers = c(1:25) # Generates a series of powers to test for future soft thresholding; soft thresholding allows for continuous data to be described
  sft = pickSoftThreshold(input, powerVector = c(1:20, seq(22, 30, by = 2)),verbose = 5, networkType = "signed", RsquaredCut = 0.8)
  #sft = pickSoftThreshold(input, corOptions = list(use = 'p'), powerVector = powers, verbose = 5) # Tests the powers against the data to see which power fits the best
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2]) # Plots the test to see if they worked (where it flattens is where it 'worked')
  
  p = -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
  p = sft$fitIndices[,1][which(p %in% max(-sign(sft$fitIndices[,3])*sft$fitIndices[,2]))] # Automatically generating the most 'optimal' power
  p = sft$powerEstimate ##always adjusting the R squared cut so that this will give the power value corresponding to the max R square
  block= blockwiseModules(input, power = p, networkType = "signed",TOMType = sign, minModuleSize = 5,
                          numericLabels = TRUE)
  #block = blockwiseModules(input, power = p, TOMType = sign, minModuleSize = 5,
                           #numericLabels = TRUE) # Assigns the data to modules by transforming the data with the calculated power
  
  moduleLabels = block$colors # Saving labels from the modules
  moduleColors = labels2colors(block$colors) # Saves colors of the modules
  MEs = block$MEs # Saving module eigenvalues which contain the information used to designate modules
  taxTree = block$dendrograms[[1]] # Saving the dendrograms which were generated from TOM data created during module assignment
  
  plotDendroAndColors(taxTree, moduleColors[block$blockGenes[[1]]], "Module colors",
                      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05) # Plots modules and dendrograms
  
  plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap", plotDendrograms = FALSE, xLabelsAngle = 90) # Plots relationships between eigenvectors to see what modules are most closely related
  
  plot(hclust(as.dist(1-cor(MEs, use="pairwise.complete.obs")), method = "average")) # Just a plot which can help to see if any module merging can occur
  
  # TOM calculations
  TOM = TOMsimilarity(adjacency(input, power = p, type = sign), useInternalMatrixAlgebra = TRUE) # Calculates the "Topographcial Overlap Matrix"; a measure of how 'related' two nodes are, what their connections are, how adjacent, etc.
  dissTOM = 1 - TOM # Changes similarity matrix to dissimilarity matrix
  plotTOM = dissTOM^7 # Raises the dis. matrix to a power in order to bring out the effects of 'moderate' effects
  diag(plotTOM) = NA # Sets the diagonal of the data to NA's so that there isn't any weird scaling (they are normally one due to self-comparisons)
  
  # # The TOM plot isn't really necessary
  # TOMplot(plotTOM, taxTree, moduleColors) # Generates the plot which looks at relationships between each constituent
  
  # Meta-data for exporting
  w = apply(input, 2, median)
  w = cbind(w, moduleColors)
  
  dimnames(TOM) = list(colnames(input), colnames(input)) # Sets the dimension names to over this matrix to the taxa it was derived from
  quantile(TOM)
  hist(TOM)
  cyt = exportNetworkToCytoscape(TOM, edgeFile = paste0("Cytoscape Files signed WLE transition 100 perm/TOM threshold zero/DNA-ICR_Edges_", thresh, "thresh_", sign, ".txt"), 
                                 nodeFile = paste0("Cytoscape Files signed WLE transition 100 perm/TOM threshold zero/DNA-ICR_Nodes_", thresh, "thresh_", sign, ".txt"), 
                                 threshold = thresh,
                                 nodeAttr = moduleColors) #w) # Threshold is the minimum allowable value for in the TOM matrix
  
  # Creating a module list
  colors = data.frame(Members = colnames(input), colors = moduleColors, hexColor = col2hex(moduleColors), type = "Microbes")
  colors$type[which(colors$Members %in% row.names(icr.bnti))] = "Metabolites"
  color.mol = mol[,c("formula","Class_detailed", "NOSC", "AImod", "HC", "OC")]; colnames(color.mol)[1] = "Members"
  color.tax = data.frame(Members = row.names(data), Order = data$Order)
  colors = colors %>% left_join(color.tax, by = "Members") %>% left_join(color.mol, by = "Members")
  write.csv(colors, paste0("Cytoscape Files signed WLE transition 100 perm/TOM threshold zero/DNA-ICR_Modules_", thresh, "thresh_", sign, ".csv"), row.names = F, quote = F)
}


# pulling out modules for correlation calcs (module eigenvalue correlation with feature bnti)

data<- data %>% mutate(OTUID=row.names(data))
tax.new<- data %>% mutate(query=OTUID)
icr.meta.new<- mol %>% mutate(query=formula)
ME.trans = block$MEs
input_t<- t(input)

# compiling eigengene information
ME.trans.meta = data.frame(query = names(block$colors),
                           Module = gsub("^", "ME", block$colors)) %>%
  left_join(tax.new, by = "query") %>% left_join(icr.meta.new, by="query")

# identifying strongest module membership
ME.strength = NULL

for(curr.mod in colnames(ME.trans)){
  # select matching data
  temp.meta = ME.trans.meta %>%
    filter(Module %in% curr.mod) %>%
    filter(!duplicated(query))
  
  # temporary genes from modules
  temp.data = input_t[which(row.names(input_t) %in% 
                                    temp.meta$query),]
  
  # correlate genes to module
  temp.corr = cor(t(temp.data),  
                  ME.trans[,which(colnames(ME.trans) %in% curr.mod)], use='pairwise.complete.obs')
  
  # convert to data frame
  temp.corr = data.frame(query = row.names(temp.corr), r_value = temp.corr)
  
  # merge into empty object
  ME.strength = rbind(ME.strength, temp.corr)
  
  # clean up
  #rm(temp.meta, temp.data, temp.corr)
  
}

# add in module correlations to the ME.meta
ME.trans.meta = ME.trans.meta %>% left_join(ME.strength)

# clean up
rm(block, sft, p, curr.mod, ME.strength)

write.csv( ME.trans.meta, "correlation_feature_module_WLE_transition.csv")

# Evaluating significant module characteristics

### r wgcna-module-analysis
# selecting the significant modules or all modules
all.mod = ME.trans.meta #%>%
 # filter(Module %in% c("ME28"))

# network stats on modules
network.stats = NULL # empty objects

input = cbind(dna.bnti.transpose, icr.bnti.transpose) #bnti.icr.wide.final)
input_t<- t(input)

for(curr.mod in unique(all.mod$Module)){
  # specify input
  input = t(input_t[which(row.names(input_t) %in% 
                                  all.mod$query[which(all.mod$Module %in%curr.mod)]),])
  
  # calculate TOM
  TOM = TOMsimilarity(adjacency(input, power = 6, type = "signed"), ##change power as per transect/location
                      TOMType = "signed")
  
  # set dimnames
  dimnames(TOM) = list(colnames(input), colnames(input))
  
  # export subnetwork
 exportNetworkToCytoscape(TOM,
                          edgeFile = paste0("module specific nodes and edges zero/Cytoscape_Edges_",
                                            curr.mod, "_0-thresh_signed",
                                             ".txt"),
                           nodeFile = paste0("module specific nodes and edges zero/Cytoscape_Nodes_",
                                             curr.mod, "_0-thresh_signed",
                                             ".txt"),
                           threshold = 0)
    
  # setting up adjacency for igraph
  TOM[TOM > 0.1] = 1 # same threshold for export
  TOM[!TOM == 1] = 0
  
  # network
  network = TOM %>% 
    graph_from_adjacency_matrix() %>%
    simplify()
  network = delete_vertices(network, v = degree(network)==0)
  
  # calculating statistics
  network = data.frame(query = names(degree(network)),
                       Degree = degree(network),
                       Closeness = closeness(network),
                       Eigen = eigen_centrality(network)$vector,
                       Betweenness = betweenness(network),
                       Hub_Score = hub_score(network)$vector)
  
  # add to empty object
  network.stats = rbind(network.stats, network)
  
  # clean up
  #rm(input, TOM, network)
}

# add network stats to object
all.mod = all.mod %>%
  left_join(network.stats)

# clean up
rm(network.stats)

write.csv(all.mod, "network_stats_correlation_by_module_wle_transition.csv")

#did not use code for hub gene below
# plotting hub genes from each module
sig.metaT %>% 
  filter(!duplicated(query)) %>%
  filter(Hub_Score >= 0.9) %>%
  group_by(Module, user_genome) %>%
  summarise(Count = n()) %>%
  left_join(tax, by = "user_genome") %>%
  mutate(Phylum = gsub("p__", "", Phylum)) %>%
  ggplot(aes(x = user_genome, y = Count))+
  geom_bar(stat = "identity", aes(fill = Phylum, alpha = Module))+
  geom_hline(yintercept = 10, color = "red", lty = 2)+
  xlab(NULL)+
  facet_grid(~Phylum, scales = "free_x", space = "free_x")+
  scale_fill_viridis_d(option="plasma")+
  scale_alpha_manual(values = c(0.4, 0.6, 0.8, 1))+
  guides(fill=guide_legend(ncol=2))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())
ggsave("A2-3_MAG_Module_Hub_Genes.pdf", height = 6, width = 10)

#did not use code for metabolite module below
# visualizing composition of metabolite module 
sig.metab %>% 
  mutate(Superclass = case_when(is.na(Superclass) ~ 
                                  "Lipids and lipid-like molecules",
                                .default = Superclass)) %>%
  group_by(Superclass) %>%
  summarise(Count = n()) %>%
  ggplot(aes(x = "", y = Count, fill = Superclass))+
  geom_bar(stat = "identity", color = "white")+
  coord_polar("y", start = 0)+
  scale_fill_viridis_d()+
  theme_void()
ggsave("A2-3_Metabolite_ME2_Composition.pdf")

######################
library(network)

nodes<- read.csv("DNA-ICR_Nodes_0.4thresh_signed.txt", sep="")
edges<- read.csv("DNA-ICR_Edges_0.4thresh_signed.txt", sep="")
write.csv(nodes, "nodes_signed.csv")
write.csv(edges, "edges_signed.csv")

network_basic <- network(edges,
                         vertex.attr = nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)

#####################################################################################################
#####code below was used to match dna and icr sample names when 16S CB and ICR CB+WLE was being used with icr having different sample names
#######################################################################################################

##matching the bnti dna and bnti icr sample names
library(tidyverse)
icr.metadata<- read.csv("icr_long_all_samples.csv")
icr.metadata.unique<- icr.metadata %>% select(!c(formula, n))
icr.metadata.unique<- unique(icr.metadata.unique)
icr.bnti.formula<- icr.bnti %>% mutate(formula=row.names(icr.bnti))
icr.bnti.long<- pivot_longer(icr.bnti.formula, !formula, names_to = "analysis_ID", values_to="value")
icr.bnti.merge <- left_join(icr.bnti.long, icr.metadata.unique, by="analysis_ID")
icr.bnti.merge.cb<- icr.bnti.merge %>% filter(region=="CB")

otu_metadata<- read.csv("sample-metadata-cb-wle.csv")
otu_metadata<- otu_metadata %>% mutate(sample_label=Sample.description)
otu_fticr_merge_names<- left_join(icr.bnti.merge.cb, otu_metadata, by="sample_label")
otu_fticr_merge_names<- otu_fticr_merge_names %>% filter(!sample.id=="NA")
bnti.icr.wide.final<- pivot_wider(otu_fticr_merge_names, id_cols=formula, names_from="sample.id", values_from="value")  
row.names(bnti.icr.wide.final) <- bnti.icr.wide.final$formula
class(bnti.icr.wide.final)
row.names<- row.names(bnti.icr.wide.final)
bnti.icr.wide.final<- bnti.icr.wide.final %>% select(!c("formula")) 
row.names(bnti.icr.wide.final) <- row.names

##removing one sample from bnti.dna because the sample is missing in bnti.icr
##s20230217034

dna.bnti<- dna.bnti %>% select(!s20230217034)
icr.bnti.final <- icr.bnti[, order((colnames(icr.bnti)))]
##arranging column order in both to match them 

bnti.icr.wide.final<- bnti.icr.wide.final %>% select(names(dna.bnti))
test <- match(colnames(dna.bnti), colnames(bnti.icr.wide.final))
##trying to understand the order of samples in icr dataframe
bnti.icr.wide.final <- bnti.icr.wide.final[, order((colnames(bnti.icr.wide.final)))]
row.names(bnti.icr.wide.final) <- row.names

##removing same number 117 from bnti.dna matrix
dna.bnti<- dna.bnti %>% select(!s20230217117)

bnti.icr.wide.final<- bnti.icr.wide.final %>% select (!s20230217078)
row.names(bnti.icr.wide.final) <- row.names

###############################################################
####subsetting icr_wide1.csv to include only CB samples
##keeping the fticr tree the same as before including all known peaks 
#################################################################


icr_wide<- read.csv("icr_wide1.csv", row.names=1)
icr_wide_transpose<- t(icr_wide)
row.names<- row.names(icr_wide_transpose)
icr_wide_transpose<- as.data.frame(icr_wide_transpose)
icr_wide_transpose<- icr_wide_transpose %>% mutate(formula=row.names(icr_wide_transpose))
row.names(icr_wide_transpose)<- NULL
icr_long<- icr_wide_transpose %>% pivot_longer(!formula, names_to = "analysis_ID", values_to = "PA")
class(icr_long)


icr.metadata<- read.csv("icr_long_all_samples.csv")
icr.metadata.unique<- icr.metadata %>% select(!c(formula, n))
icr.metadata.unique<- unique(icr.metadata.unique)
icr.merge <- left_join(icr_long, icr.metadata.unique, by="analysis_ID")
icr.merge.cb<- icr.merge %>% filter(region=="CB")

otu_metadata<- read.csv("sample-metadata-cb-wle.csv")
otu_metadata<- otu_metadata %>% mutate(sample_label=Sample.description)
otu_fticr_merge_names<- left_join(icr.merge.cb, otu_metadata, by="sample_label")
otu_fticr_merge_names<- otu_fticr_merge_names %>% filter(!sample.id=="NA")
icr.wide.final<- pivot_wider(otu_fticr_merge_names, id_cols=formula, names_from="sample.id", values_from="PA")  
row.names(icr.wide.final) <- icr.wide.final$formula
class(icr.wide.final)
row.names<- row.names(icr.wide.final)
icr.wide.final<- icr.wide.final %>% select(!c("formula")) 
row.names(icr.wide.final) <- row.names
icr.wide.final.mat<- as.matrix(icr.wide.final)
icr.wide.final.mat.unlist<- unlist(icr.wide.final.mat)
#icr.wide.final.final<- matrix(as.numeric(icr.wide.final.mat.unlist), ncol=116)


icr.wide.final.mat.unlist %>% data.frame() %>% mutate(across(where(is.character), as.numeric)) %>% as.matrix() -> mat.num
print(which(is.na(mat.num)))
mat.num.noNA<- ifelse(is.na(mat.num), 0, mat.num)
print(which(is.na(mat.num.noNA)))

icr.wide.final1 = mat.num.noNA[rowSums(mat.num.noNA[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples
View(icr.wide.final1) 

##TO CREATE NULLS WE NEED TO TRANSPOSE THE DATA

icr_wide_final_t<- t(icr.wide.final1)
write.csv( icr_wide_final_t, "icr_wide_cb.csv")





