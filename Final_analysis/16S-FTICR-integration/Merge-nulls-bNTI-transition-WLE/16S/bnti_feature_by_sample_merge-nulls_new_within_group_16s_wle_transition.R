Sample_Name = "bNTI_16S_OTU_Synoptic_WLE_transition_allsites_16S_100" # Input sample name
Factor = "transect"
Level = "transition"

# Switches for script behaviors
rm.conspec = F # Remove conspecifics
abund.weig = T # Weight values by relative abundances
#type = F # This has limited functionality specific to my project, can be safely ignored
rm.tax = T # This configures whether the assemblage data has taxonomy

### Load in necessary libraries
require(Rfast) # For faster variant of finding column minimum
require(dplyr) # For inner joining (faster than merge)
require(abind) # For list->array function
require(picante) # For match.phylo.data
require(phytools) # For midpoint.root
require(ggtree); require(ggplot2); require(ggstance); require(reshape2) # For plotting the results

# Combine matrices as array
acomb = function(...) abind(..., along = 3)


# ################## #
#### Load in data ####
# ################## #

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cleaning and organizing /Create nulls - 16S WLE/transition WLE")
data<- read.csv("OTU_WLE_rarefied20k.csv", row.names=1, header=TRUE)
#data<- t(data)
tree = read.tree("tree_wle_full_rooted_16s.nwk")
meta<- read.csv("sample_metadata_wle.csv")

# #################### #
#### Pre-processing ####
# #################### #

### Subsetting data based upon the group
# Ensure meta and data are in the same order
data = data[,which(colnames(data) %in% meta$sampleid)]
meta = meta[which(meta$sampleid %in% colnames(data)),]

if(!identical(meta$sampleid, colnames(data))){
  print("Your metadata/factor sheet doesn't match the provided data, attepmting to fix it")
  data = data[,meta$sampleid]
  
  if(!identical(meta$sampleid, colnames(data))){
    stop("Your sample names couldn't be easily fixed, so I'm stopping to prevent damage")
  }
  
  print("Names were fixed")
}

# Selecting sample subset
fac.col = which(colnames(meta) %in% Factor)
fac.samp = meta$sampleid[which(meta[,fac.col] %in% Level)]
data = data[,which(colnames(data) %in% fac.samp)]

if(min(rowSums(data)) == 0){
  data = data[-which(rowSums(data) == 0),]
}

rm(fac.col, fac.samp, meta)

data<- t(data)
phylo = match.phylo.comm(tree, data) # Matching ICR dataset to the tree

data = phylo$comm ##otus reduce from 5171 (or something) to 4837, but why?
tree = phylo$phy

rm("phylo")

# Storing dimensions
samp.num = dim(data)[1]
mem.num = dim(data)[2]

# ####################### #
#### Merging null reps ####
# ####################### #

# Merging the separate bMNTD files
files = list.files(path = "transect-transition/",
                   pattern = "Feature_Samp_Null", full.names = T) # Listing files

null.by.samp = NULL # Dummy object

for(curr.file in files){
  temp = as.data.frame(read.csv(curr.file, row.names = 1))
  null.by.samp = c(null.by.samp, list(temp))
} # Merging individual 

null.by.samp = do.call(acomb, list(null.by.samp))
rm("curr.file")

# ###################################### #
#### Determine observed bNTI(feature) ####
# ###################################### #

# Measuring distances
library(iCAMP)
coph = pdist.big(tree, output=TRUE)

# Removing conspecifics, if desired
if(rm.conspec){
  coph[coph == 0] = NA
}

# Creating empty object 
min.neigh = NULL # Creating an empty matrix to store pairwise results
comp.names = NULL

# Running through the pairwise comparisons
for (i in 1:(samp.num - 1)) {
  for (j in (i + 1):samp.num) {
    
    # Selecting members in samples and subsetting coph. correspondingly
    samp1 = colnames(data[i, data[i, ] > 0, drop = FALSE]) # Members in the first sample
    samp2 = colnames(data[j, data[j, ] > 0, drop = FALSE]) # Members in the second sample
    
    pair.dist = coph[samp1, samp2, drop = FALSE]
    
    # First sample minimums
    if(length(which(is.na(pair.dist[,1]))) > 0){
      min.dist1 = apply(pair.dist, 1, min, na.rm = T)
    } else {
      min.dist1 = rowMins(pair.dist, value = T)
    } # There is a bug with Rfast minimum calculations that causes it to report NA if it is the first value
    
    names(min.dist1) = row.names(pair.dist)
    
    # Second sample minimums
    if(length(which(is.na(pair.dist[1,]))) > 0){
      min.dist2 = apply(pair.dist, 2, min, na.rm = T)
    } else {
      min.dist2 = colMins(pair.dist, value = T)
    } # There is a bug with Rfast minimum calculations that causes it to report NA if it is the first value
    
    names(min.dist2) = colnames(pair.dist)
    
    # If desired, weighting distances by relative abundance
    if(abund.weig){
      min.dist1 = min.dist1*(data[i, samp1]/sum(data[i, samp1]))
      min.dist2 = min.dist2*(data[j, samp2]/sum(data[j, samp2]))
    }
    
    min.dist = c(min.dist1, min.dist2) # Combining those minimum distances
    min.dist = data.frame(Names = names(min.dist), Dist = min.dist) # Converting to a data frame for easier aggregation
    
    if(rm.conspec){
      # If there are duplicates when conspecifics have been removed, these are averaged as they will have different values
      min.dist = aggregate(Dist~Names, min.dist, FUN = mean) 
    } else {
      # If conspecifics haven't been removed, duplicates will by necessity be 0 - if a community member 
      # is in both samples, it has to be its own nearest neighbor
      if(length(which(duplicated(min.dist$Names) %in% TRUE)) > 0){
        if(mean(min.dist$Dist[duplicated(min.dist$Names)]) == 0){
          min.dist = min.dist[!duplicated(min.dist$Names),]
        } else {
          # If the duplicates aren't 0, something is wrong in the data - stopping the script
          stop("Something odd is happening with your duplicated values. Check that out.")
        }
      }
    } # Checking to ensure duplicates exist at all
    
    # Creating an object with all commiunity members to merge in those which were present in these two samples
    merge.dist = data.frame(Names = colnames(data))
    
    # Adding observed community member minimum distances to 
    merge.dist = left_join(x = merge.dist, y = min.dist, by = "Names")
    row.names(merge.dist) = merge.dist$Names
    
    # Removing names column
    merge.dist = merge.dist[,-which(colnames(merge.dist) %in% "Names"),drop = F]
    
    # Adding this pairwise comparison to the overall object
    min.neigh = cbind(min.neigh, as.matrix(merge.dist))
    comp.names = c(comp.names, paste0(row.names(data)[i], "-", row.names(data)[j]))
    
    # Clocking
    print(c(i, j, date()))
  }
} # Loop works through the pairwise comparisons

# Clean-up
rm(i, j, samp1, samp2, min.dist, min.dist1, min.dist2, merge.dist)

# Setting column names
colnames(min.neigh) = comp.names

# Finding average minimum neighbor distance by sample
min.by.samp = matrix(data = NA, nrow = nrow(min.neigh), ncol = nrow(data))
row.names(min.by.samp) = row.names(min.neigh)
colnames(min.by.samp) = row.names(data)

for(i in 1:samp.num){
  # Selecting current sample
  curr.samp = grep(paste0("^",row.names(data)[i], "-|", "-", row.names(data)[i], "$"), colnames(min.neigh))
  
  # Selecting members in current sample; min.neigh row.names and column names on data are identical
  curr.mem = which(row.names(min.by.samp) %in% names(data[i,which(data[i,] > 0)]))
  
  # Creating temp object
  temp = min.neigh[,curr.samp]
  
  # Setting all values not for the current sample to NA
  temp[-curr.mem,] = NA
  
  # Adding to final output matrix
  min.by.samp[,i] = rowMeans(temp, na.rm = F)
  
  # Clean-up
  rm(temp, curr.mem, curr.samp)
} 


# ############################# #
#### Calculate bNTI(feature) ####
# ############################# #

# Measuring the difference between observed and null distances
if(!identical(row.names(min.by.samp), row.names(null.by.samp))){
  stop("Your null and observed results do not match. Please double check them.")
}

feat.by.samp = matrix(c(NA), nrow = nrow(null.by.samp), ncol = ncol(null.by.samp))
row.names(feat.by.samp) = row.names(min.by.samp)
colnames(feat.by.samp) = colnames(min.by.samp)

for(i in 1:ncol(null.by.samp)){
  for(j in 1:nrow(null.by.samp)){
    m = null.by.samp[j,i,] # Just setting all the randomizations for a given comparison to a matrix
    feat.by.samp[j,i] = ((min.by.samp[j,i]-mean(m))/sd(m)) # The bNTI calculation
  }
}

write.csv(feat.by.samp, paste0(Sample_Name, "_bNTI_feature_by_samp_", Factor, "-", Level, "_", length(files), "rep.csv"), quote = F)


# ################## #
#### Plot results ####
# ################## #

feat.by.samp<- read.csv("bNTI_16S_OTU_Synoptic_WLE_transition_allsites_16S_100_bNTI_feature_by_samp_transect-transition_100rep.csv", row.names=1)

# Creating data frame to make plotting easier
feat.samp.frame = data.frame(Member = row.names(feat.by.samp), Value = feat.by.samp[,1], Direction = "Insignificant", stringsAsFactors = F)
feat.samp.frame$Direction[which(feat.samp.frame$Value <= -1)] = "Low"
feat.samp.frame$Direction[which(feat.samp.frame$Value >= 1)] = "High"
feat.samp.frame$Direction[which(feat.samp.frame$Value <= -2)] = "Sig. Low"
feat.samp.frame$Direction[which(feat.samp.frame$Value >= 2)] = "Sig. High"

feat.samp.frame$Direction = factor(feat.samp.frame$Direction, levels = c("Sig. High", "High", "Insignificant", "Low", "Sig. Low"))

# Plotting a standard tree with a point graph to see the scale of differences
# Specified labels and their new names
specified_labels <- c("8f0fb68673a8b3c26f3fd210e29b7916", "d829bee4984f82ffc2453212157caf96")
label_replacements <- data.frame(
  original_label = specified_labels,
  display_label = c("Clostridiales", "Rhizobiales")
)

# Prepare tree and visually check layout
p <- ggtree(tree)

# Add labels to specified tips directly using geom_tiplab
p1 <- p + geom_tiplab(aes(label = ifelse(label %in% label_replacements$original_label, 
                                         label_replacements$display_label[match(label, label_replacements$original_label)], 
                                         NA)), size = 3, color = "blue", offset=0.2)

p1

##labeling the point data also 
specified_labels <- c("8f0fb68673a8b3c26f3fd210e29b7916", "d829bee4984f82ffc2453212157caf96")
highlight_labels <- c("Clostridiales", "Rhizobiales")

# Creating a mapping for labeling
label_map <- data.frame(Member = specified_labels, Label = highlight_labels)

# Add label mapping to feat.samp.frame with a left join
feat.samp.frame <- left_join(feat.samp.frame, label_map, by = "Member")

# Use facet_plot to add both points and labels in the same panel
# Add points to the facet
p2 <- facet_plot(p, panel = "Sample Nearest Neighbor", data = feat.samp.frame, 
                 geom = geom_point, mapping = aes(x = Value, y = y, color = Direction), size=3)

# Add labels within the same facet
p3 <- facet_plot(p2, panel = "Sample Nearest Neighbor", data = feat.samp.frame, 
                 geom = geom_text_repel, mapping = aes(x = Value, y = y, label = Label),
                 hjust = -0.5, color = "blue", size=4, na.rm = TRUE, max.overlaps=18)

p3

p4<- p3 + scale_color_manual(values = c("firebrick4", "firebrick1", "goldenrod1", "dodgerblue1", "dodgerblue4"), drop = F) + 
  theme_tree2()+theme(legend.position = "left")+
  theme(axis.text.x = element_text(size=16, face="bold"),
        legend.text = element_text(size = 16),  # Increase legend labels size
        legend.title = element_text(size = 16) ,
        strip.text = element_text(size = 16, face="bold"))


p4

ggsave("bnti_wle_transition_all_sites_16s_100perm_06.05.25_labels_ggrepel.TIFF", plot=p4, dpi=300, height=10, width=20, units=c("in"))

##old plot
# Plotting a standard tree with a point graph to see the scale of differences
p = ggtree(tree)
p = facet_plot(p, panel = "Sample Nearest Neighbor", data = feat.samp.frame, 
               geom = geom_point, mapping = aes(x = Value, y = y, color = Direction))
p + scale_color_manual(values = c("firebrick4", "firebrick1", "goldenrod1", "dodgerblue1", "dodgerblue4"), drop = F) + 
  theme_tree2()+theme(legend.position = "left")

# Calculating percentages
feat.samp.frame = data.frame(Member = row.names(feat.by.samp), Sig_High = NA, High = NA, Insignificant = NA, 
                             Low = NA, Sig_Low = NA, stringsAsFactors = F)
feat.samp.frame$Sig_High = (apply(feat.by.samp, 1, function(x) length(which(x >= 2)))/apply(feat.by.samp, 1, function(x) length(which(!is.na(x)))))*100
feat.samp.frame$High = (apply(feat.by.samp, 1, function(x) length(which(x >= 1 & x < 2)))/apply(feat.by.samp, 1, function(x) length(which(!is.na(x)))))*100
feat.samp.frame$Insignificant = (apply(feat.by.samp, 1, function(x) length(which(abs(x) < 1)))/apply(feat.by.samp, 1, function(x) length(which(!is.na(x)))))*100
feat.samp.frame$Low = (apply(feat.by.samp, 1, function(x) length(which(x <= -1 & x > -2)))/apply(feat.by.samp, 1, function(x) length(which(!is.na(x)))))*100
feat.samp.frame$Sig_Low = (apply(feat.by.samp, 1, function(x) length(which(x <= -2)))/apply(feat.by.samp, 1, function(x) length(which(!is.na(x)))))*100

feat.samp.frame = melt(feat.samp.frame, id.vars = "Member")
feat.samp.frame$variable = factor(feat.samp.frame$variable, 
                                  levels = c("Sig_High", "High", "Insignificant", "Low", "Sig_Low"))


p = ggtree(tree)
p = facet_plot(p, panel = "Percentage", data = feat.samp.frame, geom = geom_barh,
               mapping = aes(x = value, fill = as.factor(variable)), stat = "identity")
# p = facet_plot(p, panel = 'Significantly High', data = feat.samp.frame, geom = geom_segment, 
#                mapping = aes(x = 0, xend = Sig_High, y = y, yend = y), color = "firebrick4")
# p = facet_plot(p, panel = 'High', data = feat.samp.frame, geom = geom_segment, 
#                mapping = aes(x = 0, xend = High, y = y, yend = y), color = "firebrick1")
# p = facet_plot(p, panel = 'Insignificant', data = feat.samp.frame, geom = geom_segment, 
#                mapping = aes(x = 0, xend = Insignificant, y = y, yend = y), color = "goldenrod1")
# p = facet_plot(p, panel = 'Low', data = feat.samp.frame, geom = geom_segment, 
#                mapping = aes(x = 0, xend = Low, y = y, yend = y), color = "dodgerblue1")
# p = facet_plot(p, panel = 'Significantly Low', data = feat.samp.frame, geom = geom_segment, 
#                mapping = aes(x = 0, xend = Sig_Low, y = y, yend = y), color = "dodgerblue4")
p + scale_fill_manual(values = c("firebrick4", "firebrick1", "goldenrod1", "dodgerblue1", "dodgerblue4"), drop = F) + 
  theme_tree2()+theme(legend.position = "left")


# Plotting variable members
feat.order = order(apply(feat.by.samp, 1, mean, na.rm = F), decreasing = T) # Needed to find most variable members
feat.melt = melt(feat.by.samp[feat.order[1:5],]) # Melting most variable members for ggplot

ggplot(data = feat.melt, aes(x = Var2, y = value, group = as.character(Var1)))+
  geom_point(aes(color = as.character(Var1)))+
  geom_line(aes(color = as.character(Var1)))+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(colour = "black", angle = -40, vjust = 1, hjust = 0),
        axis.text.y = element_text(colour = "black"),
        panel.border = element_rect(size = 1, colour = "black"),
        panel.grid = element_blank())


###################
##FEMS Eco Evo edits
###################

# ################## #
#### Plot results ####
# ################## #

feat.by.samp<- read.csv("bNTI_16S_OTU_Synoptic_WLE_transition_allsites_16S_100_bNTI_feature_by_samp_transect-transition_100rep.csv", row.names=1)

# Creating data frame to make plotting easier
#feat.samp.frame = data.frame(Member = row.names(feat.by.samp), Value = feat.by.samp[,1], Direction = "Insignificant", stringsAsFactors = F)
library(tidyr)
library(dplyr)
library(tibble)
library(tidyverse)

# The code to convert from wide to long
feat.samp.frame <- feat.by.samp %>%
  rownames_to_column(var = "Member") %>%   # 1. Turn row names into a real column called "Member"
  pivot_longer(
    cols = -Member,                       # 2. Pivot all columns EXCEPT the "Member" column
    names_to = "Sample",                  # 3. The names of the old columns go into a new column called "Sample"
    values_to = "Value"                   # 4. The values from those columns go into a new column called "Value"
  )

# View the new long-format data frame
print(feat.samp.frame)
View(feat.samp.frame)

feat.samp.frame<- feat.samp.frame %>% mutate(Direction="Insignificant")

##change the direction column
feat.samp.frame$Direction[which(feat.samp.frame$Value <= -1)] = "Low"
feat.samp.frame$Direction[which(feat.samp.frame$Value >= 1)] = "High"
feat.samp.frame$Direction[which(feat.samp.frame$Value <= -2)] = "Sig. Low"
feat.samp.frame$Direction[which(feat.samp.frame$Value >= 2)] = "Sig. High"

feat.samp.frame$Direction = factor(feat.samp.frame$Direction, levels = c("Sig. High", "High", "Insignificant", "Low", "Sig. Low"))

# Plotting a standard tree with a point graph to see the scale of differences
# Prepare tree and visually check layout
p <- ggtree(tree)
p

##labeling the point data also 
specified_labels <- c("8f0fb68673a8b3c26f3fd210e29b7916", "d829bee4984f82ffc2453212157caf96")
highlight_labels <- c("Clostridiales", "Rhizobiales")

# Creating a mapping for labeling
label_map <- data.frame(Member = specified_labels, Label = highlight_labels)

# Add label mapping to feat.samp.frame with a left join
feat.samp.frame <- left_join(feat.samp.frame, label_map, by = "Member")

# Use facet_plot to add both points and labels in the same panel
# Add points to the facet
p2 <- facet_plot(p, panel = "Sample Nearest Neighbor", data = feat.samp.frame, 
                 geom = geom_point, mapping = aes(x = Value, y = y, color = Direction), size=0.5)
p2
# Add labels within the same facet
## since I have many frankias and rhiziobiales labels that match the same OTU in different samples these will be duplicated so i want to keeep only one instance of it. 

# Create a new data frame specifically for plotting with single labels
feat.samp.frame_labeled <- feat.samp.frame %>%
  group_by(Member) %>%
  mutate(Label = if_else(row_number() == 1, Label, NA_character_)) %>%
  ungroup() # It's good practice to ungroup after the operation


library(ggrepel)
p3 <- facet_plot(p2, panel = "Sample Nearest Neighbor", data = feat.samp.frame_labeled, 
                 geom = geom_text_repel, mapping = aes(x = Value, y = y, label = Label),
                 hjust = -0.5, color = "blue", size=4, na.rm = TRUE)

p3

p4<- p3 + scale_color_manual(values = c("firebrick4", "firebrick1", "goldenrod1", "dodgerblue1", "dodgerblue4"), drop = F) + 
  theme_tree2()+theme(legend.position = "left")+
  theme(axis.text.x = element_text(size=16, face="bold"),
        legend.text = element_text(size = 16),  # Increase legend labels size
        legend.title = element_text(size = 16) ,
        strip.text = element_text(size = 16, face="bold"))


p4

ggsave("bnti_wle_transition_all_sites_16s_100perm_03.12.26_labels_ggrepel.TIFF", plot=p4, dpi=300, height=10, width=20, units=c("in"))


##merge the dataframe with bnti values with the sample metadata
metadata<- read.csv("sample-metadata-cb-wle.csv")
tax<- read.csv("taxonomy-dn-99.csv")
feat.samp.frame.metadata<- feat.samp.frame_labeled %>% left_join(metadata, by="Sample")
feat.samp.frame.metadata <- feat.samp.frame.metadata %>% mutate(OTUID=Member)
feat.samp.frame.metadata.tax<- feat.samp.frame.metadata %>% left_join(tax, by="OTUID")

feat.samp.frame.metadata.tax <- as.data.frame(feat.samp.frame.metadata.tax)
feat.samp.frame.metadata.tax <- feat.samp.frame.metadata.tax %>% mutate(Value_final=ifelse(is.na(Value), 0, Value))
class(feat.samp.frame.metadata.tax$Value_final)

##generate an average of all the bnti values by feature 
feat.metadata.tax.mean.bnti<- feat.samp.frame.metadata.tax %>% group_by(OTUID) %>%
  mutate(mean = mean(Value_final)) %>%
  ungroup() # It's good practice to ungroup after the operation

write.csv(feat.samp.frame.metadata, "feat.samp.frame.metadata.wle.transition.csv")
write.csv(feat.samp.frame.metadata.tax, "feat.samp.frame.metadata.tax.wle.transition.csv")
write.csv(feat.metadata.tax.mean.bnti, "feat.metadata.tax.mean.bnti.wle.transition.csv")



##keep only the orders in the network 
feat.metadata.tax.mean.bnti.subset <- feat.metadata.tax.mean.bnti %>% filter(Order==c("o__Rhizobiales", "o__Clostridiales"))
feat.metadata.tax.mean.bnti.subset

##i want to label the specific otus in the network
otus_to_label <- c("8f0fb68673a8b3c26f3fd210e29b7916", "d829bee4984f82ffc2453212157caf96")

# 2. Create a new data frame with a dedicated label column
data_for_plotting <- feat.metadata.tax.mean.bnti.subset %>%
  mutate(
    # Create the new column 'point_label'
    point_label = if_else(
      # Condition: is the OTUID in our target list?
      OTUID %in% otus_to_label, 
      # If TRUE, use the OTUID as the label
      as.character(OTUID),
      # If FALSE, use NA
      NA_character_
    )
  )

# Let's look at a few rows to see the new column
# Use glimpse() to see the structure and some data
glimpse(data_for_plotting)

data_for_plotting_edit <- data_for_plotting %>% filter(point_label!="NA")

### Create a new data frame with ONE row per 'Order'
plot_summary_data <- data_for_plotting_edit %>%
  # Group by the variable you want on your x-axis
  group_by(Order) %>% 
  
  # Summarize to get a single row for each Order
  summarize(
    # We use first() because the value is the same for all rows in the group
    mean_value = first(mean), 
    #site = first(site),
    transect = first(transect),
    point_label = first(point_label) # This also ensures you get only one label per group
  )

# View the new, clean data for plotting
print(plot_summary_data)

write.csv(plot_summary_data, "plot_summary_data_16s_wle_transition.csv")

bnti_plot<- ggplot(plot_summary_data, aes(Order, 
                                          mean_value, fill=transect
))+
  coord_cartesian(ylim = c(-20, 20)) +
  # Add the lines BEFORE your data points so they are in the background.
  # The color is SET, not mapped, so it will not affect the legend.
  
  # Red dashed lines at +2 and -2
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "red", size = 0.8) +
  
  # Blue dashed lines at +1 and -1
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue", size = 0.8) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "blue", size = 0.8) +
  
  
  # geom_boxplot(aes(color=transect),outlier.shape = NA,
  #   fill = "white", 
  #     outlier.colour = NA, 
  #  position = position_dodge(width = 0.9),
  #  show.legend = FALSE ) + 
  geom_point(aes(fill=transect),size=4, shape=21)+#position = position_jitterdodge(dodge.width = 0.9)) +
  #geom_text_repel(
  #aes(label = point_label),  # Use the new label column
  # position = position_jitterdodge(dodge.width = 0.9), # Dodge labels with points
  # size = 4,
  #color = "black",
  #box.padding = 0.5,         # Increase space between label and point
  # max.overlaps = Inf         # Ensure all labels are shown
  #) +
  #facet_grid(~site)+
  scale_fill_manual(values=c("#E69F00"))+
  #scale_color_manual(values=c("#CC79A7",  "#E69F00", "#0072B2"))+
  #scale_shape_manual(values=c(21,22,24))+
  labs(
    y = "βNTIfeat",
    x = "Order")+ # Correct x-axis label
  #  fill = "transect")+
  #ylab("βNTIfeat") + #change index name as needed
  #xlab("Site")+
  #guides(fill=guide_legend(override.aes = list(shape=21), order=1),
  #  shape = guide_legend(override.aes = list(fill = "black"), order=2)) +
  #labs(fill="transect", shape="site")+
  theme(axis.text.x = element_text(colour = "black", size = 12,angle = 90, vjust = 0.7, hjust=0.5),
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title.y= element_text(size=14, color="black"),
        axis.title.x= element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.title = element_text(size =16, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")) +
  theme(plot.title=element_text(size=16, face="bold", hjust=0.5))

bnti_plot


ggsave("bnti_wle_transition_all_sites_16s_100perm_03.12.26_points.TIFF", plot=bnti_plot, dpi=300, height=5, width=6, units=c("in"))









