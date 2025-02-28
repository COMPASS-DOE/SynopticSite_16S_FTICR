---
title: "Core analysis"
author: "Marta Nierychlo and Morten Simonsen Dueholm" 
date: "2021-08-25"
edits: "Sreejata Bandopadhyay"
edit_date: "2023-11-21"
---
#Load packages

install.packages("remotes")
remotes::install_github("kasperskytte/ampvis2")
library(ampvis2)
library(data.table)
library(tidyverse)
library(patchwork)
library(viridis)

########################
##Lake Erie#############
#######################

#Load data

otutab<- read.csv("OTU_clean_noneg_rarefied20k.csv", row.names=1)
taxonomy<- read.csv("taxonomy-dn-99.csv", row.names=1)
metadata<- read.csv("sample-metadata-cb-wle.csv")
v4 <- amp_load(otutab = otutab, 
                taxonomy = taxonomy,
                metadata = metadata)

v4n <- amp_subset_samples(v4, minreads = 20000, normalise = TRUE)

#Subset 
v4n_subset_wle <- amp_subset_samples(v4n, region == "Lake Erie ") #change location as needed


### Genus-level core and conditional abundant taxa for V1V3 data

data <- v4n_subset_wle
group_by <- "sample.id"

#Create taxonomy
tax <- v4n_subset_wle$tax[1:7]
tax$Kingdom <- gsub("d__","",tax$Kingdom)
tax$Phylum <- gsub("p__","",tax$Phylum)
tax$Class <- gsub("c__","",tax$Class)
tax$Order <- gsub("o__","",tax$Order)
tax$Family <- gsub("f__","",tax$Family)
tax$Genus <- gsub("g__","",tax$Genus)
tax$Species <- gsub("s__","",tax$Species)
tax <- distinct(tax)

d <- amp_export_long(
  data,
  #metadata_vars = group_by, keeping this as the default which includes all variables
  tax_levels = c("Genus"))

#group up and summarise for core taxa
gg <- d[, .(sum = sum(count)), by = c("Genus", group_by)]   
setorderv(gg, c(group_by, "sum"), order = -1)
#calculate proportion % abundance of each ASV
gg[, Genusprop := sum / sum(sum) * 100, by = group_by]
#calculate how many times given ASV is observed in plants (>0)
gg <- gg[Genus!=""]
gg[, nObs := sum(sum > 0), by = Genus]
#calculate how many times (in how many plants) ASV is >0.1% abundant
gg[, nCore := sum(Genusprop > 0.1), by = Genus]
#add column with sum of all the reads
gg[, totalReads := sum(sum)]
gg_summary <- gg[
  , 
  .(nGenus = uniqueN(Genus)),
  by = .(nObs, nCore)]

### Plot for core genera
#generate main plot
p1m <- ggplot(gg_summary,aes(
  x = as.numeric(nObs), #factors to align correctly with margin plots
  y = as.numeric(nCore), #factors to align correctly with margin plots
  weight = nGenus)) +
  geom_hex(bins = 100)+
  scale_x_continuous(breaks = seq(0,n_distinct(gg$sample.id),25), limits=c(-1,n_distinct(gg$sample.id)+1)) + #for LE changing the intervals to 25 instead of 50 for CB
  scale_y_continuous(breaks = seq(0,n_distinct(gg$sample.id),25), limits=c(-1,n_distinct(gg$sample.id)+1)) +#for LE changing the intervals to 25 instead of 50 for CB
  xlab(paste0("Observed in n samples")) +
  ylab(paste0("More than 0.1% abundant \nin n samples")) +
  scale_fill_viridis(option="turbo", trans = "log10", breaks = c(1, 10, 100, 1000), limits=c(1,1000)) + 
  theme_bw() +
  geom_hline(yintercept = n_distinct(gg$sample.id)*0.2, linetype="dashed", color = "darkgreen") +
  geom_hline(yintercept = n_distinct(gg$sample.id)*0.5, linetype="dashed", color = "darkgreen") +
  geom_hline(yintercept = n_distinct(gg$sample.id)*0.8, linetype="dashed", color = "darkgreen") +
  annotate("text", x=0, y=n_distinct(gg$sample.id)*0.25, label= "Abundant in greater than 20% of samples \n(loose core)",, fontface="bold", hjust = 0, color = "black") +
  annotate("text", x=0, y=n_distinct(gg$sample.id)*0.55, label= "Abundant in greater than 50% of samples \n(general core)",, fontface="bold", hjust = 0, color = "black") +
  annotate("text", x=0, y=n_distinct(gg$sample.id)*0.85, label= "Abundant in greater than 80% of samples \n(strict core)", fontface="bold", hjust = 0, color = "black") +
  theme(legend.position = "bottom")
p1m
#x margin plot  
p1x <- ggplot(gg[, .(nObsSum = sum(sum)/unique(totalReads)*100), by = .(nObs)][order(nObs)][,nObsSumCum:=cumsum(nObsSum)], aes(x = as.numeric(nObs), y= nObsSumCum)) +
  geom_line(width = 1) +
  ylab("Cumulative genus \nabundance (%)") +
  scale_x_continuous(breaks = seq(0,n_distinct(gg$sample.id),25), limits=c(-1,n_distinct(gg$sample.id)+1))+ #for LE changing the intervals to 25 instead of 50 for CB
  theme_bw() +
  theme(
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    panel.grid.minor.y = element_blank())
p1x

#y margin plot  
p1y <- ggplot(
  gg[, .(nCoreSum = sum(sum)/unique(totalReads)*100), by = .(nCore)][order(nCore)][,nCoreSumCum:=cumsum(nCoreSum)], aes(x = as.numeric(nCore), y = nCoreSumCum)) +
  geom_line(width=1) +
  ylab("Cumulative genus \nabundance (%)") +
  scale_x_continuous(breaks = seq(0,n_distinct(gg$sample.id),25), limits=c(-1,n_distinct(gg$sample.id)+1))+ #for LE changing the intervals to 25 instead of 50 for CB
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    #axis.ticks.y = element_blank(),
    #axis.text.y = element_blank(),
    panel.grid.minor.x = element_blank()) +
  coord_flip()+
  geom_vline(xintercept = n_distinct(gg$sample.id)*0.2, linetype="dashed", color = "darkgreen")+
  geom_vline(xintercept = n_distinct(gg$sample.id)*0.5, linetype="dashed", color = "darkgreen")+
  geom_vline(xintercept = n_distinct(gg$sample.id)*0.8, linetype="dashed", color = "darkgreen")
p1y
###Retrive core data
#subset ASVs present in the core
ggV4_1 <- gg[(gg$nCore >= n_distinct(gg$sample.id)*0.2),] %>%
  group_by(Genus) %>%
  summarise(mean_abu = mean(Genusprop)) %>%
  arrange(desc(mean_abu))

ggV4_1 <- ggV4_1[ ggV4_1$Genus != "", ]
ggV4_1$Genus <- gsub("g__","",ggV4_1$Genus)
ggV4_1[,"V4"] <- 2

ggV4_2 <- gg[(gg$nCore >= n_distinct(gg$sample.id)*0.5),] %>%
  group_by(Genus) %>%
  summarise(mean_abu = mean(Genusprop)) %>%
  arrange(desc(mean_abu))

ggV4_2 <- ggV4_2[ ggV4_2$Genus != "", ]
ggV4_2$Genus <- gsub("g__","",ggV4_2$Genus)
ggV4_2[,"V4"] <- 3

ggV4_3 <- gg[(gg$nCore >= n_distinct(gg$sample.id)*0.8),] %>%
  group_by(Genus) %>%
  summarise(mean_abu = mean(Genusprop)) %>%
  arrange(desc(mean_abu))

ggV4_3 <- ggV4_3[ ggV4_3$Genus != "", ]
ggV4_3$Genus <- gsub("g__","",ggV4_3$Genus)
ggV4_3[,"V4"] <- 4

### Plot for conditional abundant genera

#group up and summarise for conditional abundant taxa (CAT)
gg2 <- d[, .(sum = sum(count)), by = c("Genus", group_by)]   
setorderv(gg2, c(group_by, "sum"), order = -1)
#calculate proportion % abundance of each ASV
gg2[, Genusprop := sum / sum(sum) * 100, by = group_by]
#calculate how many times given ASV is observed in plants (>0)
gg2 <- gg2[Genus!=""]
gg2 <- gg2[!(substr(Genus,4,100) %in% ggV4_1$Genus)] # Remove core genera
gg2[, nObs := sum(sum > 0), by = Genus]
#calculate how many times (in how many plants) ASV is >1% abundant
gg2[, nCA := sum(Genusprop > 1), by = Genus]
#add column with sum of all the reads
gg2[, totalReads := sum(sum)]
gg2_summary <- gg2[,.(nGenus = uniqueN(Genus)), by = .(nObs, nCA)]

#generate main plot
p2m <- ggplot(gg2_summary,aes(
  x = as.numeric(nObs), #factors to align correctly with margin plots
  y = as.numeric(nCA), #factors to align correctly with margin plots
  weight = nGenus)) +
  geom_hex(bins = 100)+
  scale_x_continuous(breaks = seq(0,n_distinct(gg$sample.id),25), limits=c(-1,n_distinct(gg$sample.id)+1)) +#for LE changing the intervals to 25 instead of 50 for CB
  #scale_y_continuous(breaks = seq(0,50,5), limits=c(-1,50+1)) +
  xlab(paste0("Observed in n samples")) +
  ylab(paste0("More than 1% abundant in n samples")) +
  scale_fill_viridis(option="turbo", trans = "log10", breaks = c(1, 10, 100, 1000), limits=c(1,1000)) + 
  theme_bw() +
  geom_hline(yintercept = 1, linetype="dashed", color = "darkgreen") +
 annotate("text", x=0, y=3, label= "More than 1% abundant in at least one sample \n(conditionally rare or abundant genera)", fontface="bold", hjust = 0, color = "black") + ##changing y from 7.5 to 3 for LE
  theme(legend.position = "bottom")
p2m
#x margin plot  
p2x <- ggplot(gg2[, .(nObsSum = sum(sum)/unique(totalReads)*100), by = .(nObs)][order(nObs)][,nObsSumCum:=cumsum(nObsSum)], aes(x = as.numeric(nObs), y= nObsSumCum)) +
  geom_line() +
  ylab("Cumulative genus \nabundance (%)") +
 scale_x_continuous(breaks = seq(0,n_distinct(gg$sample.id),25), limits=c(-1,n_distinct(gg$sample.id)+1)) +
 # scale_y_continuous(limits=c(0,100)) +
  theme_bw() +
  theme(
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    panel.grid.minor.y = element_blank())
p2x
#y margin plot  
p2y <- ggplot(gg2[, .(nCASum = sum(sum)/unique(totalReads)*100), by = .(nCA)][order(nCA)][,nCASumCum:=cumsum(nCASum)], aes(x = as.numeric(nCA), y = nCASumCum)) +
  geom_line(width=1) +
  ylab("Cumulative genus \nabundance (%)") +
  #scale_x_continuous(breaks = seq(0,n_distinct(gg$sample.id),50), limits=c(-1,n_distinct(gg$sample.id)+1)) +
  #scale_x_continuous(breaks = seq(0,50,5), limits=c(-1,50+1)) +
 scale_y_continuous(limits=c(0,101)) +    ##the WLE plot was removing one data point if the limiot was not set to 101 even though max was 100 (I checked)
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(),
    #axis.ticks.y = element_blank(),
   # axis.text.y = element_blank(),
    panel.grid.minor.x = element_blank()) +
  coord_flip()+
  geom_vline(xintercept = 1, linetype="dashed", color = "darkgreen")
p2y
# Create combined plot
p1 <- p1x + plot_spacer() + p2x + plot_spacer() +
  p1m + p1y + p2m + p2y +
  plot_layout(ncol=4, widths = c(4,2,4,2), heights = c(1.5,4)) & theme(text = element_text(size=14))
p1
ggsave(filename="Core_genera_V4_WLE_new_turbo.pdf", plot=p1, width=17, height=8, useDingbats=FALSE, limitsize=FALSE)

###Retrieve core data
#subset ASVs present in the core
ggV4_4 <- gg2[(gg2$nCA >= 1)] %>%
  group_by(Genus) %>%
  summarise(mean_abu = mean(Genusprop)) %>%
  arrange(desc(mean_abu))

ggV4_4 <- ggV4_4[ ggV4_4$Genus != "", ]
ggV4_4$Genus <- gsub("g__","",ggV4_4$Genus)
ggV4_4[,"V4"] <- 1

#Merge core data for v4
ggV4_4 <- filter(ggV4_4, !(Genus %in% ggV4_3$Genus))
ggV4_1 <- filter(ggV4_1, !(Genus %in% ggV4_2$Genus))
ggV4_2 <- filter(ggV4_2, !(Genus %in% ggV4_3$Genus))


V4_core <- rbind(ggV4_3[,c(1,3)], ggV4_2[,c(1,3)]) %>%
  rbind(., ggV4_1[,c(1,3)]) %>%
  rbind(., ggV4_4[,c(1,3)])

V4_core<- V4_core%>% mutate(type=ifelse(V4==1, "Cond.Rare.Abundant", "Core"))
write.csv(V4_core, "Core_CondAbun_genera_WLE_checkingcorrectness.csv")

################
###Chesapeake bay
#################


#Subset 
v4n_subset_cb <- amp_subset_samples(v4n, region == "Chesapeake Bay") #change location as needed  


### Genus-level core and conditional abundant taxa for V1V3 data

data <- v4n_subset_cb
group_by <- "sample.id"

#Create taxonomy
tax <- v4n_subset_cb$tax[1:7]
tax$Kingdom <- gsub("d__","",tax$Kingdom)
tax$Phylum <- gsub("p__","",tax$Phylum)
tax$Class <- gsub("c__","",tax$Class)
tax$Order <- gsub("o__","",tax$Order)
tax$Family <- gsub("f__","",tax$Family)
tax$Genus <- gsub("g__","",tax$Genus)
tax$Species <- gsub("s__","",tax$Species)
tax <- distinct(tax)

d <- amp_export_long(
  data,
  #metadata_vars = group_by, keeping this as the default which includes all variables
  tax_levels = c("Genus"))

#group up and summarise for core taxa
gg <- d[, .(sum = sum(count)), by = c("Genus", group_by)]   
setorderv(gg, c(group_by, "sum"), order = -1)
#calculate proportion % abundance of each ASV
gg[, Genusprop := sum / sum(sum) * 100, by = group_by]
#calculate how many times given ASV is observed in plants (>0)
gg <- gg[Genus!=""]
gg[, nObs := sum(sum > 0), by = Genus]
#calculate how many times (in how many plants) ASV is >0.1% abundant
gg[, nCore := sum(Genusprop > 0.1), by = Genus]
#add column with sum of all the reads
gg[, totalReads := sum(sum)]
gg_summary <- gg[
  , 
  .(nGenus = uniqueN(Genus)),
  by = .(nObs, nCore)]

### Plot for core genera
#generate main plot
q1m <- ggplot(gg_summary,aes(
  x = as.numeric(nObs), #factors to align correctly with margin plots
  y = as.numeric(nCore), #factors to align correctly with margin plots
  weight = nGenus)) +
  geom_hex(bins = 100)+
  scale_x_continuous(breaks = seq(0,n_distinct(gg$sample.id),50), limits=c(-1,n_distinct(gg$sample.id)+1)) + #for LE changing the intervals to 25 instead of 50 for CB
  scale_y_continuous(breaks = seq(0,n_distinct(gg$sample.id),50), limits=c(-1,n_distinct(gg$sample.id)+1)) +#for LE changing the intervals to 25 instead of 50 for CB
  xlab(paste0("Observed in n samples")) +
  ylab(paste0("More than 0.1% abundant \nin n samples")) +
  scale_fill_viridis(option="turbo", trans = "log10", breaks = c(1, 10, 100, 1000), limits=c(1,1000)) + 
  theme_bw() +
  geom_hline(yintercept = n_distinct(gg$sample.id)*0.2, linetype="dashed", color = "darkgreen") +
  geom_hline(yintercept = n_distinct(gg$sample.id)*0.5, linetype="dashed", color = "darkgreen") +
  geom_hline(yintercept = n_distinct(gg$sample.id)*0.8, linetype="dashed", color = "darkgreen") +
  annotate("text", x=0, y=n_distinct(gg$sample.id)*0.25, label= "Abundant in greater than 20% of samples \n(loose core)",, fontface="bold", hjust = 0, color = "black") +
  annotate("text", x=0, y=n_distinct(gg$sample.id)*0.55, label= "Abundant in greater than 50% of samples \n(general core)",, fontface="bold", hjust = 0, color = "black") +
  annotate("text", x=0, y=n_distinct(gg$sample.id)*0.85, label= "Abundant in greater than 80% of samples \n(strict core)", fontface="bold", hjust = 0, color = "black") +
  theme(legend.position = "bottom")
q1m
#x margin plot  
q1x <- ggplot(gg[, .(nObsSum = sum(sum)/unique(totalReads)*100), by = .(nObs)][order(nObs)][,nObsSumCum:=cumsum(nObsSum)], aes(x = as.numeric(nObs), y= nObsSumCum)) +
  geom_line(width = 1) +
  ylab("Cumulative genus \nabundance (%)") +
  scale_x_continuous(breaks = seq(0,n_distinct(gg$sample.id),50), limits=c(-1,n_distinct(gg$sample.id)+1))+ #for LE changing the intervals to 25 instead of 50 for CB
  theme_bw() +
  theme(
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    panel.grid.minor.y = element_blank())
q1x

#y margin plot  
q1y <- ggplot(
  gg[, .(nCoreSum = sum(sum)/unique(totalReads)*100), by = .(nCore)][order(nCore)][,nCoreSumCum:=cumsum(nCoreSum)], aes(x = as.numeric(nCore), y = nCoreSumCum)) +
  geom_line(width=1) +
  ylab("Cumulative genus \nabundance (%)") +
  scale_x_continuous(breaks = seq(0,n_distinct(gg$sample.id),50), limits=c(-1,n_distinct(gg$sample.id)+1))+ #for LE changing the intervals to 25 instead of 50 for CB
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    #axis.ticks.y = element_blank(),
    #axis.text.y = element_blank(),
    panel.grid.minor.x = element_blank()) +
  coord_flip()+
  geom_vline(xintercept = n_distinct(gg$sample.id)*0.2, linetype="dashed", color = "darkgreen")+
  geom_vline(xintercept = n_distinct(gg$sample.id)*0.5, linetype="dashed", color = "darkgreen")+
  geom_vline(xintercept = n_distinct(gg$sample.id)*0.8, linetype="dashed", color = "darkgreen")
q1y
###Retrive core data
#subset ASVs present in the core
ggV4_1 <- gg[(gg$nCore >= n_distinct(gg$sample.id)*0.2),] %>%
  group_by(Genus) %>%
  summarise(mean_abu = mean(Genusprop)) %>%
  arrange(desc(mean_abu))

ggV4_1 <- ggV4_1[ ggV4_1$Genus != "", ]
ggV4_1$Genus <- gsub("g__","",ggV4_1$Genus)
ggV4_1[,"V4"] <- 2

ggV4_2 <- gg[(gg$nCore >= n_distinct(gg$sample.id)*0.5),] %>%
  group_by(Genus) %>%
  summarise(mean_abu = mean(Genusprop)) %>%
  arrange(desc(mean_abu))

ggV4_2 <- ggV4_2[ ggV4_2$Genus != "", ]
ggV4_2$Genus <- gsub("g__","",ggV4_2$Genus)
ggV4_2[,"V4"] <- 3

ggV4_3 <- gg[(gg$nCore >= n_distinct(gg$sample.id)*0.8),] %>%
  group_by(Genus) %>%
  summarise(mean_abu = mean(Genusprop)) %>%
  arrange(desc(mean_abu))

ggV4_3 <- ggV4_3[ ggV4_3$Genus != "", ]
ggV4_3$Genus <- gsub("g__","",ggV4_3$Genus)
ggV4_3[,"V4"] <- 4

### Plot for conditional abundant genera

#group up and summarise for conditional abundant taxa (CAT)
gg2 <- d[, .(sum = sum(count)), by = c("Genus", group_by)]   
setorderv(gg2, c(group_by, "sum"), order = -1)
#calculate proportion % abundance of each ASV
gg2[, Genusprop := sum / sum(sum) * 100, by = group_by]
#calculate how many times given ASV is observed in plants (>0)
gg2 <- gg2[Genus!=""]
gg2 <- gg2[!(substr(Genus,4,100) %in% ggV4_1$Genus)] # Remove core genera
gg2[, nObs := sum(sum > 0), by = Genus]
#calculate how many times (in how many plants) ASV is >1% abundant
gg2[, nCA := sum(Genusprop > 1), by = Genus]
#add column with sum of all the reads
gg2[, totalReads := sum(sum)]
gg2_summary <- gg2[,.(nGenus = uniqueN(Genus)), by = .(nObs, nCA)]

#generate main plot
q2m <- ggplot(gg2_summary,aes(
  x = as.numeric(nObs), #factors to align correctly with margin plots
  y = as.numeric(nCA), #factors to align correctly with margin plots
  weight = nGenus)) +
  geom_hex(bins = 100)+
  scale_x_continuous(breaks = seq(0,n_distinct(gg$sample.id),50), limits=c(-1,n_distinct(gg$sample.id)+1)) +#for LE changing the intervals to 25 instead of 50 for CB
  #scale_y_continuous(breaks = seq(0,50,5), limits=c(-1,50+1)) +
  xlab(paste0("Observed in n samples")) +
  ylab(paste0("More than 1% abundant in n samples")) +
  scale_fill_viridis(option="turbo", trans = "log10", breaks = c(1, 10, 100, 1000), limits=c(1,1000)) + 
  theme_bw() +
  geom_hline(yintercept = 1, linetype="dashed", color = "darkgreen") +
  annotate("text", x=0, y=7.5, label= "More than 1% abundant in at least one sample \n(conditionally rare or abundant genera)", fontface="bold", hjust = 0, color = "black") + ##changing y from 7.5 to 3 for LE
  theme(legend.position = "bottom")
q2m
#x margin plot  
q2x <- ggplot(gg2[, .(nObsSum = sum(sum)/unique(totalReads)*100), by = .(nObs)][order(nObs)][,nObsSumCum:=cumsum(nObsSum)], aes(x = as.numeric(nObs), y= nObsSumCum)) +
  geom_line() +
  ylab("Cumulative genus \nabundance (%)") +
  scale_x_continuous(breaks = seq(0,n_distinct(gg$sample.id),50), limits=c(-1,n_distinct(gg$sample.id)+1)) +
  # scale_y_continuous(limits=c(0,100)) +
  theme_bw() +
  theme(
    axis.title.y = element_text(),
    axis.title.x = element_blank(),
    #axis.ticks.x = element_blank(),
    #axis.text.x = element_blank(),
    panel.grid.minor.y = element_blank())
q2x
#y margin plot  
q2y <- ggplot(gg2[, .(nCASum = sum(sum)/unique(totalReads)*100), by = .(nCA)][order(nCA)][,nCASumCum:=cumsum(nCASum)], aes(x = as.numeric(nCA), y = nCASumCum)) +
  geom_line(width=1) +
  ylab("Cumulative genus \nabundance (%)") +
  #scale_x_continuous(breaks = seq(0,n_distinct(gg$sample.id),50), limits=c(-1,n_distinct(gg$sample.id)+1)) +
  #scale_x_continuous(breaks = seq(0,50,5), limits=c(-1,50+1)) +
  scale_y_continuous(limits=c(0,100)) +    ##the WLE plot was removing one data point if the limiot was not set to 101 even though max was 100 (I checked)
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_text(),
    #axis.ticks.y = element_blank(),
    # axis.text.y = element_blank(),
    panel.grid.minor.x = element_blank()) +
  coord_flip()+
  geom_vline(xintercept = 1, linetype="dashed", color = "darkgreen")
q2y


q1 <- q1x + plot_spacer() + q2x + plot_spacer() +
  q1m + q1y + q2m + q2y +
  plot_layout(ncol=4, widths = c(4,2,4,2), heights = c(1.5,4)) & theme(text = element_text(size=14))
q1
ggsave(filename="Core_genera_V4_CB_new_turbo.pdf", plot=q1, width=17, height=8, useDingbats=FALSE, limitsize=FALSE)


###Retrieve core data
#subset ASVs present in the core
ggV4_4 <- gg2[(gg2$nCA >= 1)] %>%
  group_by(Genus) %>%
  summarise(mean_abu = mean(Genusprop)) %>%
  arrange(desc(mean_abu))

ggV4_4 <- ggV4_4[ ggV4_4$Genus != "", ]
ggV4_4$Genus <- gsub("g__","",ggV4_4$Genus)
ggV4_4[,"V4"] <- 1

#Merge core data for v4
ggV4_4 <- filter(ggV4_4, !(Genus %in% ggV4_3$Genus))
ggV4_1 <- filter(ggV4_1, !(Genus %in% ggV4_2$Genus))
ggV4_2 <- filter(ggV4_2, !(Genus %in% ggV4_3$Genus))


V4_core <- rbind(ggV4_3[,c(1,3)], ggV4_2[,c(1,3)]) %>%
  rbind(., ggV4_1[,c(1,3)]) %>%
  rbind(., ggV4_4[,c(1,3)])

V4_core<- V4_core%>% mutate(type=ifelse(V4==1, "Cond.Rare.Abundant", "Core"))
write.csv(V4_core, "Core_CondAbun_genera_CB_checkingcorrectness.csv")


##combing the LE and CB graphs

core <- q1x + plot_spacer() + p1x + plot_spacer() +
  q1m + q1y + p1m + p1y +
  plot_layout(ncol=4, widths = c(4,2,4,2), heights = c(1.5,4)) & theme(text = element_text(size=14))+
  plot_annotation(title="Western Lake Erie", theme=)
core

ggsave(filename="Core_genera_V4_CB_WLE_new_turbo.pdf", plot=core, width=17, height=8, useDingbats=FALSE, limitsize=FALSE)


crat<- q2x +plot_spacer() + p2x +plot_spacer()+
  q2m + q2y + p2m + p2y+
  plot_layout(ncol=4, widths = c(4,2,4,2), heights = c(1.5,4)) & theme(text = element_text(size=14))
crat
ggsave(filename="CRAT_genera_V4_CB_WLE_new_turbo.pdf", plot=crat, width=17, height=8, useDingbats=FALSE, limitsize=FALSE)


####################################################################
###plotting relative abundance bubble plot using the core genera only 
###################################################################

##Chesapeake Bay

otu_final<- read.csv("OTU_CB_rarefied20k.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #5976 otus with no counts
View(otu1)
max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
metadata <- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)
sample_data_new<- metadata %>% mutate(transect1= ifelse(transect=="wetland-center", "wetland", transect))

sample_data_new$transect1 <- factor(sample_data_new$transect1,levels = c("upland", "transition", "wetland", "wetland-transition-edge"), ordered=TRUE)


sample_data_new_cb<- sample_data_new %>% filter(region=="Chesapeake Bay")

##changing name of sites within CB

sample_data_new_cb<- sample_data_new_cb %>% mutate(site_new=recode(site, "GCREW"="GCW"))

tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(sample_data_new_cb)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 24964 taxa and 117 samples ]
#sample_data() Sample Data:       [ 117 samples by 17 sample variables ]
#tax_table()   Taxonomy Table:    [ 24964 taxa by 8 taxonomic ranks ]



##barplot class level
phyloseq_class <- phyloseq_merged %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Genus) # Sort data frame alphabetically by phylum/family etc

#phyloseq_class$Genus[phyloseq_class$Abundance < 0.05] <- "< 5% abund."

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




######################
##dot plots, using the code below
###################

##for dotplots 
phyloseq_class1<- phyloseq_class %>% group_by(site, transect1, Genus) %>% mutate(mean=mean(Abundance))


#Filter out less than 1% abundance 

#phyloseq_class2<- phyloseq_class1 %>% filter(Class!="< 5% abund.")

#View(phyloseq_class2) 

#write.csv(phyloseq_class2, "Rel_abundance_CB_agglom_genus_classgreaterthan5%.csv") ##edited columns (drought and summarise new) so that pre-drought component exists for both drought and watered
#phyloseq_class2<- read.csv("Rel_abundance_CB_agglom_genus_classgreaterthan5%.csv")


##left join core taxa csv file with the above phyloseq_class2 object to plot only the strict and general core taxa. 

core_cb<- read.csv("core_crat_cb.csv", header=TRUE)
core_cb_edit<- core_cb %>% mutate(Genus_new=paste0("g__", Genus))
core_cb_edit<- core_cb_edit %>% mutate(Genus_old=Genus) %>% dplyr::select(!Genus)
core_cb_edit<- core_cb_edit %>% mutate(Genus= Genus_new) %>% dplyr::select (!Genus_new)
join_core<- phyloseq_class1 %>% right_join(core_cb_edit, by="Genus")

unique(join_core$Genus)

#Filtering out low abundance taxa
join_core_filter<- join_core%>% filter(mean>=0.01)
join_core_filter_new<- join_core_filter %>% filter(type!="Cond.Rare.Abundant")


library(ggh4x)

library(viridis)
bb<-c(0.160, 0.100,0.050,0.010) # define breaks.
ll<-c("0.160", "0.100","0.050","0.010") # labels. ##all samples representation checked for bean planted drought and watered such that the range does not miss any samples
plot_dotplot_cb_coretaxa<-
  ggplot(join_core_filter_new,aes(x=transect1,y=Genus,size=mean,color=mean))+
  theme_bw()+
  facet_nested(~type+site_new)+
  geom_point()+
  scale_x_discrete(limits = c("upland", "transition", "wetland"))+
  #facet_grid(~type+site,scales="free")+ 
  theme(strip.text.x=element_text(size=12))+
  #scale_colour_gradient(low="#132B43",
  #high="#56B1F7",limits=c(0.001,0.3),breaks=bb)+
  scale_colour_viridis(limits=c(0.010,0.160),breaks=bb)+
  scale_size_continuous(limits=c(0.010,0.160),labels=ll,breaks=bb, name="Mean Relative\nAbundance")+
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1))+
  #ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(color="Mean Relative\nAbundance")+
  ylab("Genus") + 
  theme(axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=14)) + 
  theme(legend.title = element_text(size=16))+ theme(legend.text = element_text(size=12))

plot_dotplot_cb_coretaxa

##Figure 
ggsave(filename = "dotplot_cb_coretaxa_redofinal.tiff", plot = plot_dotplot_cb_coretaxa,
       width = 30 ,
       height = 20, units = c("cm"),
       dpi = 200)

ggsave(filename = "dotplot_cb_coretaxa_redofinal_02.14.25.tiff", plot = plot_dotplot_cb_coretaxa,
       width = 30 ,
       height = 20, units = c("cm"),
       dpi = 200)


###Western Lake Erie 

otu_final<- read.csv("OTU_WLE_rarefied20k.csv", sep=",",row.names=1)

any(rowSums(otu_final[])<1)

otu1 = otu_final[rowSums(otu_final[])<1, ,drop=FALSE] #5976 otus with no counts
View(otu1)
max(rowSums(otu1))
otu_final = otu_final[rowSums(otu_final[])>0, ,drop=FALSE] ##drop taxa that are not present in any samples, 40216 taxa remains
View(otu_final)

library(phyloseq)


tax<- read.csv("taxonomy-dn-99.csv", sep=",", row.names=1)
metadata <- read.csv("sample-metadata-cb-wle.csv", sep=",", row.names=1)
sample_data_new<- metadata %>% mutate(transect1= ifelse(transect=="wetland-center"|transect=="wetland-transition-edge", "wetland", transect))

sample_data_new$transect1 <- factor(sample_data_new$transect1,levels = c("upland", "transition", "wetland"), ordered=TRUE)


sample_data_new_wle<- sample_data_new %>% filter(region=="Lake Erie ")

##changing name of sites within WLE

sample_data_new_wle<- sample_data_new_wle %>% mutate(site_new=recode(site, "PR"="PTR", "CC"="CRC"))

tax = as.matrix(tax) 

otu<- otu_table(otu_final, taxa_are_rows = TRUE)
tax<- tax_table(tax)
metadata<- sample_data(sample_data_new_wle)

colnames(tax)
tax
phyloseq_merged<- merge_phyloseq(otu, tax, metadata)
phyloseq_merged

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 19262 taxa and 82 samples ]
#sample_data() Sample Data:       [ 82 samples by 17 sample variables ]
#tax_table()   Taxonomy Table:    [ 19262 taxa by 8 taxonomic ranks ]


##barplot class level
phyloseq_class <- phyloseq_merged %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at family/phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa
  arrange(Genus) # Sort data frame alphabetically by phylum/family etc

#phyloseq_class$Genus[phyloseq_class$Abundance < 0.05] <- "< 5% abund."

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


##for dotplots 
phyloseq_class1<- phyloseq_class %>% group_by(site, transect1, Genus) %>% mutate(mean=mean(Abundance))


##left join core taxa csv file with the above phyloseq_class2 object to plot only the strict and general core taxa. 

core_wle<- read.csv("core_crat_wle.csv", header=TRUE)
core_wle_edit<- core_wle %>% mutate(Genus_new=paste0("g__", Genus))
core_wle_edit<- core_wle_edit %>% mutate(Genus_old=Genus) %>% dplyr::select(!Genus)
core_wle_edit<- core_wle_edit %>% mutate(Genus= Genus_new) %>% dplyr::select (!Genus_new)
join_core<- phyloseq_class1 %>% right_join(core_wle_edit, by="Genus")

unique(join_core$Genus)
#Filtering out low abundance taxa
join_core_filter<- join_core%>% filter(mean>=0.01) ##changed from 0.01 to match the 0.1% abundance on the other core taxa figure but its too many to plot
join_core_filter_new<- join_core_filter %>% filter(type!="Cond.Rare.Abundant")

unique(join_core_filter$Genus)
unique(join_core_filter_new$Genus)

library(viridis)
bb<-c(0.120, 0.100,0.050,0.010) # define breaks.
ll<-c("0.120", "0.100","0.050","0.010") # labels. 

##not using these breaks below at 0.001 because its too many to plot
#bb<-c(0.120, 0.100,0.075, 0.050,0.025, 0.001) # define breaks.
#ll<-c("0.120","0.100","0.075", "0.050","0.025", "0.001") # labels. ##all samples representation checked for bean planted drought and watered such that the range does not miss any samples
plot_dotplot_wle_coretaxa<-
  ggplot(join_core_filter_new,aes(x=transect1,y=Genus,size=mean,color=mean))+
  theme_bw()+
  facet_nested(~type+site_new)+
  geom_point()+
  scale_x_discrete(limits = c("upland", "transition", "wetland"))+
  #facet_grid(~type+site,scales="free")+ 
  theme(strip.text.x=element_text(size=12))+
  #scale_colour_gradient(low="#132B43",
  #high="#56B1F7",limits=c(0.001,0.3),breaks=bb)+
  scale_colour_viridis(limits=c(0.010,0.120),breaks=bb)+
  scale_size_continuous(limits=c(0.010,0.120),labels=ll,breaks=bb, name="Mean Relative\nAbundance")+
  theme(axis.title.x = element_blank()) + theme(axis.text.x = element_text(size=12, angle = 90, vjust = 0.5, hjust=1))+
  #ggtitle("")+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(color="Mean Relative\nAbundance")+
  ylab("Genus") + 
  theme(axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=14)) + 
  theme(legend.title = element_text(size=16))+ theme(legend.text = element_text(size=12))

plot_dotplot_wle_coretaxa

##Figure 
ggsave(filename = "dotplot_wle_coretaxa_final.tiff", plot = plot_dotplot_wle_coretaxa,
       width = 30 ,
       height = 20, units = c("cm"),
       dpi = 200)

ggsave(filename = "dotplot_wle_coretaxa_final_02.14.25.tiff", plot = plot_dotplot_wle_coretaxa,
       width = 30 ,
       height = 20, units = c("cm"),
       dpi = 200)
