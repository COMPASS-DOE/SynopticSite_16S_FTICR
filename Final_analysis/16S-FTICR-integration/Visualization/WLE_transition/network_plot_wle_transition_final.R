library(network)
library(tidyverse)
library(igraph)
library(GGally)

devtools::install_github("briatte/ggnet")
library(ggnet)
library(sna)
library(ggplot2)


########################################################################
#generating module figures selected from stats (r2 etc) and plotting only module relevant nodes and edges
#########################################################################

##first filter out grey color modules or ME0
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE transition 100 perm/TOM threshold zero")
mod.stats<- read.csv("network_stats_correlation_by_module_WLE_transition.csv") ##wle is not case sensitive
mod.stats.filter<- mod.stats %>% filter(Module!="ME0")


##EXAMINE THE MOD.FILTER FILE AND SELECT MODULES BASED ON R SQUARE VALUES 
##picking modules 33,8,26 because these modules have ASVs in them along with organic matter

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE transition 100 perm/TOM threshold zero/module specific nodes and edges")

mod33.nodes<- read.csv("Cytoscape_Nodes_ME33_0-thresh_signed.txt", sep="")

mod33.edges<- read.csv("Cytoscape_Edges_ME33_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod33.edges,
                         vertex.attr = mod33.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod33.edges, directed=FALSE)##changing directed to TRUE
x=network.vertex.names(net)
x=ifelse(nchar(x)>=25, "OTU", "molecule")
net %v% "type" = x
#net %v% "nodecolor" = c(OTU="yellow", molecule="blue")
#y2 = RColorBrewer::brewer.pal(9, "Set1")[ c(3, 1, 9,6, 7) ]

#names(y2)<- levels(x)
ggnet2(net, node.size =3, color = "type", palette = c("OTU" = "gold", "molecule" = "black"))
#ggnet2(net, node.size =1, color ="black")
y=network.vertex.names(net)
x.data<- as.data.frame(x)
y.data<- as.data.frame(y)
merge.x.y<- cbind(x.data, y.data)

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE transition 100 perm/TOM threshold zero")
tax<- read.csv("taxonomy-dn-99.csv")
icr_meta<- read.csv("icr_meta_final.csv")

merge.x.y<- merge.x.y %>% mutate(OTUID=y) 

merge.metadata<- left_join(merge.x.y, tax, by="OTUID") ##do this merge if OTU is present else skip
merge.metadata<- merge.metadata %>% mutate(formula=y) 
merge.metadata.icr<- left_join(merge.metadata, icr_meta, by="formula")


library(dplyr)
library(tidyverse)
library(readr)
library(tidyr)
merge.metadata.otu.icr<- merge.metadata.icr %>% mutate(value=ifelse(is.na(Order)==TRUE, Class_detailed, Order))
#merge.metadata.otu.icr<- filter(merge.metadata.otu.icr, !(Kingdom=="d__Bacteria" && Phylum==""))
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod33.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c("#a35151",'#636bb7',"#111b77")
# "#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
#"#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
#  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
#   "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
# "darkslategrey", "hotpink3", "indianred1")

#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_33_wle_transition<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_33_wle_transition



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod33.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod33.edges<- mod33.edges %>% mutate(y=toNode)
mod33.edges.edit<- mod33.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod33.edges.new<- arrange(mod33.edges.edit, toNode, desc(toNode))
unique(mod33.edges.new$toNode)

net<- network(mod33.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod33.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Micrococcales'='#a35151','aromatic'='#636bb7', 'unsaturated/lignin'='#111b77'))
mod33.final

ggsave("module_33_wle_transition.TIFF", plot=mod33.final, dpi=300, height=5, width=6, units=c("in"))

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE transition 100 perm/TOM threshold zero/module specific nodes and edges")

mod8.nodes<- read.csv("Cytoscape_Nodes_ME8_0-thresh_signed.txt", sep="")

mod8.edges<- read.csv("Cytoscape_Edges_ME8_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod8.edges,
                         vertex.attr = mod8.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod8.edges, directed=FALSE)##changing directed to TRUE
x=network.vertex.names(net)
x=ifelse(nchar(x)>=25, "OTU", "molecule")
net %v% "type" = x
#net %v% "nodecolor" = c(OTU="yellow", molecule="blue")
#y2 = RColorBrewer::brewer.pal(9, "Set1")[ c(3, 1, 9,6, 7) ]

#names(y2)<- levels(x)
ggnet2(net, node.size =3, color = "type", palette = c("OTU" = "gold", "molecule" = "black"))
#ggnet2(net, node.size =1, color ="black")
y=network.vertex.names(net)
x.data<- as.data.frame(x)
y.data<- as.data.frame(y)
merge.x.y<- cbind(x.data, y.data)

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE transition 100 perm/TOM threshold zero")
tax<- read.csv("taxonomy-dn-99.csv")
icr_meta<- read.csv("icr_meta_final.csv")

merge.x.y<- merge.x.y %>% mutate(OTUID=y) 

merge.metadata<- left_join(merge.x.y, tax, by="OTUID") ##do this merge if OTU is present else skip
merge.metadata<- merge.metadata %>% mutate(formula=y) 
merge.metadata.icr<- left_join(merge.metadata, icr_meta, by="formula")


library(dplyr)
library(tidyverse)
library(readr)
library(tidyr)
merge.metadata.otu.icr<- merge.metadata.icr %>% mutate(value=ifelse(is.na(Order)==TRUE, Class_detailed, Order))
#merge.metadata.otu.icr<- filter(merge.metadata.otu.icr, !(Kingdom=="d__Bacteria" && Phylum==""))
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod8.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c("darkslategrey","#111b77")
# "#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
#"#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
#  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
#   "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
# "darkslategrey", "hotpink3", "indianred1")

#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_8_wle_transition<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_8_wle_transition



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod8.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod8.edges<- mod8.edges %>% mutate(y=toNode)
mod8.edges.edit<- mod8.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod8.edges.new<- arrange(mod8.edges.edit, toNode, desc(toNode))
unique(mod8.edges.new$toNode)

net<- network(mod8.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod8.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Clostridiales'='darkslategrey','unsaturated/lignin'='#111b77'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod8.final

ggsave("module_8_wle_transition_weighted.TIFF", plot=mod8.final, dpi=300, height=5, width=6, units=c("in"))


##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE transition 100 perm/TOM threshold zero/module specific nodes and edges")

mod26.nodes<- read.csv("Cytoscape_Nodes_ME26_0-thresh_signed.txt", sep="")

mod26.edges<- read.csv("Cytoscape_Edges_ME26_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod26.edges,
                         vertex.attr = mod26.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod26.edges, directed=FALSE)##changing directed to TRUE
x=network.vertex.names(net)
x=ifelse(nchar(x)>=25, "OTU", "molecule")
net %v% "type" = x
#net %v% "nodecolor" = c(OTU="yellow", molecule="blue")
#y2 = RColorBrewer::brewer.pal(9, "Set1")[ c(3, 1, 9,6, 7) ]

#names(y2)<- levels(x)
ggnet2(net, node.size =3, color = "type", palette = c("OTU" = "gold", "molecule" = "black"))
#ggnet2(net, node.size =1, color ="black")
y=network.vertex.names(net)
x.data<- as.data.frame(x)
y.data<- as.data.frame(y)
merge.x.y<- cbind(x.data, y.data)

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE transition 100 perm/TOM threshold zero")
tax<- read.csv("taxonomy-dn-99.csv")
icr_meta<- read.csv("icr_meta_final.csv")

merge.x.y<- merge.x.y %>% mutate(OTUID=y) 

merge.metadata<- left_join(merge.x.y, tax, by="OTUID") ##do this merge if OTU is present else skip
merge.metadata<- merge.metadata %>% mutate(formula=y) 
merge.metadata.icr<- left_join(merge.metadata, icr_meta, by="formula")


library(dplyr)
library(tidyverse)
library(readr)
library(tidyr)
merge.metadata.otu.icr<- merge.metadata.icr %>% mutate(value=ifelse(is.na(Order)==TRUE, Class_detailed, Order))
#merge.metadata.otu.icr<- filter(merge.metadata.otu.icr, !(Kingdom=="d__Bacteria" && Phylum==""))
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod26.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c("#560d0d","#111b77","#636bb7")
# "#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
#"#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
#  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
#   "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
# "darkslategrey", "hotpink3", "indianred1")

#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_26_wle_transition<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_26_wle_transition



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod26.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod26.edges<- mod26.edges %>% mutate(y=toNode)
mod26.edges.edit<- mod26.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod26.edges.new<- arrange(mod26.edges.edit, toNode, desc(toNode))
unique(mod26.edges.new$toNode)

net<- network(mod26.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod26.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod26.final

ggsave("module_26_wle_transition_weighted.TIFF", plot=mod26.final, dpi=300, height=5, width=6, units=c("in"))




