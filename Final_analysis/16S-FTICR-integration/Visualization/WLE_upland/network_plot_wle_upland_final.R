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
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero")
mod.stats<- read.csv("network_stats_correlation_by_module_WLE_upland.csv")
mod.stats.filter<- mod.stats %>% filter(Module!="ME0")


##EXAMINE THE MOD.FILTER FILE AND SELECT MODULES BASED ON R SQUARE VALUES 
##picking modules, 32,44,63,33,15,56,67 because these modules have ASVs in them along with organic matter

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero/module specific nodes and edges")

mod32.nodes<- read.csv("Cytoscape_Nodes_ME32_0-thresh_signed.txt", sep="")

mod32.edges<- read.csv("Cytoscape_Edges_ME32_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod32.edges,
                         vertex.attr = mod32.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod32.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod32.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c("#1e0047","#111b77","#636bb7")
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

module_32_wle_upland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_32_wle_upland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod32.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod32.edges<- mod32.edges %>% mutate(y=toNode)
mod32.edges.edit<- mod32.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod32.edges.new<- arrange(mod32.edges.edit, toNode, desc(toNode))
unique(mod32.edges.new$toNode)

net<- network(mod32.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod32.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Chthoniobacterales'='#1e0047','unsaturated/lignin'='#111b77','aromatic'='#636bb7'))
mod32.final

ggsave("module_32_wle_upland.TIFF", plot=mod32.final, dpi=300, height=5, width=6, units=c("in"))

##module 44

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero/module specific nodes and edges")

mod44.nodes<- read.csv("Cytoscape_Nodes_ME44_0-thresh_signed.txt", sep="")

mod44.edges<- read.csv("Cytoscape_Edges_ME44_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod44.edges,
                         vertex.attr = mod44.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod44.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod44.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c("#a0fffc","#111b77")
# "#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
# "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
#"#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
#   "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
# "darkslategrey", "hotpink3", "indianred1")

#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_44_wle_upland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_44_wle_upland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod44.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod44.edges<- mod44.edges %>% mutate(y=toNode)
mod44.edges.edit<- mod44.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod44.edges.new<- arrange(mod44.edges.edit, toNode, desc(toNode))
unique(mod44.edges.new$toNode)

net<- network(mod44.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod44.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Corynebacteriales'='#a0fffc','unsaturated/lignin'='#111b77'))
mod44.final

ggsave("module_44_wle_upland.TIFF", plot=mod44.final, dpi=300, height=5, width=6, units=c("in"))


##module 63


##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero/module specific nodes and edges")

mod63.nodes<- read.csv("Cytoscape_Nodes_ME63_0-thresh_signed.txt", sep="")

mod63.edges<- read.csv("Cytoscape_Edges_ME63_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod63.edges,
                         vertex.attr = mod63.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod63.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod63.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c("#D14285","#111b77","#636bb7")
# "#560d0d","#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
# "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
#  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
#   "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
# "darkslategrey", "hotpink3", "indianred1")

#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_63_wle_upland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_63_wle_upland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod63.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod63.edges<- mod63.edges %>% mutate(y=toNode)
mod63.edges.edit<- mod63.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod63.edges.new<- arrange(mod63.edges.edit, toNode, desc(toNode))
unique(mod63.edges.new$toNode)

net<- network(mod63.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod63.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Streptomycetales'='#D14285', 'unsaturated/lignin'='#111b77',  'aromatic'='#636bb7'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod63.final

ggsave("module_63_wle_upland_weighted.TIFF", plot=mod63.final, dpi=300, height=5, width=6, units=c("in"))



##module 33


##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero/module specific nodes and edges")

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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero")
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

palette14= c("#560d0d","#111b77","#636bb7")
# "#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
# "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
#  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
#   "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
# "darkslategrey", "hotpink3", "indianred1")

#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_33_wle_upland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_33_wle_upland



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

mod33.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod33.final

ggsave("module_33_wle_upland_weighted.TIFF", plot=mod33.final, dpi=300, height=5, width=6, units=c("in"))



##module 15



##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero/module specific nodes and edges")

mod15.nodes<- read.csv("Cytoscape_Nodes_ME15_0-thresh_signed.txt", sep="")

mod15.edges<- read.csv("Cytoscape_Edges_ME15_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod15.edges,
                         vertex.attr = mod15.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod15.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod15.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c("forestgreen","#111b77","#636bb7")
#"#560d0d", "#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
# "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
# "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
  #"#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
#"darkslategrey", "hotpink3", "indianred1")

#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_15_wle_upland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_15_wle_upland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod15.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod15.edges<- mod15.edges %>% mutate(y=toNode)
mod15.edges.edit<- mod15.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod15.edges.new<- arrange(mod15.edges.edit, toNode, desc(toNode))
unique(mod15.edges.new$toNode)

net<- network(mod15.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod15.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__KD4-96'='forestgreen','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7'))
mod15.final

ggsave("module_15_wle_upland.TIFF", plot=mod15.final, dpi=300, height=5, width=6, units=c("in"))



##module 56


##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero/module specific nodes and edges")

mod56.nodes<- read.csv("Cytoscape_Nodes_ME56_0-thresh_signed.txt", sep="")

mod56.edges<- read.csv("Cytoscape_Edges_ME56_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod56.edges,
                         vertex.attr = mod56.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod56.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod56.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c("#560d0d","#111b77")
# "#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
# "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
#  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
#   "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
# "darkslategrey", "hotpink3", "indianred1")

#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_56_wle_upland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_56_wle_upland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod56.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod56.edges<- mod56.edges %>% mutate(y=toNode)
mod56.edges.edit<- mod56.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod56.edges.new<- arrange(mod56.edges.edit, toNode, desc(toNode))
unique(mod56.edges.new$toNode)

net<- network(mod56.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod56.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77'))
mod56.final

ggsave("module_56_wle_upland.TIFF", plot=mod56.final, dpi=300, height=5, width=6, units=c("in"))


##module 67


##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero/module specific nodes and edges")

mod67.nodes<- read.csv("Cytoscape_Nodes_ME67_0-thresh_signed.txt", sep="")

mod67.edges<- read.csv("Cytoscape_Edges_ME67_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod67.edges,
                         vertex.attr = mod67.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod67.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE upland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod67.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c("#560d0d","#636bb7","#111b77")
# "#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
# "#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
#  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
#   "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
# "darkslategrey", "hotpink3", "indianred1")

#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_67_wle_upland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_67_wle_upland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod67.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod67.edges<- mod67.edges %>% mutate(y=toNode)
mod67.edges.edit<- mod67.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod67.edges.new<- arrange(mod67.edges.edit, toNode, desc(toNode))
unique(mod67.edges.new$toNode)

net<- network(mod67.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod67.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7', 'condensed aromatic'='#521899'))
mod67.final

ggsave("module_67_wle_upland.TIFF", plot=mod67.final, dpi=300, height=5, width=6, units=c("in"))














