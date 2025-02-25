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
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero")
mod.stats<- read.csv("network_stats_correlation_by_module_CB_wetland.csv")
mod.stats.filter<- mod.stats %>% filter(Module!="ME0")


##EXAMINE THE MOD.FILTER FILE AND SELECT MODULES BASED ON R SQUARE VALUES 
##picking modules, 2, 3, 5, 6, 7, 8, 11, 17, 18, 26, 30 because these modules have ASVs in them along with organic matter

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero/module specific nodes and edges")

mod2.nodes<- read.csv("Cytoscape_Nodes_ME2_0-thresh_signed.txt", sep="")

mod2.edges<- read.csv("Cytoscape_Edges_ME2_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod2.edges,
                         vertex.attr = mod2.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod2.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod2.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c("#560d0d","#111b77","#636bb7", "#521899")
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

module_2_cb_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_2_cb_wetland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod2.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod2.edges<- mod2.edges %>% mutate(y=toNode)
mod2.edges.edit<- mod2.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod2.edges.new<- arrange(mod2.edges.edit, toNode, desc(toNode))
unique(mod2.edges.new$toNode)

net<- network(mod2.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod2.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7', 'condensed aromatic'='#521899'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod2.final

ggsave("module_2_cb_wetland_weighted.TIFF", plot=mod2.final, dpi=300, height=5, width=6, units=c("in"))


##module 3

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero/module specific nodes and edges")

mod3.nodes<- read.csv("Cytoscape_Nodes_ME3_0-thresh_signed.txt", sep="")

mod3.edges<- read.csv("Cytoscape_Edges_ME3_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod3.edges,
                         vertex.attr = mod3.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)



##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod3.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod3.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c("#fa7efc","#636bb7", "#521899", "#111b77")
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

module_3_cb_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_3_cb_wetland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod3.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod3.edges<- mod3.edges %>% mutate(y=toNode)
mod3.edges.edit<- mod3.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod3.edges.new<- arrange(mod3.edges.edit, toNode, desc(toNode))
unique(mod3.edges.new$toNode)

net<- network(mod3.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod3.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__B2M28'='#fa7efc', 'aromatic'='#636bb7','condensed aromatic'='#521899', 'unsaturated/lignin'='#111b77'))
mod3.final

ggsave("module_3_cb_wetland.TIFF", plot=mod3.final, dpi=300, height=5, width=6, units=c("in"))

##module 5
##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero/module specific nodes and edges")

mod5.nodes<- read.csv("Cytoscape_Nodes_ME5_0-thresh_signed.txt", sep="")

mod5.edges<- read.csv("Cytoscape_Edges_ME5_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod5.edges,
                         vertex.attr = mod5.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod5.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod5.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c("indianred1","#560d0d", "#111b77", "#636bb7")
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

module_5_cb_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_5_cb_wetland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod5.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod5.edges<- mod5.edges %>% mutate(y=toNode)
mod5.edges.edit<- mod5.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod5.edges.new<- arrange(mod5.edges.edit, toNode, desc(toNode))
unique(mod5.edges.new$toNode)

net<- network(mod5.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod5.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Pirellulales'='indianred1', 'o__Rhizobiales'='#560d0d', 'unsaturated/lignin'='#111b77' ,'aromatic'='#636bb7'))
mod5.final

ggsave("module_5_cb_wetland.TIFF", plot=mod5.final, dpi=300, height=5, width=6, units=c("in"))

##module 7

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero/module specific nodes and edges")

mod7.nodes<- read.csv("Cytoscape_Nodes_ME7_0-thresh_signed.txt", sep="")

mod7.edges<- read.csv("Cytoscape_Edges_ME7_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod7.edges,
                         vertex.attr = mod7.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod7.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod7.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c ("#599861","#111b77", "#636bb7","#521899")
# "#560d0d","#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
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

module_7_cb_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_7_cb_wetland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod7.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod7.edges<- mod7.edges %>% mutate(y=toNode)
mod7.edges.edit<- mod7.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod7.edges.new<- arrange(mod7.edges.edit, toNode, desc(toNode))
unique(mod7.edges.new$toNode)

net<- network(mod7.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod7.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c("o__Desulfobacterales"="#599861",'unsaturated/lignin'='#111b77','aromatic'='#636bb7', 'condensed aromatic'='#521899'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod7.final

ggsave("module_7_cb_wetland_weighted.TIFF", plot=mod7.final, dpi=300, height=5, width=6, units=c("in"))


##module 8

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero/module specific nodes and edges")

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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero")
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

palette14= c ("#599861","#111b77", "#636bb7","#521899")
# "#560d0d","#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
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

module_8_cb_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_8_cb_wetland



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

mod8.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c("o__Desulfobacterales"="#599861",'unsaturated/lignin'='#111b77','aromatic'='#636bb7', 'condensed aromatic'='#521899'))
mod8.final

ggsave("module_8_cb_wetland.TIFF", plot=mod8.final, dpi=300, height=5, width=6, units=c("in"))


##module 11

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero/module specific nodes and edges")

mod11.nodes<- read.csv("Cytoscape_Nodes_ME11_0-thresh_signed.txt", sep="")

mod11.edges<- read.csv("Cytoscape_Edges_ME11_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod11.edges,
                         vertex.attr = mod11.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod11.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod11.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c ("#CBD588", "#599861","#111b77", "#636bb7","#521899")
# "#560d0d","#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
 #"#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
#"#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
#   "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
# "darkslategrey", "hotpink3", "indianred1")

#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_11_cb_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_11_cb_wetland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod11.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod11.edges<- mod11.edges %>% mutate(y=toNode)
mod11.edges.edit<- mod11.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod11.edges.new<- arrange(mod11.edges.edit, toNode, desc(toNode))
unique(mod11.edges.new$toNode)

net<- network(mod11.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod11.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__KD4-96'='#CBD588','o__Anaerolineales'='#599861', 'unsaturated/lignin'='#111b77','aromatic'='#636bb7', 'condensed aromatic'='#521899'))
mod11.final

ggsave("module_11_cb_wetland.TIFF", plot=mod11.final, dpi=300, height=5, width=6, units=c("in"))


##module 17

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero/module specific nodes and edges")

mod17.nodes<- read.csv("Cytoscape_Nodes_ME17_0-thresh_signed.txt", sep="")

mod17.edges<- read.csv("Cytoscape_Edges_ME17_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod17.edges,
                         vertex.attr = mod17.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod17.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod17.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c ("#560d0d","#111b77","#521899", "#636bb7")
# "#560d0d","#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
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

module_17_cb_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_17_cb_wetland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod17.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod17.edges<- mod17.edges %>% mutate(y=toNode)
mod17.edges.edit<- mod17.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod17.edges.new<- arrange(mod17.edges.edit, toNode, desc(toNode))
unique(mod17.edges.new$toNode)

net<- network(mod17.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod17.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77','condensed aromatic'='#521899','aromatic'='#636bb7'))
mod17.final

ggsave("module_17_cb_wetland.TIFF", plot=mod17.final, dpi=300, height=5, width=6, units=c("in"))

##module 18

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero/module specific nodes and edges")

mod18.nodes<- read.csv("Cytoscape_Nodes_ME18_0-thresh_signed.txt", sep="")

mod18.edges<- read.csv("Cytoscape_Edges_ME18_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod18.edges,
                         vertex.attr = mod18.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod18.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod18.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c ("pink", "#636bb7", "#111b77")
# "#560d0d","#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
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

module_18_cb_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_18_cb_wetland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod18.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod18.edges<- mod18.edges %>% mutate(y=toNode)
mod18.edges.edit<- mod18.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod18.edges.new<- arrange(mod18.edges.edit, toNode, desc(toNode))
unique(mod18.edges.new$toNode)

net<- network(mod18.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod18.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Run-SP154'='pink', 'aromatic'='#636bb7', 'unsaturated/lignin'='#111b77'))
mod18.final

ggsave("module_18_cb_wetland.TIFF", plot=mod18.final, dpi=300, height=5, width=6, units=c("in"))


##module 26

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero/module specific nodes and edges")

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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero")
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

palette14= c ("mediumorchid1","#111b77")
# "#560d0d","#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
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

module_26_cb_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_26_cb_wetland



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

mod26.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Ignavibacteriales'='mediumorchid1','unsaturated/lignin'='#111b77'))
mod26.final

ggsave("module_26_cb_wetland.TIFF", plot=mod26.final, dpi=300, height=5, width=6, units=c("in"))



##module 30

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero/module specific nodes and edges")

mod30.nodes<- read.csv("Cytoscape_Nodes_ME30_0-thresh_signed.txt", sep="")

mod30.edges<- read.csv("Cytoscape_Edges_ME30_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod30.edges,
                         vertex.attr = mod30.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod30.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed CB wetland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod30.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c ("#CD9BCD","#111b77")
# "#560d0d","#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
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

module_30_cb_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_30_cb_wetland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod30.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod30.edges<- mod30.edges %>% mutate(y=toNode)
mod30.edges.edit<- mod30.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod30.edges.new<- arrange(mod30.edges.edit, toNode, desc(toNode))
unique(mod30.edges.new$toNode)

net<- network(mod30.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod30.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Solirubrobacterales'='#CD9BCD','unsaturated/lignin'='#111b77'))
mod30.final

ggsave("module_30_cb_wetland.TIFF", plot=mod30.final, dpi=300, height=5, width=6, units=c("in"))




