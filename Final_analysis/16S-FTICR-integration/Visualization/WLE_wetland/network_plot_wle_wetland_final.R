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
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE wetland 100 perm/TOM threshold zero")
mod.stats<- read.csv("network_stats_correlation_by_module_WLE_wetland.csv") ##wle is not case sensitive
mod.stats.filter<- mod.stats %>% filter(Module!="ME0")


##EXAMINE THE MOD.FILTER FILE AND SELECT MODULES BASED ON R SQUARE VALUES 
##picking modules 1,2,3,4,5 because these modules have ASVs in them. Interestingly no organic matter correlated with these modules.

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE wetland 100 perm/TOM threshold zero/module specific nodes and edges")

mod1.nodes<- read.csv("Cytoscape_Nodes_ME1_0-thresh_signed.txt", sep="")

mod1.edges<- read.csv("Cytoscape_Edges_ME1_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod1.edges,
                         vertex.attr = mod1.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod1.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE wetland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod1.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c('#924900','#599861','#b66dff','#CD9BCD','#920000','#ff6db6','#648fff','#ffb6db',"goldenrod",'#b6dbff','#004949','#560d0d','#009292')
 #"#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
#"#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
 # "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
 #  "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
# "darkslategrey", "hotpink3", "indianred1")

###COLORS BY ORDER


##'o__Chthoniobacterales'='#1e0047'
##'o__Corynebacteriales'='#a0fffc'
#'o__Streptomycetales'='#D14285'
#'o__Rhizobiales'='#560d0d'
#'o__KD4-96'='forestgreen'
# 'o__Micrococcales'='#a35151'
#'o__Clostridiales'='darkslategrey'
#'o__Frankiales'='#fcb067'
#'o__B2M28'='#fa7efc'
#o__Pirellulales'='indianred1'
#"o__Desulfobacterales"="#599861"
#'o__KD4-96'='#CBD588'##need to keep this either forestgreen or this color
#'o__Anaerolineales'='#599861'##need to change this color
#'o__Run-SP154'='pink'
#'o__Ignavibacteriales'='mediumorchid1'
#'o__Solirubrobacterales'='#CD9BCD'
#'o__Bacillales'='#000000'
#'o__bacteriap25'='#004949'
#'o__Bacteroidales'='#009292'
#'o__Burkholderiales'='#ff6db6'
#'o__Desulfobaccales'='#ffb6db'
#'o__Gaiellales'='#490092'
#'o__Geobacterales'='#006edb'
#'o__MB-A2-108'='#b66dff'
#'o__Methylococcales'='#6db6ff'
#'o__Nannocystales'='#b6dbff'
#'o__Pedosphaerales'='#920000'
#'o__PeM15'='#924900'
#'o__Peptostreptococcales-Tissierellales'='#db6d00'
#'o__Subgroup_17'='#24ff24'
#'o__Syntrophales'='#ffff6d'
#'o__Vicinamibacterales'='#648fff'

#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_1_wle_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_1_wle_wetland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod1.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod1.edges<- mod1.edges %>% mutate(y=toNode)
mod1.edges.edit<- mod1.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod1.edges.new<- arrange(mod1.edges.edit, toNode, desc(toNode))
unique(mod1.edges.new$toNode)

net<- network(mod1.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod1.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__PeM15'='#924900','o__Anaerolineales'='#599861','o__MB-A2-108'='#b66dff',          
                                                                              'o__Solirubrobacterales'='#CD9BCD','o__Pedosphaerales'='#920000','o__Burkholderiales'='#ff6db6',   
                                                                              'o__Vicinamibacterales'='#648fff','o__Desulfobaccales'='#ffb6db','o__Desulfobacterales'='goldenrod',
                                                                              'o__Nannocystales'='#b6dbff','o__bacteriap25'='#004949','o__Rhizobiales'='#560d0d',        
                                                                              'o__Bacteroidales'='#009292'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod1.final

ggsave("module_1_wle_wetland_weighted.TIFF", plot=mod1.final, dpi=300, height=5, width=6, units=c("in"))

##go to module specific nodes and edges folder and select nodes and edges from the desired modules
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE wetland 100 perm/TOM threshold zero/module specific nodes and edges")

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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE wetland 100 perm/TOM threshold zero")
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

palette14= c('#599861','#CBD588','#CD9BCD','darkslategrey','#db6d00','#000000','#1e0047','goldenrod','#560d0d','mediumorchid1')
# "#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
#"#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
#  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
#   "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
# "darkslategrey", "hotpink3", "indianred1")


###COLORS BY ORDER


##'o__Chthoniobacterales'='#1e0047'
##'o__Corynebacteriales'='#a0fffc'
#'o__Streptomycetales'='#D14285'
#'o__Rhizobiales'='#560d0d'
#'o__KD4-96'='forestgreen'
# 'o__Micrococcales'='#a35151'
#'o__Clostridiales'='darkslategrey'
#'o__Frankiales'='#fcb067'
#'o__B2M28'='#fa7efc'
#o__Pirellulales'='indianred1'
#"o__Desulfobacterales"="goldenrod"
#'o__KD4-96'='#CBD588'##need to keep this either forestgreen or this color
#'o__Anaerolineales'='#599861'##need to change this color
#'o__Run-SP154'='pink'
#'o__Ignavibacteriales'='mediumorchid1'
#'o__Solirubrobacterales'='#CD9BCD'
#'o__Bacillales'='#000000'
#'o__bacteriap25'='#004949'
#'o__Bacteroidales'='#009292'
#'o__Burkholderiales'='#ff6db6'
#'o__Desulfobaccales'='#ffb6db'
#'o__Gaiellales'='#490092'
#'o__Geobacterales'='#006edb'
#'o__MB-A2-108'='#b66dff'
#'o__Methylococcales'='#6db6ff'
#'o__Nannocystales'='#b6dbff'
#'o__Pedosphaerales'='#920000'
#'o__PeM15'='#924900'
#'o__Peptostreptococcales-Tissierellales'='#db6d00'
#'o__Subgroup_17'='#24ff24'
#'o__Syntrophales'='#ffff6d'
#'o__Vicinamibacterales'='#648fff'








#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_2_wle_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_2_wle_wetland



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

mod2.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Anaerolineales'='#599861','o__KD4-96'='#CBD588'  ,                           
                                                                              'o__Solirubrobacterales'='#CD9BCD','o__Clostridiales'='darkslategrey' ,              
                                                                              'o__Peptostreptococcales-Tissierellales'='#db6d00','o__Bacillales'='#000000',                        
                                                                              'o__Chthoniobacterales'='#1e0047'  ,               'o__Desulfobacterales'='goldenrod',               
                                                                              'o__Rhizobiales'='#560d0d',                         'o__Ignavibacteriales'='mediumorchid1'       ), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod2.final

ggsave("module_2_wle_wetland_weighted.TIFF", plot=mod2.final, dpi=300, height=5, width=6, units=c("in"))

#################
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE wetland 100 perm/TOM threshold zero/module specific nodes and edges")

mod3.nodes<- read.csv("Cytoscape_Nodes_ME3_0-thresh_signed.txt", sep="")

mod3.edges<- read.csv("Cytoscape_Edges_ME3_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod3.edges,
                         vertex.attr = mod3.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE wetland 100 perm/TOM threshold zero")
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

palette14= c('#599861','#490092','#CD9BCD','darkslategrey','#920000','#ff6db6','#560d0d')
# "#283dff","#bfc5ff","#014443","#195637","#117744","#60ffaf","#b7ffdb","#825121",
#"#a35151", "#dba4a4", "pink",
#"#ea7f17","#fcb067","#ffe8d3","#d8d6d4","#82807f", "#3f3e3d","#000000",
#"#5b5b19","#ae09ea","#521899","#1e0047","mediumorchid1", "#fcfc00",
#  "#CBD588", "#599861", "#508578","#FFD700","#FFA500","#a0fffc",
#"#DA5724","#fff0d6","#8A7C64", "#AD6F3B", "#CD9BCD","#fa7efc",
#   "#D14285", "darkorchid1", "deepskyblue4", "forestgreen","goldenrod", "khaki4",
# "darkslategrey", "hotpink3", "indianred1")


###COLORS BY ORDER


#'o__Chthoniobacterales'='#1e0047'
##'o__Corynebacteriales'='#a0fffc'
#'o__Streptomycetales'='#D14285'
#'o__Rhizobiales'='#560d0d'
#'o__KD4-96'='forestgreen'
# 'o__Micrococcales'='#a35151'
#'o__Clostridiales'='darkslategrey'
#'o__Frankiales'='#fcb067'
#'o__B2M28'='#fa7efc'
#o__Pirellulales'='indianred1'
#"o__Desulfobacterales"="goldenrod"
#'o__KD4-96'='#CBD588'##need to keep this either forestgreen or this color
#'o__Anaerolineales'='#599861'##need to change this color for desulfobacterales
#'o__Run-SP154'='pink'
#'o__Ignavibacteriales'='mediumorchid1'
#'o__Solirubrobacterales'='#CD9BCD'
#'o__Bacillales'='#000000'
#'o__bacteriap25'='#004949'
#'o__Bacteroidales'='#009292'
#'o__Burkholderiales'='#ff6db6'
#'o__Desulfobaccales'='#ffb6db'
#'o__Gaiellales'='#490092'
#'o__Geobacterales'='#006edb'
#'o__MB-A2-108'='#b66dff'
#'o__Methylococcales'='#6db6ff'
#'o__Nannocystales'='#b6dbff'
#'o__Pedosphaerales'='#920000'
#'o__PeM15'='#924900'
#'o__Peptostreptococcales-Tissierellales'='#db6d00'
#'o__Subgroup_17'='#24ff24'
#'o__Syntrophales'='#ffff6d'
#'o__Vicinamibacterales'='#648fff'





#names(net)<- levels(palette)
library(RColorBrewer)
library(ggplot2)
library(grid)
library(patchwork)

module_3_wle_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_3_wle_wetland



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

mod3.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Anaerolineales'='#599861' ,     'o__Gaiellales'='#490092'   ,       'o__Solirubrobacterales'='#CD9BCD' ,'o__Clostridiales'='darkslategrey',
                                                                              'o__Pedosphaerales'='#920000'    ,  'o__Burkholderiales'='#ff6db6'  ,   'o__Rhizobiales'='#560d0d' ))
mod3.final

ggsave("module_3_wle_wetland.TIFF", plot=mod3.final, dpi=300, height=5, width=6, units=c("in"))


#########
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE wetland 100 perm/TOM threshold zero/module specific nodes and edges")

mod4.nodes<- read.csv("Cytoscape_Nodes_ME4_0-thresh_signed.txt", sep="")

mod4.edges<- read.csv("Cytoscape_Edges_ME4_0-thresh_signed.txt", sep="")


##network plot option 1
network_basic <- network(mod4.edges,
                         vertex.attr = mod4.nodes,
                         matrix.type = "edgelist",
                         ignore.eval = FALSE)
plot(network_basic)


####
##trying with ggnet

library(viridis)
pal=viridis(n=100,option="D")

##network plot option 2
net<- network(mod4.edges, directed=FALSE)##changing directed to TRUE
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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE wetland 100 perm/TOM threshold zero")
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
write.csv(merge.metadata.otu.icr, "merge.metadata.otu.icr.order.mod4.csv")

new_var<- merge.metadata.otu.icr$value
net %v% "Order/Class of feature" = new_var

#palette options

library(randomcoloR)
n <- 129
palette <- distinctColorPalette(n)

palette14= c('#6db6ff','#006edb','#ffff6d','#560d0d','mediumorchid1')
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

module_4_wle_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_4_wle_wetland



####setting colors a different way by naming the color for each level, this works better

merge.metadata.otu.icr<- read.csv("merge.metadata.otu.icr.order.mod4.csv")#remove empty taxa order rows if needed
merge.metadata.otu.icr.unique<- as.data.frame(unique(merge.metadata.otu.icr$value))
palette.data<- as.data.frame(palette14)
merge.otu.icr.palette<- cbind(merge.metadata.otu.icr.unique, palette.data)
merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(value=`unique(merge.metadata.otu.icr$value)`)
join.otu.icr.palette<- left_join(merge.metadata.otu.icr, merge.otu.icr.palette, by="value")

merge.otu.icr.palette<- merge.otu.icr.palette %>% mutate(palette_new=paste0("'", value, "'", "=","'",palette14, "'"))
merge.otu.icr.palette$palette_new

##matching the order of vertex names in edges_select file with join.otu.icr.palette
mod4.edges<- mod4.edges %>% mutate(y=toNode)
mod4.edges.edit<- mod4.edges %>% right_join(join.otu.icr.palette, by="y")  %>% filter(!is.na(toNode)==TRUE)


join.otu.icr.palette.new<- arrange(join.otu.icr.palette, y, desc(y))
mod4.edges.new<- arrange(mod4.edges.edit, toNode, desc(toNode))
unique(mod4.edges.new$toNode)

net<- network(mod4.edges.new, directed=FALSE)##changing directed to TRUE
z=network.vertex.names(net)
z.data<- as.data.frame(z)
join.otu.icr.palette.new.1 <- join.otu.icr.palette.new[match(z.data$z, join.otu.icr.palette.new$y),]

new_var<- join.otu.icr.palette.new.1$value
net %v% "Class" = new_var
#col<- join.otu.icr.palette$palette
#net %v% "Class_col" <- col
type<- join.otu.icr.palette.new.1$x
net %v% "type" <- type

mod4.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Methylococcales'='#6db6ff'    ,     'o__Geobacterales'='#006edb',           'o__Syntrophales'='#ffff6d' ,          
                                                                              'o__Rhizobiales'='#560d0d'    ,         'o__Ignavibacteriales'='mediumorchid1'))
mod4.final

ggsave("module_4_wle_wetland.TIFF", plot=mod4.final, dpi=300, height=5, width=6, units=c("in"))


##############
setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE wetland 100 perm/TOM threshold zero/module specific nodes and edges")

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

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE wetland 100 perm/TOM threshold zero")
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

palette14= c('#490092','#ff6db6','#648fff','#24ff24','#006edb','#560d0d')
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

module_5_wle_wetland<- ggnet2(net, node.size =3, color = "Order/Class of feature") #palette="Class_col") #color most likely here is in grayscale
module_5_wle_wetland



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

mod5.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Gaiellales'='#490092' ,'o__Burkholderiales'='#ff6db6','o__Vicinamibacterales'='#648fff','o__Subgroup_17'='#24ff24'  ,     
                                                                              'o__Geobacterales'='#006edb','o__Rhizobiales'='#560d0d'))
mod5.final

ggsave("module_5_wle_wetland.TIFF", plot=mod5.final, dpi=300, height=5, width=6, units=c("in"))


##matching the wetland indicators with the otus found in wetland modules

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Cytoscape Files signed WLE wetland 100 perm/TOM threshold zero")
mod.stats<- read.csv("network_stats_correlation_by_module_WLE_wetland.csv") ##wle is not case sensitive
mod.stats.filter<- mod.stats %>% filter(Module!="ME0")
View(mod.stats.filter)

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Indicator species and modules")
ind_wle<- read.csv("indicator_transect_WLE.csv", header=TRUE)
ind_wle_edit<- ind_wle %>% mutate(query=OTU)
ind_mod<- mod.stats.filter %>% left_join(ind_wle_edit, by="query")

setwd("/Users/band985/Library/CloudStorage/OneDrive-PNNL/Documents/16S synoptic site/16S-FTICR-combined-analysis/Analysis after first meeting/Indicator species and modules/wle wetland")
write.csv(ind_mod,"indicator_module_wle_wetland.csv")





