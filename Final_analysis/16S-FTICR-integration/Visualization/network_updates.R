#CB upland module 2 edit 

mod2.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod2.final
color_mapping <- c('o__Rhizobiales'='#560d0d', 'unsaturated/lignin'='#111b77', 'aromatic'='#636bb7')
shape_mapping <- c('OTU' = 17, 'molecule' = 16) # 24 for triangle, 21 for circle
mod2.final <- mod2.final + 
    scale_shape_manual(values = shape_mapping, guide=guide_legend(title="Type", override.aes = list(size=7)))+
    scale_color_manual(
       values = color_mapping,
       breaks= c('o__Rhizobiales','unsaturated/lignin', 'aromatic'),
       guide = guide_legend(
       title = "OTU Order/ \nFTICR feature Class", 
       override.aes = list(
       shape = c(17,16,16),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
       size=7
          )
         )
    )

mod2.final
ggsave("module_2_cb_upland_weighted_legendupdate_04.04.25.TIFF", plot=mod2.final, dpi=300, height=5, width=6, units=c("in"))

##CB upland module 11 edit
mod11.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Frankiales'='#fcb067','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod11.final
color_mapping <- c('o__Frankiales'='#fcb067','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7')
shape_mapping <- c('OTU' = 17, 'molecule' = 16) # 24 for triangle, 21 for circle
mod11.final <- mod11.final + 
  scale_shape_manual(values = shape_mapping, guide=guide_legend(title="Type", override.aes = list(size=7)))+
  scale_color_manual(
    values = color_mapping,
    breaks= c('o__Frankiales','unsaturated/lignin', 'aromatic'),
    guide = guide_legend(
      title = "OTU Order/ \nFTICR feature Class", 
      override.aes = list(
        shape = c(17,16,16),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
        size=7
      )
    )
  )

mod11.final
ggsave("module_11_cb_upland_weighted_legendupdate_05.12.25.TIFF", plot=mod11.final, dpi=300, height=5, width=6, units=c("in"))


##cb wetland module 2 edit
mod2.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7', 'condensed aromatic'='#521899'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod2.final
color_mapping <- c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7', 'condensed aromatic'='#521899')
shape_mapping <- c('OTU' = 17, 'molecule' = 16) # 24 for triangle, 21 for circle
mod2.final <- mod2.final + 
  scale_shape_manual(values = shape_mapping, guide=guide_legend(title="Type", override.aes = list(size=7)))+
  scale_color_manual(
    values = color_mapping,
    breaks= c('o__Rhizobiales','unsaturated/lignin', 'aromatic', 'condensed aromatic'),
    guide = guide_legend(
      title = "OTU Order/ \nFTICR feature Class", 
      override.aes = list(
        shape = c(17,16,16,16),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
        size=7
      )
    )
  )
mod2.final

ggsave("module_2_cb_wetland_weighted_legendupdate_05.12.25.TIFF", plot=mod2.final, dpi=300, height=5, width=6, units=c("in"))





##cb wetland module 7 edit
mod7.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c("o__Desulfobacterales"="#599861",'unsaturated/lignin'='#111b77','aromatic'='#636bb7', 'condensed aromatic'='#521899'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod7.final

color_mapping <- c("o__Desulfobacterales"="#599861",'unsaturated/lignin'='#111b77','aromatic'='#636bb7', 'condensed aromatic'='#521899')
shape_mapping <- c('OTU' = 17, 'molecule' = 16) # 24 for triangle, 21 for circle
mod7.final <- mod7.final + 
  scale_shape_manual(values = shape_mapping, guide=guide_legend(title="Type", override.aes = list(size=7)))+
  scale_color_manual(
    values = color_mapping,
    breaks= c('o__Desulfobacterales','unsaturated/lignin','aromatic', 'condensed aromatic'),
    guide = guide_legend(
      title = "OTU Order/ \nFTICR feature Class", 
      override.aes = list(
        shape = c(17,16,16,16),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
        size=7
      )
    )
  )
mod7.final

ggsave("module_7_cb_wetland_weighted_legendupdate_04.04.25.TIFF", plot=mod7.final, dpi=300, height=5, width=6, units=c("in"))


##wle upland 

##module 33 edit
mod33.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod33.final

color_mapping <- c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7')
shape_mapping <- c('OTU' = 17, 'molecule' = 16) # 24 for triangle, 21 for circle
mod33.final <- mod33.final + 
  scale_shape_manual(values = shape_mapping, guide=guide_legend(title="Type", order=1, override.aes = list(size=7)))+
  scale_color_manual(
    values = color_mapping,
    breaks= c('o__Rhizobiales', 'unsaturated/lignin',  'aromatic'),
    guide = guide_legend(
      title = "OTU Order/ \nFTICR feature Class", 
      order=2,
      override.aes = list(
        shape = c(17,16,16),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
        size=7
      )
    )
  )
mod33.final

ggsave("module_33_wle_upland_weighted_legendupdate_05.12.25.TIFF", plot=mod33.final, dpi=300, height=5, width=6, units=c("in"))



#module 63 edit
mod63.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Streptomycetales'='#D14285', 'unsaturated/lignin'='#111b77',  'aromatic'='#636bb7'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod63.final

color_mapping <- c('o__Streptomycetales'='#D14285', 'unsaturated/lignin'='#111b77',  'aromatic'='#636bb7')
shape_mapping <- c('OTU' = 17, 'molecule' = 16) # 24 for triangle, 21 for circle
mod63.final <- mod63.final + 
  scale_shape_manual(values = shape_mapping, guide=guide_legend(title="Type", order=1, override.aes = list(size=7)))+
  scale_color_manual(
    values = color_mapping,
    breaks= c('o__Streptomycetales', 'unsaturated/lignin',  'aromatic'),
    guide = guide_legend(
      title = "OTU Order/ \nFTICR feature Class", 
      order=2,
      override.aes = list(
        shape = c(17,16,16),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
        size=7
      )
    )
  )
mod63.final


ggsave("module_63_wle_upland_weighted_legendupdate_04.04.25.TIFF", plot=mod63.final, dpi=300, height=5, width=6, units=c("in"))


##wle transition

##module 26 edit

mod26.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod26.final
color_mapping <- c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7')
shape_mapping <- c('OTU' = 17, 'molecule' = 16) # 24 for triangle, 21 for circle
mod26.final <- mod26.final + 
  scale_shape_manual(values = shape_mapping, guide=guide_legend(title="Type", order=1, override.aes = list(size=7)))+
  scale_color_manual(
    values = color_mapping,
    breaks= c('o__Rhizobiales','unsaturated/lignin', 'aromatic'),
    guide = guide_legend(
      title = "OTU Order/ \nFTICR feature Class", 
      order=2,
      override.aes = list(
        shape = c(17,16,16),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
        size=7
      )
    )
  )

mod26.final

ggsave("module_26_wle_transition_weighted_legendupdate_05.12.25.TIFF", plot=mod26.final, dpi=300, height=5, width=6, units=c("in"))

##module 8 edit
mod8.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Clostridiales'='darkslategrey','unsaturated/lignin'='#111b77'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod8.final

color_mapping <- c('o__Clostridiales'='darkslategrey','unsaturated/lignin'='#111b77')
shape_mapping <- c('OTU' = 17, 'molecule' = 16) # 24 for triangle, 21 for circle
mod8.final <- mod8.final + 
  scale_shape_manual(values = shape_mapping, guide=guide_legend(title="Type", order=1, override.aes = list(size=7)))+
  scale_color_manual(
    values = color_mapping,
    breaks= c('o__Clostridiales','unsaturated/lignin'),
    guide = guide_legend(
      title = "OTU Order/ \nFTICR feature Class", 
      order=2,
      override.aes = list(
        shape = c(17,16),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
        size=7
      )
    )
  )
mod8.final

ggsave("module_8_wle_transition_weighted_legendupdate_04.04.25.TIFF", plot=mod8.final, dpi=300, height=5, width=6, units=c("in"))



##wle wetland
##module 1 edit

mod1.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__PeM15'='#924900','o__Anaerolineales'='goldenrod','o__MB-A2-108'='#b66dff',          
                                                                              'o__Solirubrobacterales'='#CD9BCD','o__Pedosphaerales'='#920000','o__Burkholderiales'='khaki',   
                                                                              'o__Vicinamibacterales'='#648fff','o__Desulfobaccales'='#ffb6db','o__Desulfobacterales'='#599861',
                                                                              'o__Nannocystales'='#b6dbff','o__bacteriap25'='darkorange','o__Rhizobiales'='#560d0d',        
                                                                              'o__Bacteroidales'='#009292'), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
  #scale_shape_manual(values = 17,name="type")+guides(shape=guide_legend(override.aes = list(size = 6)))

mod1.final

color_mapping <- c('o__PeM15'='#924900','o__Anaerolineales'='goldenrod','o__MB-A2-108'='#b66dff',          
                   'o__Solirubrobacterales'='#CD9BCD','o__Pedosphaerales'='#920000','o__Burkholderiales'='khaki',   
                   'o__Vicinamibacterales'='#648fff','o__Desulfobaccales'='#ffb6db','o__Desulfobacterales'='#599861',
                   'o__Nannocystales'='#b6dbff','o__bacteriap25'='darkorange','o__Rhizobiales'='#560d0d',        
                   'o__Bacteroidales'='#009292')
shape_mapping <- c('OTU' = 17) # 24 for triangle, 21 for circle

##sort OTUs alphabetically
sorted_OTUs <- sort(names(color_mapping))

mod1.final <- mod1.final + 
  scale_shape_manual(values = shape_mapping, guide=guide_legend(title="Type", order=1, override.aes = list(size=7)))+
  scale_color_manual(
    values = color_mapping,
    breaks= sorted_OTUs,
    guide = guide_legend(
      title = "OTU Order/ \nFTICR feature Class", 
      order=2,
      override.aes = list(
        shape = rep(17, length(sorted_OTUs)),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
        size=7
      )
    )
  )
mod1.final

ggsave("module_1_wle_wetland_weighted_legendupdate_04.04.25.TIFF", plot=mod1.final, dpi=300, height=5, width=6, units=c("in"))


##module 2 edit
mod2.final<- ggnet2(net, node.size=7, color ="Class", shape="type", palette=c('o__Anaerolineales'='goldenrod','o__KD4-96'='#CBD588'  ,                           
                                                                              'o__Solirubrobacterales'='#CD9BCD','o__Clostridiales'='darkslategrey' ,              
                                                                              'o__Peptostreptococcales-Tissierellales'='#db6d00','o__Bacillales'='#000000',                        
                                                                              'o__Chthoniobacterales'='#1e0047'  ,               'o__Desulfobacterales'='#599861',               
                                                                              'o__Rhizobiales'='#560d0d',                         'o__Ignavibacteriales'='mediumorchid1'       ), edge.size="weight", legend.size=12, color.legend="OTU Order/ \nFTICR feature Class")
mod2.final

color_mapping <- c('o__Anaerolineales'='goldenrod','o__KD4-96'='#CBD588'  ,                           
                   'o__Solirubrobacterales'='#CD9BCD','o__Clostridiales'='darkslategrey' ,              
                   'o__Peptostreptococcales-Tissierellales'='#db6d00','o__Bacillales'='#000000',                        
                   'o__Chthoniobacterales'='#1e0047'  ,               'o__Desulfobacterales'='#599861',               
                   'o__Rhizobiales'='#560d0d',                         'o__Ignavibacteriales'='mediumorchid1'      )
shape_mapping <- c('OTU' = 17) # 24 for triangle, 21 for circle

##sort OTUs alphabetically
sorted_OTUs <- sort(names(color_mapping))

mod2.final <- mod2.final + 
  scale_shape_manual(values = shape_mapping, guide=guide_legend(title="Type", order=1, override.aes = list(size=7)))+
  scale_color_manual(
    values = color_mapping,
    breaks= sorted_OTUs,
    guide = guide_legend(
      title = "OTU Order/ \nFTICR feature Class", 
      order=2,
      override.aes = list(
        shape = rep(17, length(sorted_OTUs)),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
        size=7
      )
    )
  )
mod2.final

ggsave("module_2_wle_wetland_weighted_legendupdate_05.12.25.TIFF", plot=mod2.final, dpi=300, height=5, width=6, units=c("in"))



