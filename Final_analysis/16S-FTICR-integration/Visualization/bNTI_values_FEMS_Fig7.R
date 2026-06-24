
library(tidyverse)
library(ggplot2)

######
#####cb upland
#####

cb_upland<- read.csv("cb_upland_summary_16s_icr.csv")

cb_upland$Color_Label <- factor(
  cb_upland$Color_Label, 
  levels = c('o__Frankiales', 'o__Rhizobiales','unsaturated/lignin','aromatic'), ordered = TRUE
)

bnti_plot_cb_upland<- ggplot(cb_upland, aes(x = Type, y = mean_value, color=Color_Label, shape=Type))+
  coord_cartesian(ylim = c(-20, 3)) +
  # Add the lines BEFORE your data points so they are in the background.
  # The color is SET, not mapped, so it will not affect the legend.
  
  # Red dashed lines at +2 and -2
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "red", size = 0.8) +
  
  # Blue dashed lines at +1 and -1
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue", size = 0.8) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "blue", size = 0.8) +
  geom_point(size=4, position = position_jitterdodge(dodge.width = 0.9))+
  #geom_point(size=4, position = position_jitter(width = 0.2, height=0))+
  scale_color_manual(values=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7','o__Frankiales'='#fcb067'),
                     guide = guide_legend(
                       title = "OTU Order/FTICR feature Class", 
                       override.aes = list(
                         shape = c(17,17,16,16),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
                         size=5)))+
  scale_shape_manual(values = c("OTU" = 17, "molecule" = 16), guide=guide_legend(title="Type", override.aes = list(size=5))) +
  labs(
    y = "βNTIfeat",
    shape="Type")+
  theme(axis.text.x = element_text(colour = "black", size = 12,angle = 0, vjust = 0.7, hjust=0.5),
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title.y= element_text(size=14, color="black"),
        axis.title.x= element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.title = element_text(size =12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")) +
  theme(plot.title=element_text(size=16, face="bold", hjust=0.5))
  
bnti_plot_cb_upland

ggsave("bnti_cb_upland_03.13.26_points_new.TIFF", plot=bnti_plot_cb_upland, dpi=300, height=4, width=5, units=c("in"))

ggsave("bnti_cb_upland_03.13.26_points_new.TIFF", plot=bnti_plot_cb_upland, dpi=300, height=4, width=7, units=c("in"))



######
###cb wetland
######

cb_wetland<- read.csv("cb_wetland_summary_16s_icr.csv")

cb_wetland$Color_Label <- factor(
  cb_wetland$Color_Label, 
  levels = c('o__Desulfobacterales', 'o__Rhizobiales','unsaturated/lignin','aromatic','condensed aromatic'), ordered = TRUE
)

bnti_plot_cb_wetland<- ggplot(cb_wetland, aes(x = Type, y = mean_value, color=Color_Label, shape=Type))+
  coord_cartesian(ylim = c(-20, 3)) +
  # Add the lines BEFORE your data points so they are in the background.
  # The color is SET, not mapped, so it will not affect the legend.
  
  # Red dashed lines at +2 and -2
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "red", size = 0.8) +
  
  # Blue dashed lines at +1 and -1
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue", size = 0.8) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "blue", size = 0.8) +
  #geom_point(size=4, position = position_jitter(width = 0.2, height=0))+
  geom_point(size=4, position = position_jitterdodge(dodge.width = 0.9))+
  scale_color_manual(values=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7', 'condensed aromatic'='#521899',"o__Desulfobacterales"="#599861"),
                     guide = guide_legend(
                       title = "OTU Order/FTICR feature Class", 
                       override.aes = list(
                         shape = c(17,17,16,16, 16),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
                         size=5)))+
  scale_shape_manual(values = c("OTU" = 17, "molecule" = 16), guide=guide_legend(title="Type", override.aes = list(size=5))) +
  labs(
    y = "βNTIfeat",
    shape="Type")+
  theme(axis.text.x = element_text(colour = "black", size = 12,angle = 0, vjust = 0.7, hjust=0.5),
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title.y= element_text(size=14, color="black"),
        axis.title.x= element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.title = element_text(size =12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")) +
  theme(plot.title=element_text(size=16, face="bold", hjust=0.5))

bnti_plot_cb_wetland

ggsave("bnti_cb_wetland_03.13.26_points_new.TIFF", plot=bnti_plot_cb_wetland, dpi=300, height=4, width=5, units=c("in"))

ggsave("bnti_cb_wetland_03.13.26_points_new.TIFF", plot=bnti_plot_cb_wetland, dpi=300, height=4, width=7, units=c("in"))

###
###wle upland
####


wle_upland<- read.csv("wle_upland_summary_16s_icr.csv")

wle_upland$Color_Label <- factor(
  wle_upland$Color_Label, 
  levels = c('o__Rhizobiales','o__Streptomycetales', 'unsaturated/lignin','aromatic'), ordered = TRUE
)

bnti_plot_wle_upland<- ggplot(wle_upland, aes(x = Type, y = mean_value, color=Color_Label, shape=Type))+
  coord_cartesian(ylim = c(-20, 3)) +
  # Add the lines BEFORE your data points so they are in the background.
  # The color is SET, not mapped, so it will not affect the legend.
  
  # Red dashed lines at +2 and -2
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "red", size = 0.8) +
  
  # Blue dashed lines at +1 and -1
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue", size = 0.8) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "blue", size = 0.8) +
 #geom_point(size=4, position = position_jitter(width = 0.2, height=0))+
  geom_point(size=4, position = position_jitterdodge(dodge.width = 0.9))+
  scale_color_manual(values=c('o__Rhizobiales'='#560d0d','o__Streptomycetales'='#D14285','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7'),
                     guide = guide_legend(
                       title = "OTU Order/FTICR feature Class", 
                       override.aes = list(
                         shape = c(17,17,16,16),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
                         size=5)))+
  scale_shape_manual(values = c("OTU" = 17, "molecule" = 16), guide=guide_legend(title="Type", override.aes = list(size=5))) +
  labs(
    y = "βNTIfeat",
    shape="Type")+
  theme(axis.text.x = element_text(colour = "black", size = 12,angle = 0, vjust = 0.7, hjust=0.5),
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title.y= element_text(size=14, color="black"),
        axis.title.x= element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.title = element_text(size =12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")) +
  theme(plot.title=element_text(size=16, face="bold", hjust=0.5))

bnti_plot_wle_upland

ggsave("bnti_wle_upland_03.13.26_points_new.TIFF", plot=bnti_plot_wle_upland, dpi=300, height=4, width=5, units=c("in"))

ggsave("bnti_wle_upland_03.13.26_points_new.TIFF", plot=bnti_plot_wle_upland, dpi=300, height=4, width=7, units=c("in"))


####
####wle transition
####


wle_transition<- read.csv("wle_transition_summary_16s_icr.csv")

wle_transition$Color_Label <- factor(
  wle_transition$Color_Label, 
  levels = c('o__Clostridiales','o__Rhizobiales', 'unsaturated/lignin','aromatic'), ordered = TRUE
)

bnti_plot_wle_transition<- ggplot(wle_transition, aes(x = Type, y = mean_value, color=Color_Label, shape=Type))+
  coord_cartesian(ylim = c(-20, 3)) +
  # Add the lines BEFORE your data points so they are in the background.
  # The color is SET, not mapped, so it will not affect the legend.
  
  # Red dashed lines at +2 and -2
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "red", size = 0.8) +
  
  # Blue dashed lines at +1 and -1
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue", size = 0.8) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "blue", size = 0.8) +
  #geom_point(size=4, position = position_jitter(width = 0.2, height=0))+
  geom_point(size=4, position = position_jitterdodge(dodge.width = 0.9))+
  scale_color_manual(values=c('o__Rhizobiales'='#560d0d','unsaturated/lignin'='#111b77', 'aromatic'='#636bb7','o__Clostridiales'='darkslategrey'),
                     guide = guide_legend(
                       title = "OTU Order/FTICR feature Class", 
                       override.aes = list(
                         shape = c(17,17,16,16),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
                         size=5)))+
  scale_shape_manual(values = c("OTU" = 17, "molecule" = 16), guide=guide_legend(title="Type", override.aes = list(size=5))) +
  labs(
    y = "βNTIfeat",
    shape="Type")+
  theme(axis.text.x = element_text(colour = "black", size = 12,angle = 0, vjust = 0.7, hjust=0.5),
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title.y= element_text(size=14, color="black"),
        axis.title.x= element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.title = element_text(size =12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")) +
  theme(plot.title=element_text(size=16, face="bold", hjust=0.5))

bnti_plot_wle_transition

ggsave("bnti_wle_transition_03.13.26_points_new.TIFF", plot=bnti_plot_wle_transition, dpi=300, height=4, width=5, units=c("in"))

ggsave("bnti_wle_transition_03.13.26_points_new.TIFF", plot=bnti_plot_wle_transition, dpi=300, height=4, width=7, units=c("in"))

###
###wle wetland
###

wle_wetland<- read.csv("wle_wetland_summary_16s.csv")

wle_wetland$Color_Label <- factor(
  wle_wetland$Color_Label, 
  levels = c('o__Anaerolineales','o__Bacillales', 'o__bacteriap25', 'o__Bacteroidales', 'o__Burkholderiales', 'o__Chthoniobacterales','o__Clostridiales', 'o__Desulfobaccales', 
             'o__Desulfobacterales', 'o__Ignavibacteriales','o__KD4-96', 'o__MB-A2-108', 'o__Nannocystales', 'o__Pedosphaerales', 'o__PeM15', 'o__Peptostreptococcales-Tissierellales',
             'o__Rhizobiales', 'o__Solirubrobacterales', 'o__Vicinamibacterales'), ordered = TRUE
)

bnti_plot_wle_wetland<- ggplot(wle_wetland, aes(x = Type, y = mean_value, color=Color_Label, shape=Type))+
  coord_cartesian(ylim = c(-20, 3)) +
  # Add the lines BEFORE your data points so they are in the background.
  # The color is SET, not mapped, so it will not affect the legend.
  
  # Red dashed lines at +2 and -2
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = -2, linetype = "dashed", color = "red", size = 0.8) +
  
  # Blue dashed lines at +1 and -1
  geom_hline(yintercept = 1, linetype = "dashed", color = "blue", size = 0.8) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "blue", size = 0.8) +
  #geom_point(size=4, position = position_jitter(width = 0.2, height=0))+
  geom_point(size=4, position = position_jitterdodge(dodge.width = 0.9))+
  scale_color_manual(values=c('o__Anaerolineales'='goldenrod', 'o__Bacillales'='#000000','o__bacteriap25'='darkorange', 'o__Bacteroidales'='#009292', 'o__Burkholderiales'='khaki', 
                             'o__Chthoniobacterales'='#1e0047'  , 'o__Clostridiales'='darkslategrey' , 'o__Desulfobaccales'='#ffb6db' , 'o__Desulfobacterales'='#599861', 
                             'o__Ignavibacteriales'='mediumorchid1', 'o__KD4-96'='#CBD588'  , 'o__MB-A2-108'='#b66dff','o__Nannocystales'='#b6dbff', 'o__Pedosphaerales'='#920000', 
                             'o__PeM15'='#924900', 'o__Peptostreptococcales-Tissierellales'='#db6d00','o__Rhizobiales'='#560d0d','o__Solirubrobacterales'='#CD9BCD',    
                             'o__Vicinamibacterales'='#648fff' ),
                     guide = guide_legend(
                       title = "OTU Order/FTICR feature Class", 
                       override.aes = list(
                         shape = c(17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17),# ensuring shapes match those in the plot; triangle for OTUs and circles for Features
                         size=5)) )+
  scale_shape_manual(values = c("OTU" = 17, "molecule" = 16), guide=guide_legend(title="Type", override.aes = list(size=5))) +
  labs(
    y = "βNTIfeat",
    shape="Type")+
  theme(axis.text.x = element_text(colour = "black", size = 12,angle = 0, vjust = 0.7, hjust=0.5),
        axis.text.y = element_text(size = 12, colour = "black"), 
        axis.title.y= element_text(size=14, color="black"),
        axis.title.x= element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"), 
        legend.title = element_text(size =12, colour = "black"),
        legend.text = element_text(size = 12, colour = "black")) +
  theme(plot.title=element_text(size=16, face="bold", hjust=0.5))

bnti_plot_wle_wetland

ggsave("bnti_wle_wetland_03.13.26_points_new.TIFF", plot=bnti_plot_wle_wetland, dpi=300, height=6, width=7, units=c("in"))


ggsave("bnti_wle_wetland_03.13.26_points_new.TIFF", plot=bnti_plot_wle_wetland, dpi=300, height=4, width=8, units=c("in"))


