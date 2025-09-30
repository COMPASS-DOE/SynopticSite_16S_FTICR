arrow_map <- aes(xend = xvar, yend = yvar, x = 0, y = 0, 
                 color = NULL)

label_map <- aes(x = xvar, y = yvar, ##0.02 IN CB, 0.05 IN WLE
                 color = NULL, label = varname)

arrowhead = arrow(length = unit(0.02, "npc"), type = "closed")


##new plot FINAL
biplot_wle<-ggbiplot(pca_wle$pca_int, obs.scale = 1, var.scale = 1,
                     #groups = pca_wle$grp$transect, 
                     ellipse = FALSE, circle = FALSE, var.axes = FALSE, varname.size = 5, varname.adjust = 2.2)+
  #geom_segment(mapping = arrow_map, linewidth = 1.2, data = arrowdf, color = "black", arrow = arrowhead) + 
  #geom_segment(mapping = arrow_map, linewidth = 0.5, data = arrowdf, color = "white", arrow = arrowhead) + 
  #geom_label_repel(mapping = label_map, data = filter(arrowdf, xend < 0), show.legend = FALSE, size=7*6/14, hjust=0.5, fill="white", color="black") + 
  #geom_label_repel(mapping = label_map, data = filter(arrowdf, xend > 0), show.legend = FALSE, size=7*6/14, hjust=0.6, fill="white", color="black") + 
  #title = "Erie") + 
  scale_fill_manual(values = c("#CC79A7", "#E69F00", "#0072B2"))+
  scale_shape_manual(values=c(21,22,24))+
  geom_point(aes(shape = pca_wle$grp$site, fill = pca_wle$grp$transect), size = 6, color="black") +
  #geom_point(colour = "white", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))+labs(color="transect")+theme_classic()
  ggtitle("Erie")+
  theme_classic()+
  xlim(-4,4)+
  ylim(-4,4)+
  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14))+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))+
  theme(legend.text = element_text(size=12),  legend.title = element_text(size=16))+
  labs(fill="transect", shape="site")+
  guides(fill=guide_legend(override.aes = list(shape=21), order=1),
         shape = guide_legend(override.aes = list(fill = "black"), order=2))+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16 ))


biplot_wle

#getting the arrowdf from prcomp loadings

# Extracting variable loadings for arrows
loadings <- pca_wle$pca_int$rotation
loadings

## Create the arrowdf
## Assume you want to scale the arrows properly; adjust the scaling factor accordingly
scaling_factor <- 3 # This could be tuned
arrowdf <- as.data.frame(loadings) * scaling_factor
arrowdf$varname <- rownames(arrowdf) # Add variable names

# Naming the columns for ggplot compatibility
colnames(arrowdf) <- c("xvar", "yvar", "zvar", "wvar", "varname")

# Remove unused dimensions (if only the first two PCs are considered)
arrowdf <- arrowdf %>% select(xvar, yvar, varname) 

# Verifying arrowdf content
print(arrowdf)


##code to bring arrows and text to the front

# Add arrow and text layers on top
biplot_wle1 <- biplot_wle +
  geom_segment(mapping = arrow_map, linewidth = 1.5, data = arrowdf, color = "black", arrow = arrowhead) +
 # geom_segment(mapping = arrow_map, linewidth = 0.5, data = arrowdf, color = "white", arrow = arrowhead) +
  geom_label_repel(mapping = label_map, data = filter(arrowdf, xvar < 0), 
                   show.legend = FALSE, size = 11 * 6 / 14, nudge_x = ifelse(arrowdf$xvar < 0, -0.15, 0.15), 
                   nudge_y = ifelse(arrowdf$yvar < 0, 0.8, -0.8), fill = "white", color = "black") +
  geom_label_repel(mapping = label_map, data = filter(arrowdf, xvar > 0), 
                   show.legend = FALSE, size = 11 * 6 / 14, nudge_x = ifelse(arrowdf$xvar < 0, -0.15, 0.15), 
                   nudge_y = ifelse(arrowdf$yvar < 0, -0.6, 0.15), fill = "white", color = "black")

# Print the final plot
print(biplot_wle1)

ggsave(file="biplot_wle_pca_newlabels_04.04.25.tiff", plot = biplot_wle1, width = 7, height = 6)

#######################
####scrap code not used
######################

###################
###start scrap code
###################

# Optional: Modifying the arrow thickness if needed if using ggbiplot arrows
seg <- which(sapply(biplot_wle$layers, function(x) class(x$geom)[1] == 'GeomSegment'))
biplot_wle$layers[[seg[[1]]]]$aes_params$colour <- 'black'
biplot_wle$layers[[seg[[1]]]]$aes_params$size <- 1.5


##the below code is using layers within ggbiplot, i did not use this because i could not control the plot aesthetics well
seg <- which(sapply(biplot_wle$layers, function(x) class(x$geom)[1] == 'GeomSegment'))
txt <- which(sapply(biplot_wle$layers, function(x) class(x$geom)[1] == 'GeomText'))

#fixing arrow width
biplot_wle$layers[[seg]]$aes_params$colour <- 'black'
biplot_wle$layers[[seg]]$aes_params$size <- 1.5

biplot_wle

##making plot without text
biplot_wle$layers[[txt]]$aes_params$colour <- 'white'

ggsave(file="biplot_wle_pca_nolabel.tiff", plot = biplot_wle, width = 7, height = 6)

#geom_segment(mapping = arrow_map, linewidth = 1.2, data = arrowdf, color = "black", arrow = arrowhead) 
#geom_segment(mapping = arrow_map, linewidth = 0.5, data = arrowdf, color = "white", arrow = arrowhead) 
#geom_label_repel(mapping = label_map, data = filter(biplot_wle$layers[[txt]]$data, xend < 0), show.legend = FALSE, size=7*6/14, hjust=0.5, fill="white", color="black") 
#geom_label_repel(mapping = label_map, data = filter(biplot_wle$layers[[txt]]$data, xend > 0), show.legend = FALSE, size=7*6/14, hjust=0.6, fill="white", color="black") 

library(ggrepel)
#biplot_wle$layers[[txt]] <- geom_label_repel(mapping = label_map, data = filter(biplot_wle$layers[[txt]]$data, xvar < 0), show.legend = FALSE, size=6, hjust=0.9, vjust=1, fill="white", color="black", force = 2)  
#biplot_wle$layers[[txt]] <- geom_label_repel(mapping = label_map, data = filter(biplot_wle$layers[[txt]]$data, xvar > 0), show.legend = FALSE, size=6, hjust=1, vjust=1, fill="white", color="black") 

#fixing text parameters
biplot_wle$layers[[txt]] <- geom_label_repel(mapping = label_map, data = filter(biplot_wle$layers[[txt]]$data), show.legend = FALSE, size=6, hjust=1, vjust=0.48,  fill="white", color="black") 
biplot_wle

#fixing arrow parameters
biplot_wle$layers[[seg]] <- geom_segment(mapping = arrow_map, linewidth = 1.2, data = filter(biplot_wle$layers[[seg]]$data), color = "black", arrow = arrowhead) 

# biplot_wle$layers[[seg]] <- geom_segment(mapping = arrow_map, linewidth = 0.5, data = filter(biplot_wle$layers[[seg]]$data), fill = "white", arrow = arrowhead) 


#biplot_wle$layers[[txt]] <- geom_label(aes(x = xvar, y = yvar, label = varname,
                                      #     angle = angle, hjust = hjust), 
                                   #    label.size = NA,
                                    #   data = biplot_wle$layers[[txt]]$data, 
                                    #   fill = '#dddddd80')

biplot_wle

##################
####end scrap code
#################

##chesapeake 


arrow_map <- aes(xend = xvar, yend = yvar, x = 0, y = 0, 
                 color = NULL)

label_map <- aes(x = xvar, y = yvar, ##0.02 IN CB, 0.05 IN WLE
                 color = NULL, label = varname)

arrowhead = arrow(length = unit(0.02, "npc"), type = "closed")

##new plot FINAL
biplot_cb<-ggbiplot(pca_cb$pca_int, obs.scale = 1, var.scale = 1,
                    #groups = pca_cb$grp$transect, 
                    ellipse = FALSE, circle = FALSE, varname.size = 5, varname.adjust = 2.2, var.axes = FALSE) +
  #title = "Chesapeake") + 
  scale_fill_manual(values = c("#CC79A7", "#E69F00", "#0072B2"))+
  scale_shape_manual(values=c(21,22,24))+
  geom_point(aes(shape = pca_cb$grp$site, fill = pca_cb$grp$transect), size = 6, color="black") +
  #geom_point(colour = "white", size = 1.5) + theme(plot.title = element_text(hjust = 0.5))+labs(color="transect")+theme_classic()
  ggtitle("Chesapeake")+
  theme_classic()+
  xlim(-4,4)+
  ylim(-4,4)+
  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14))+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))+
  theme(legend.text = element_text(size=12),  legend.title = element_text(size=16))+
  labs(fill="transect", shape="site")+
  guides(fill=guide_legend(override.aes = list(shape=21), order=1),
         shape = guide_legend(override.aes = list(fill = "black"), order=2))+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16 ))

biplot_cb

#getting the arrowdf from prcomp loadings

# Extracting variable loadings for arrows
loadings <- pca_cb$pca_int$rotation
loadings

## Create the arrowdf
## Assume you want to scale the arrows properly; adjust the scaling factor accordingly
scaling_factor <- 3 # This could be tuned
arrowdf <- as.data.frame(loadings) * scaling_factor
arrowdf$varname <- rownames(arrowdf) # Add variable names

# Naming the columns for ggplot compatibility
colnames(arrowdf) <- c("xvar", "yvar", "zvar", "wvar", "varname")

# Remove unused dimensions (if only the first two PCs are considered)
arrowdf <- arrowdf %>% select(xvar, yvar, varname) 

# Verifying arrowdf content
print(arrowdf)


##code to bring arrows and text to the front

# Add arrow and text layers on top
biplot_cb1 <- biplot_cb +
  geom_segment(mapping = arrow_map, linewidth = 1.5, data = arrowdf, color = "black", arrow = arrowhead) +
  # geom_segment(mapping = arrow_map, linewidth = 0.5, data = arrowdf, color = "white", arrow = arrowhead) +
  geom_label_repel(mapping = label_map, data = filter(arrowdf, xvar < 0), 
                   show.legend = FALSE, size = 11 * 6 / 14, nudge_x = ifelse(arrowdf$xvar < 0, -0.15, 0.15), 
                   nudge_y = ifelse(arrowdf$yvar < 0, 0.7, -0.2), fill = "white", color = "black") +
  geom_label_repel(mapping = label_map, data = filter(arrowdf, xvar > 0), 
                   show.legend = FALSE, size = 11 * 6 / 14, nudge_x = ifelse(arrowdf$xvar < 0, -0.15, 0.15), 
                   nudge_y = ifelse(arrowdf$yvar < 0, -0.8, 0.7), fill = "white", color = "black")

# Print the final plot
print(biplot_cb1)

ggsave(file="biplot_cb_pca_newlabels_04.04.25.tiff", plot = biplot_cb1, width = 7, height = 6)

########start scrap code
##the below code is using layers within ggbiplot, i did not use this because i could not control the plot aesthetics well

seg <- which(sapply(biplot_cb$layers, function(x) class(x$geom)[1] == 'GeomSegment'))
txt <- which(sapply(biplot_cb$layers, function(x) class(x$geom)[1] == 'GeomText'))

#fixing arrow width
biplot_cb$layers[[seg]]$aes_params$colour <- 'black'
biplot_cb$layers[[seg]]$aes_params$size <- 1.5
biplot_cb
##making plot without text
biplot_cb$layers[[txt]]$aes_params$colour <- 'white'
biplot_cb
ggsave(file="biplot_cb_pca_nolabel.tiff", plot = biplot_cb, width = 7, height = 6)

#geom_segment(mapping = arrow_map, linewidth = 1.2, data = arrowdf, color = "black", arrow = arrowhead) 
#geom_segment(mapping = arrow_map, linewidth = 0.5, data = arrowdf, color = "white", arrow = arrowhead) 
#geom_label_repel(mapping = label_map, data = filter(biplot_wle$layers[[txt]]$data, xend < 0), show.legend = FALSE, size=7*6/14, hjust=0.5, fill="white", color="black") 
#geom_label_repel(mapping = label_map, data = filter(biplot_wle$layers[[txt]]$data, xend > 0), show.legend = FALSE, size=7*6/14, hjust=0.6, fill="white", color="black") 

library(ggrepel)
#biplot_wle$layers[[txt]] <- geom_label_repel(mapping = label_map, data = filter(biplot_wle$layers[[txt]]$data, xvar < 0), show.legend = FALSE, size=6, hjust=0.9, vjust=1, fill="white", color="black", force = 2)  
#biplot_wle$layers[[txt]] <- geom_label_repel(mapping = label_map, data = filter(biplot_wle$layers[[txt]]$data, xvar > 0), show.legend = FALSE, size=6, hjust=1, vjust=1, fill="white", color="black") 

#fixing text parameters
biplot_wle$layers[[txt]] <- geom_label_repel(mapping = label_map, data = filter(biplot_wle$layers[[txt]]$data), show.legend = FALSE, size=6, hjust=1, vjust=0.48,  fill="white", color="black") 
biplot_wle

#fixing arrow parameters
biplot_wle$layers[[seg]] <- geom_segment(mapping = arrow_map, linewidth = 1.2, data = filter(biplot_wle$layers[[seg]]$data), color = "black", arrow = arrowhead) 

# biplot_wle$layers[[seg]] <- geom_segment(mapping = arrow_map, linewidth = 0.5, data = filter(biplot_wle$layers[[seg]]$data), fill = "white", arrow = arrowhead) 


#biplot_wle$layers[[txt]] <- geom_label(aes(x = xvar, y = yvar, label = varname,
#     angle = angle, hjust = hjust), 
#    label.size = NA,
#   data = biplot_wle$layers[[txt]]$data, 
#   fill = '#dddddd80')

biplot_wle

##end scrap code
