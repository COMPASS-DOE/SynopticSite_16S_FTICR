library(tidyverse)
icr_data_long<- read.csv("icr_long_all_samples.csv")
icr_meta<- read.csv("icr_meta.csv")
sample_key<- read.csv("sample-metadata-cb-wle.csv")

##run this function to calculate icr relative abundances
compute_icr_relabund = function(icr_data_long, icr_meta){
  
  icr_data_long %>% 
    # filter(horizon != "B") %>% 
    # add the Class column to the data
    left_join(dplyr::select(icr_meta, formula, Class), by = "formula") %>% 
    # calculate abundance of each Class as the sum of all counts
    group_by(sample_label, Class) %>%
    dplyr::summarise(abund = sum(presence)) %>%
    ungroup %>% 
    # create a new column for total counts per core assignment
    # and then calculate relative abundance  
    group_by(sample_label) %>% 
    dplyr::mutate(total = sum(abund),
                  relabund  = ((abund/total)*100))
}


##run the function on data
icr_relabund<- compute_icr_relabund(icr_data_long, icr_meta)

## run functions below to make fit_pca_function_icr work
reorder_transect = function(dat){
  dat %>% 
    mutate(transect = factor(transect, levels = c("upland", "transition", "wte", "wc", "wetland")))
}
reorder_site = function(dat){
  dat %>% 
    mutate(site = factor(site, levels = c("CC", "PR", "OWC", "GCREW", "MSM", "GWI")))
}
reorder_horizon = function(dat){
  dat %>% 
    mutate(horizon = factor(horizon, levels = c("O", "A", "B")))
}
##run function below for principal components analysis
fit_pca_function_icr = function(icr_relabund_samples, sample_key){
  relabund_pca =
    icr_relabund_samples %>% 
    left_join(sample_key) %>% 
    drop_na() %>% 
    ungroup %>% 
    dplyr::select(-c(abund, total)) %>% 
    spread(Class, relabund) %>% 
    filter(!is.na(region)) %>% 
    replace(.,is.na(.),0)
  
  num = 
    relabund_pca %>% 
    dplyr::select(where(is.numeric))
  
  grp = 
    relabund_pca %>% 
    dplyr::select(where(is.character)) %>% 
    dplyr::mutate(row = row_number()) %>% 
    reorder_horizon()
  
  pca_int = prcomp(num, scale. = T)
  
  list(num = num,
       grp = grp,
       pca_int = pca_int)
}




##not using these functions
#fit_pca_function<- fit_pca_function_icr(icr_relabund, sample_key)
#compute_icr_pca<- compute_icr_pca(icr_relabund, sample_key)

##breaking up the inputs to fit_pca_function_icr to check the data products created within the function

relabund_pca =
  icr_relabund %>% 
  left_join(sample_key) %>% 
  drop_na() %>% 
  ungroup %>% 
  dplyr::select(-c(abund, total)) %>% 
  spread(Class, relabund) %>% 
  filter(!is.na(region)) %>% 
  replace(.,is.na(.),0)

num = 
  relabund_pca %>% 
  dplyr::select(where(is.numeric))

grp = 
  relabund_pca %>% 
  dplyr::select(where(is.character)) %>% 
  dplyr::mutate(row = row_number()) %>% 
  reorder_horizon()

pca_int = prcomp(num, scale. = T)

list(num = num,
     grp = grp,
     pca_int = pca_int)

###############################
##running the compute_icr_pca function step by step
###############################
sample_key = 
  sample_key %>% 
  mutate(transect = recode(transect, "wc"= "wetland", "wte"="wetland"))%>%
  mutate(region=recode(region, "CB"="Chesapeake", "WLE"="Erie"))%>%
  mutate(site=recode(site, "GCREW"="GCW", "CC"="CRC", "PR"="PTR"))
#sample_key$transect <- factor(sample_key$transect,levels = c("upland", "transition", "wetland"), ordered=TRUE)

pca_overall = fit_pca_function_icr(icr_relabund, sample_key)
pca_overall$grp$transect <- factor(pca_overall$grp$transect, levels=c("upland", "transition", "wetland"), ordered=TRUE)
pca_overall$grp$region <- factor(pca_overall$grp$region, levels=c("Chesapeake", "Erie"), ordered=TRUE)

pca_wle = fit_pca_function_icr(icr_relabund, sample_key %>% filter(region == "Erie"))
pca_wle$grp$transect <- factor(pca_wle$grp$transect, levels=c("upland", "transition", "wetland"), ordered=TRUE)
pca_wle$grp$site <- factor(pca_wle$grp$site, levels=c("CRC", "OWC", "PTR"), ordered=TRUE)

pca_cb = fit_pca_function_icr(icr_relabund, sample_key %>% filter(region == "Chesapeake"))
pca_cb$grp$transect <- factor(pca_cb$grp$transect, levels=c("upland", "transition", "wetland"), ordered=TRUE)
pca_cb$grp$site <- factor(pca_cb$grp$site, levels=c("GCW", "GWI", "MSM"), ordered=TRUE)

# PCA biplots
biplot_all = 
  ggbiplot(pca_overall$pca_int, obs.scale = 1, var.scale = 1,
           groups = pca_overall$grp$region,
           ellipse = FALSE, circle = FALSE, var.axes = TRUE, alpha = 0) +
  geom_point(size=3,stroke=1, alpha = 1,
             aes(shape = pca_overall$grp$transect,
                 color = groups))+
  scale_color_manual(values=c("#009E73", "#D55E00"))+
  scale_shape_manual(#breaks = (pca_overall$grp$transect, c("upland", "transition", "wetland")),
    #labels = c("upland", "transition", "wetland"),
    # values = c("upland"=15,"transition"=16,"wetland"=17))+
    values=c(15,16,17))+
  scale_shape(solid = FALSE)+
  xlim(-6,6)+
  ylim(-6, 6)+
  labs(color = "region", shape = "transect")+
  #labs(title = "FTICR: Chesapeake Bay and Lake Erie",
  #color = "", shape = "")+
  theme_classic()
#theme(legend.position = "top", legend.box = "vertical")

biplot_all

library(ggbiplot)

##NEW FINAL PLOT
biplot_all<-ggbiplot(pca_overall$pca_int, obs.scale = 1, var.scale = 1,
                     #groups = pca_overall$grp$region,
                     ellipse = FALSE, circle = FALSE, var.axes = TRUE, alpha = 0) +
  #title = "PCoA WLE CB") + 
  # scale_color_manual(values = c("#CC79A7", "#E69F00", "#0072B2", "black"))+
  scale_shape_manual(values=c(21,22, 24))+
  scale_fill_manual(values=c("#009E73", "#D55E00"))+
  geom_point(aes(shape=pca_overall$grp$transect, fill = pca_overall$grp$region), size = 6, color="black") +
  #geom_point(colour = "white", size = 1.5) + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme_classic()+
  xlim(-6,6)+
  ylim(-6, 6)+
  theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14))+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))+
  theme(legend.text = element_text(size=12),  legend.title = element_text(size=16))+
  labs(fill="region", shape="transect")+
  guides(fill=guide_legend(override.aes = list(shape=21), order=1),
         shape = guide_legend(override.aes = list(fill = "black"), order=2))

biplot_all


ggsave(file="biplots_cb_wle_colblind_07.05.24-ordered.tiff", plot = biplot_all, width = 8, height = 5)

ggsave(file="biplots_cb_wle_colblind_08.07.24-ordered-withBhor.tiff", plot = biplot_all, width = 8, height = 5)

##making it squared
ggsave(file="biplots_cb_wle_colblind_02.18.25-ordered-withBhor.tiff", plot = biplot_all, width = 7, height = 6)

library(viridis)

plasma_pal <- c("red", viridis::plasma(n = 6))

plasma_pal
#> [1] "red"       "#0D0887FF" "#6A00A8FF" "#B12A90FF" "#E16462FF" "#FCA636FF"
#> [7] "#F0F921FF"

biplot_wle = 
  ggbiplot(pca_wle$pca_int, obs.scale = 1, var.scale = 1,
           groups = pca_wle$grp$transect, 
           ellipse = FALSE, circle = FALSE, var.axes = TRUE, alpha = 0) +
  geom_point(size=3,stroke=1, alpha = 1,
             aes(shape = pca_wle$grp$site,
                 color = groups))+
  scale_color_manual(breaks = c("upland", "transition", "wetland"), 
                     values=c("#CC79A7", "#E69F00", "#0072B2")) +
  scale_shape_manual(values=c(15,16,17))+
  scale_shape(solid = FALSE)+
  labs(color = "transect", shape = "site")+
  xlim(-4,4)+
  ylim(-3.5,3.5)+
  theme_classic()
# theme(legend.position = "top", legend.box = "vertical")

biplot_wle

##new plot FINAL
biplot_wle<-ggbiplot(pca_wle$pca_int, obs.scale = 1, var.scale = 1,
                     #groups = pca_wle$grp$transect, 
                     ellipse = FALSE, circle = FALSE, var.axes = TRUE, alpha = 0) +
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


ggsave(file="biplot_wle_col_blind_07.05.24.tiff", plot = biplot_wle, width = 8, height = 5) ##this file got deleted by mistake, the older one without B horizon

ggsave(file="biplot_wle_col_blind_08.07.24_withBhor.tiff", plot = biplot_wle, width = 8, height = 5)

##making it squared
ggsave(file="biplot_wle_col_blind_02.18.25_withBhor.tiff", plot = biplot_wle, width = 7, height = 6)

biplot_cb = 
  ggbiplot(pca_cb$pca_int, obs.scale = 1, var.scale = 1,
           groups = pca_cb$grp$transect, 
           ellipse = FALSE, circle = FALSE, var.axes = TRUE, alpha = 0) +
  geom_point(size=3,stroke=1, alpha = 1,
             aes(shape = pca_cb$grp$site,
                 color = groups))+
  scale_color_manual(breaks = c("upland", "transition", "wetland"), 
                     values=c("#CC79A7", "#E69F00", "#0072B2"))+
  scale_shape_manual(values=c(15,16,17))+
  scale_shape(solid = FALSE)+
  labs(color = "transect", shape = "site")+
  xlim(-4,4)+
  ylim(-3.5,2)+
  theme_classic()
#theme(legend.position = "top", legend.box = "vertical")

biplot_cb

##new plot FINAL
biplot_cb<-ggbiplot(pca_cb$pca_int, obs.scale = 1, var.scale = 1,
                    #groups = pca_cb$grp$transect, 
                    ellipse = FALSE, circle = FALSE, var.axes = TRUE, alpha = 0) +
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


ggsave(file="biplot_cb_col_blind_07.05.24.tiff", plot = biplot_cb, width = 8, height = 5)
ggsave(file="biplot_cb_col_blind_08.07.24_withBhor.tiff", plot = biplot_cb, width = 8, height = 5)

##making it squared
ggsave(file="biplot_cb_col_blind_02.18.25_withBhor.tiff", plot = biplot_cb, width = 7, height = 6)


#library(patchwork)
#biplot_all = biplot_all+
 # biplot_cb + biplot_wle

#list(biplot_all = biplot_all,
#biplot_regions = biplot_regions)

##not using this combination approach
#ggsave(file="biplots_all_KFPcode_editedpanel_ChangingSiteAbrvtns.tiff", plot = biplot_all, width = 10, height = 7)

#####################################################
##Van krevlen by region and transect; total and unique
#####################################################

icr_data_trt<- read.csv("icr_long_treatments.csv")
icr_meta<- read.csv("icr_meta.csv")

##run function below

##gg_vankrev
gg_vankrev <- function(data,mapping){
  ggplot(data,mapping) +
    # plot points
    geom_point(size=0.5, alpha = 0.5) + # set size and transparency
    # axis labels
    ylab("H/C") +
    xlab("O/C") +
    # axis limits
    xlim(0,1.25) +
    ylim(0,2.5) +
    # add boundary lines for Van Krevelen regions
    geom_segment(x = 0.0, y = 1.5, xend = 1.2, yend = 1.5,color="black",linetype="longdash") +
    geom_segment(x = 0.0, y = 0.7, xend = 1.2, yend = 0.4,color="black",linetype="longdash") +
    geom_segment(x = 0.0, y = 1.06, xend = 1.2, yend = 0.51,color="black",linetype="longdash") +
    guides(colour = guide_legend(override.aes = list(alpha=1, size = 1)))
}

##not using chunk below
#plot_vk<- plot_vankrevelen(icr_data_trt, icr_meta)
#plot_vk
#plot_vk_unique<-plot_vankrevelen_unique(icr_data_trt, icr_meta)
#plot_vk_unique


##start here
library(patchwork)
library(ggpubr)


##ran the plots individually within the function to make the below code work

vk_domains = 
  icr_meta %>% 
  dplyr::select(formula, HC, OC, Class_detailed, Class) %>% 
  gg_vankrev(aes(x = OC, y = HC, color = Class))+
  scale_color_manual(breaks=c("aliphatic", "aromatic", "condensed aromatic", "unsaturated/lignin" ),
                     values= c("#B2DF8A", "#FEE391", "#810F7C", "#238443"))+
  ggtitle("FTICR: Domains in Chesapeake and Erie")+
  theme_classic()+
  theme(aspect.ratio = 1)##this makes the plots squared regardless of axes limits


data_hcoc = 
  icr_data_trt %>% 
  #filter(horizon != "B") %>% 
  left_join(icr_meta %>% dplyr::select(formula, HC, OC)) %>% 
  mutate(transect = recode(transect, "wc" = "wetland", "wte"="wetland")) %>% 
  mutate(site=recode(site,"GCREW"="GCW", "CC"="CRC", "PR"="PTR" )) %>%
  reorder_horizon() %>% reorder_transect()

vk_wle = 
  data_hcoc %>% 
  filter(region == "WLE") %>% 
  gg_vankrev(aes(x = OC, y = HC, color = transect))+
  stat_ellipse(level = 0.90, show.legend = F)+
  scale_color_manual(values = c("#CC79A7", "#E69F00", "#0072B2"))+
  labs(color = "transect",
       title = "FTICR: all peaks (blank corrected)",
       subtitle = "Erie")+
  facet_wrap(~site)+
  theme_classic()+
  theme(aspect.ratio = 1)+
  NULL

vk_cb = 
  data_hcoc %>% 
  filter(region == "CB" ) %>% 
  gg_vankrev(aes(x = OC, y = HC, color = transect))+
  stat_ellipse(level = 0.90, show.legend = F)+
  scale_color_manual(breaks = c("upland", "transition", "wetland"),
                     values = c("#CC79A7", "#E69F00", "#0072B2"))+
  labs(color = "transect",
       title = "FTICR: all peaks (blank corrected)",
       subtitle = "Chesapeake")+
  facet_wrap(~site)+
  theme_classic()+
  theme(aspect.ratio = 1)+
  NULL


list(vk_domains = vk_domains,
     vk_wle = vk_wle,
     vk_cb = vk_cb)



unique_hcoc = 
  icr_data_trt %>% 
  # filter(horizon != "B") %>% 
  #filter(!(site == "MSM" & horizon == "A")) %>% 
  group_by(formula, region, site, horizon) %>% 
  dplyr::mutate(n = n()) %>% 
  filter(n == 1) %>% 
  left_join(icr_meta %>% dplyr::select(formula, HC, OC)) %>% 
  mutate(transect = recode(transect, "wc" = "wetland", "wte"="wetland")) %>% 
  mutate(site=recode(site,"GCREW"="GCW", "CC"="CRC", "PR"="PTR" )) %>%
  reorder_horizon() %>% reorder_transect()


vk_unique_wle = 
  unique_hcoc %>% 
  filter(region == "WLE") %>% 
  gg_vankrev(aes(x = OC, y = HC, color = transect))+
  stat_ellipse(level = 0.90, show.legend = F)+
  scale_color_manual(breaks = c("upland", "transition", "wetland"),
                     values = c("#CC79A7", "#E69F00", "#0072B2"))+
  labs(color = "transect")+
  facet_grid(~site)+
  labs(title = "FTICR Unique Peaks",
       subtitle = "Erie")+
  theme_classic()+
  theme(aspect.ratio = 1)


vk_unique_cb = 
  unique_hcoc %>% 
  filter(region == "CB") %>% 
  gg_vankrev(aes(x = OC, y = HC, color = transect))+
  stat_ellipse(level = 0.90, show.legend = F)+
  scale_color_manual(breaks = c("upland", "transition", "wetland"),
                     values = c("#CC79A7", "#E69F00", "#0072B2"))+
  labs(color = "transect")+
  facet_grid(~site)+
  labs(title = "FTICR Unique Peaks",
       subtitle = "Chesapeake")+
  theme_classic()+
  theme(aspect.ratio = 1)


plot_vk= ggarrange(vk_domains, vk_cb, vk_wle, nrow=3)
plot_vk
plot_vk_unique= ggarrange(vk_unique_cb, vk_unique_wle, nrow=2)
plot_vk_unique
ggsave(file="VK_unique_cb_wle_07.05.24.tiff", plot = plot_vk_unique, width = 13, height = 7)
ggsave(file="VK_domains_cb_wle_all_07.05.24.tiff", plot = plot_vk, width = 20, height = 7)

##making it squared
ggsave(file="VK_unique_cb_wle_02.18.25.tiff", plot = plot_vk_unique, width = 10, height = 8)
ggsave(file="VK_domains_cb_wle_all_02.18.25.tiff", plot = plot_vk, width = 8, height = 10)



##permanova 


# permanova -----------------------------------------------------------

##function combining region
compute_permanova = function(icr_relabund_samples){
  relabund_wide = 
    icr_relabund_samples %>% 
    left_join(sample_key) %>% 
    # filter(horizon != "B") %>% 
    ungroup %>% 
    dplyr::select(-c(abund, total)) %>% 
    spread(Class, relabund) %>% 
    filter(!is.na(region)) %>% 
    replace(.,is.na(.),0)
  
  permanova_fticr_all = 
    adonis(relabund_wide %>% dplyr::select(where(is.numeric)) ~ (region + transect + site + horizon)^2, 
           data = relabund_wide)
  broom::tidy(permanova_fticr_all$aov.tab)
}

#function within region for CB
compute_permanova_within_cb = function(icr_relabund_samples){
  relabund_wide = 
    icr_relabund_samples %>% 
    left_join(sample_key) %>% 
    # filter(horizon != "B") %>% 
    ungroup %>% 
    dplyr::select(-c(abund, total)) %>% 
    spread(Class, relabund) %>% 
    filter(!is.na(region)) %>% 
    replace(.,is.na(.),0)
  
  permanova_fticr_all = 
    adonis(relabund_wide %>% dplyr::select(where(is.numeric)) ~ (transect + site + horizon)^2, 
           data = relabund_wide)
  broom::tidy(permanova_fticr_all$aov.tab)
}

#function within region for WLE
compute_permanova_within_wle = function(icr_relabund_samples){
  relabund_wide = 
    icr_relabund_samples %>% 
    left_join(sample_key) %>% 
    # filter(horizon != "B") %>% 
    ungroup %>% 
    dplyr::select(-c(abund, total)) %>% 
    spread(Class, relabund) %>% 
    filter(!is.na(region)) %>% 
    replace(.,is.na(.),0)
  
  permanova_fticr_all = 
    adonis(relabund_wide %>% dplyr::select(where(is.numeric)) ~ (transect + site)^2, 
           data = relabund_wide)
  broom::tidy(permanova_fticr_all$aov.tab)
}

icr_relabund_region<- icr_relabund %>% left_join(sample_key, by="sample_label")
icr_relabund_cb<- icr_relabund_region %>% filter(region=="Chesapeake Bay")
icr_relabund_wle<- icr_relabund_region %>% filter(region=="Lake Erie")

##combining regions
permanova<- compute_permanova(icr_relabund)
permanova

term                df SumsOfSqs  MeanSqs F.Model     R2 p.value
<chr>            <dbl>     <dbl>    <dbl>   <dbl>  <dbl>   <dbl>
  1 region               1    0.0896  0.0896    41.3  0.0601   0.001
2 transect             2    0.0607  0.0304    14.0  0.0408   0.001
3 site                 4    0.410   0.102     47.3  0.275    0.001
4 horizon              2    0.104   0.0521    24.0  0.0700   0.001
5 region:transect      2    0.162   0.0812    37.4  0.109    0.001
6 transect:site        8    0.153   0.0191     8.80 0.102    0.001
7 transect:horizon     2    0.0974  0.0487    22.5  0.0654   0.001
8 site:horizon         2    0.0255  0.0128     5.89 0.0171   0.015
9 Residuals          179    0.388   0.00217   NA    0.260   NA    
10 Total              202    1.49   NA         NA    1       NA  


##CB
permanova<- compute_permanova_within_cb(icr_relabund_cb)
permanova

# A tibble: 8 × 7
term                df SumsOfSqs  MeanSqs F.Model     R2 p.value
<chr>            <dbl>     <dbl>    <dbl>   <dbl>  <dbl>   <dbl>
  1 transect             2    0.118   0.0591    19.7  0.102    0.001
2 site                 2    0.387   0.193     64.5  0.334    0.001
3 horizon              2    0.145   0.0725    24.2  0.125    0.001
4 transect:site        4    0.0815  0.0204     6.79 0.0704   0.001
5 transect:horizon     2    0.0974  0.0487    16.2  0.0841   0.001
6 site:horizon         2    0.0255  0.0128     4.25 0.0220   0.012
7 Residuals          101    0.303   0.00300   NA    0.262   NA    
8 Total              115    1.16   NA         NA    1       NA    




#WLE
permanova<- compute_permanova_within_wle(icr_relabund_wle)
permanova

# A tibble: 5 × 7
term             df SumsOfSqs  MeanSqs F.Model    R2 p.value
<chr>         <dbl>     <dbl>    <dbl>   <dbl> <dbl>   <dbl>
  1 transect          2    0.0603  0.0302     27.7 0.248   0.001
2 site              2    0.0267  0.0133     12.3 0.110   0.001
3 transect:site     4    0.0710  0.0178     16.3 0.292   0.001
4 Residuals        78    0.0849  0.00109    NA   0.349  NA    
5 Total            86    0.243  NA          NA   1      NA    


