indicator_cb <- read.csv("indicator_transect_CB.csv")
indicator_wle<- read.csv("indicator_transect_WLE.csv")
network_features<- read.csv("network_features_otu.csv")

library(tidyverse)
 
upland_cb_ind<- indicator_cb %>% filter(s.upland==1)
transition_cb_ind<- indicator_cb %>% filter(s.transition==1)
wetland_cb_ind<- indicator_cb %>% filter(s.wetland==1)

upland_wle_ind<- indicator_wle %>% filter(s.upland==1)
transition_wle_ind<- indicator_wle %>% filter(s.transition==1)
wetland_wle_ind<- indicator_wle %>% filter(s.wetland==1)


network_cb_upland<- network_features %>% filter (region=="Chesapeake Bay" & transect=="upland")
network_cb_wetland<- network_features %>% filter (region=="Chesapeake Bay" & transect=="wetland")
network_wle_upland<- network_features %>% filter (region=="Lake Erie" & transect=="upland")
network_wle_transition<- network_features %>% filter (region=="Lake Erie" & transect=="transition")
network_wle_wetland<- network_features %>% filter (region=="Lake Erie" & transect=="wetland")

##matching these indicators in network data
##all present as indicators
match_cb_upland<- network_cb_upland %>% inner_join(upland_cb_ind, by = "OTU")
write.csv(match_cb_upland, "match_ind_network_cb_upland.csv")

##saving the unmatched ones
nomatch_cb_upland<- network_cb_upland %>% anti_join(match_cb_upland, by="OTU")
write.csv(nomatch_cb_upland, "nomatch_with_indicator_cb_upland.csv")

#all present as indicators
match_cb_wetland<- network_cb_wetland %>% inner_join(wetland_cb_ind, by = "OTU")
write.csv(match_cb_wetland, "match_ind_network_cb_wetland.csv")

##saving the unmatched ones
nomatch_cb_wetland<- network_cb_wetland %>% anti_join(match_cb_wetland, by="OTU")
write.csv(nomatch_cb_wetland, "nomatch_with_indicator_cb_wetland.csv")

##all match except one. one of the rhizobiales in module 33 is not an indicator taxa, but there is another rhizobiales in module 56 which is a wle upland indicator
match_wle_upland<- network_wle_upland %>% inner_join(upland_wle_ind, by = "OTU")
write.csv(match_wle_upland, "match_ind_network_wle_upland.csv")

##saving the unmatched ones
nomatch_wle_upland<- network_wle_upland %>% anti_join(match_wle_upland, by="OTU")
write.csv(nomatch_wle_upland, "nomatch_with_indicator_wle_upland.csv")

## only one match. one rhizobiales OTU in module 26
match_wle_transition<- network_wle_transition %>% inner_join(transition_wle_ind, by = "OTU")
write.csv(match_wle_transition, "match_ind_network_wle_transition.csv")

##saving the unmatched ones
nomatch_wle_transition<- network_wle_transition %>% anti_join(match_wle_transition, by="OTU")
write.csv(nomatch_wle_transition, "nomatch_with_indicator_wle_transition.csv")

##there are 34 matches out of 46 networks OTUs
match_wle_wetland<- network_wle_wetland %>% inner_join(wetland_wle_ind, by = "OTU")
write.csv(match_wle_wetland, "match_ind_network_wle_wetland.csv")


##saving the unmatched ones
nomatch_wle_wetland<- network_wle_wetland %>% anti_join(match_wle_wetland, by="OTU")
write.csv(nomatch_wle_wetland, "nomatch_with_indicator_wle_wetland.csv")

