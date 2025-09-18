
## **********************************************************
## Load libraries and source files
## **********************************************************

rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd("microbeBiogeography/analysis/NetworkTurnover/")

source("src/chao.R")
source("src/betaNet.R")
source("src/writeResultsTable.R")
source("src/networkTurnover.R")
source("src/distDecay.R")
source("src/standardize_weights_microbes.R")

library(ggplot2)
library(lme4)
library(lmerTest)
library(igraph)
library(ggpubr)
library(emmeans)
library(bipartite)
library(dplyr)
library(fields)
library(brms)
library(tidybayes)
library(gridExtra)
library(bayesplot)
library(glmmTMB)
library(performance)
library(betapart)
library(grid)
library(gridExtra)
library(ggborderline)

#load("../../microbeBiogeographyData.Rdata") ## TODO update with correct filepath before ms submission
load("../../../skyIslands/data/spec_RBCL_16s.Rdata")

microbe_data <- prepMicrobeWeights(spec.net)

these_sites <- microbe_data %>%
  select(WeightsMicrobe, SiteYearSr) %>%
  filter(WeightsMicrobe == 1) %>%
  select(!WeightsMicrobe) %>%
  distinct()

site_names <- these_sites$SiteYearSr
#"C:/Users/rah10/Dropbox (University of Oregon)/microbeBiogeography/analysis/NetworkTurnover"

## TODO: calculate bray curtis dissim for bee community at each site
veg_df <- read.csv("../../../SkyIslands/data/veg_bloom_quad_sp.csv") %>%
  mutate(SiteYearSR = paste(Site, Year, SampleRound, sep=";")) %>%
  group_by(SiteYearSR, PlantGenusSpecies) %>%
  mutate(TotalAbund = sum(FloralAbundance)) %>%
  select(SiteYearSR, PlantGenusSpecies, TotalAbund) %>%
  filter(SiteYearSR %in% site_names) %>%
  distinct() %>%
  pivot_wider(
    names_from = PlantGenusSpecies,
    values_from = TotalAbund,
    values_fill = 0
  ) %>%
arrange(SiteYearSR) %>%
  ungroup() %>%
  select(!SiteYearSR)



# Flower dissimilarity matrix
veg_dist <- vegan::vegdist(veg_df, method = "bray")

geo <- microbe_data %>%
  filter(WeightsMicrobe == 1) %>%
  group_by(SiteYearSr) %>%
  select(Long, Lat, SiteYearSr) %>%
  distinct() %>%
  arrange(SiteYearSr) 
geo$SiteYearSr <- NULL   

d.geo <- distm(geo, fun = distHaversine)
dist.geo <- as.dist(d.geo)/1000

set.seed(1298)
dist_decay_model <- betapart::decay.model(veg_dist,
                                              dist.geo,
                                              y.type="dissimilarities",
                                              model.type = "exponential",
                                              perm=999)

floral.decay <- plot_decay_ggplot_single(dist_decay_model,
                                      xlab="Geographic Distance (km)",
                                      ylab="Floral Community Dissimilarity",
                                      mod1color='purple',
                                      col = "black",
                                      lty = "solid",
                                      lwd = 1.5,
                                      cex = 1)

floral.decay

save(floral.decay, "figures/floral_distance_decay.pdf")
#   
# get_bee_community_dissim <- function(data, host.type, sim.type){
# 
#   if(host.type=="All"){
#     abund <- data %>%
#       filter(WeightsMicrobe == 1) %>%
#       select(Site, GenusSpecies) %>%
#       group_by(Site, GenusSpecies) %>%
#       arrange(Site) %>%
#       filter(GenusSpecies != "") %>%
#       summarize(HostAbund = n()) %>%
#       pivot_wider(names_from = GenusSpecies, values_from = HostAbund) %>%
#       mutate_all(~replace(., is.na(.), 0))
#     sitenames <- abund$Site
#     abund$Site <- NULL
#     rownames(abund) <- sitenames
#     
#     #distance matrix of sites
#     geo <- data %>%
#       filter(WeightsMicrobe == 1) %>%
#       group_by(Site) %>%
#       select(Site, Long, Lat) %>%
#       distinct() %>%
#       arrange(Site)
#     rownames(geo) <- sitenames
#     geo$Site <- NULL
#   } 
#   if(host.type=="Social"){
#     abund <- data %>%
#       filter(WeightsMicrobe == 1) %>%
#       filter(Genus %in% c("Apis", "Bombus")) %>%
#       select(Site, GenusSpecies) %>%
#       group_by(Site, GenusSpecies) %>%
#       arrange(Site) %>%
#       filter(GenusSpecies != "") %>%
#       summarize(HostAbund = n()) %>%
#       pivot_wider(names_from = GenusSpecies, values_from = HostAbund) %>%
#       mutate_all(~replace(., is.na(.), 0))
#     sitenames <- abund$Site
#     abund$Site <- NULL
#     rownames(abund) <- sitenames
#     
#     #distance matrix of sites
#     geo <- data %>%
#       filter(WeightsMicrobe == 1) %>%
#       filter(Genus %in% c("Apis", "Bombus")) %>%
#       group_by(Site) %>%
#       select(Site, Long, Lat) %>%
#       distinct() %>%
#       arrange(Site)
#     rownames(geo) <- sitenames
#     geo$Site <- NULL
#   }
#   if(host.type=="Solitary"){
#     abund <- data %>%
#       filter(WeightsMicrobe == 1) %>%
#       filter(Genus %in% c("Anthophora", "Andrena", "Melissodes", "Megachile")) %>%
#       select(Site, GenusSpecies) %>%
#       group_by(Site, GenusSpecies) %>%
#       arrange(Site) %>%
#       filter(GenusSpecies != "") %>%
#       summarize(HostAbund = n()) %>%
#       pivot_wider(names_from = GenusSpecies, values_from = HostAbund) %>%
#       mutate_all(~replace(., is.na(.), 0))
#     sitenames <- abund$Site
#     abund$Site <- NULL
#     rownames(abund) <- sitenames
#     
#     #distance matrix of sites
#     geo <- data %>%
#       filter(WeightsMicrobe == 1) %>%
#       filter(Genus %in% c("Anthophora", "Andrena", "Melissodes", "Megachile")) %>%
#       group_by(Site) %>%
#       select(Site, Long, Lat) %>%
#       distinct() %>%
#       arrange(Site)
#     rownames(geo) <- sitenames
#     geo$Site <- NULL
#   }
#   
#   
#   
#   dist.abund <- vegdist(abund, method = "bray")
#   
#   if(sim.type == "similarities"){
#     dist.sim <- 1 - dist.abund
#   }
#   
#   d.geo <- distm(geo, fun = distHaversine)
#   dist.geo <- as.dist(d.geo)/1000
#   
#   #abundance vs geographic mantel test
#   set.seed(1298)
#   abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 999, na.rm = TRUE)
#   print(abund_geo)
#   
#   if(sim.type == "similarities"){
#     set.seed(1298)
#     dist_decay_model <- betapart::decay.model(dist.sim,
#                                               dist.geo,
#                                               y.type="similarities",
#                                               model.type = "exponential",
#                                               perm=999)
# 
#   }
#   if(sim.type == "dissimilarities"){
#     set.seed(1298)
#     dist_decay_model <- betapart::decay.model(dist.sim,
#                                               dist.geo,
#                                               y.type="dissimilarities",
#                                               model.type = "exponential",
#                                               perm=999)
# 
#   }
#   
#   return(dist_decay_model)
#   
# }

## prep microbe weights
spec.net <- prepMicrobeWeights(spec.net)


##
social_model <- get_bee_community_dissim(spec.net, host.type="Social", sim.type="similarities")
solitary_model <- get_bee_community_dissim(spec.net, host.type="Solitary", sim.type="similarities")

plot_decay_ggplot_combined(social_model,
                           solitary_model,
                           mod1color='blue',
                           group1label = 'Social',
                           mod2color='gold',
                           group2label = 'Solitary',
                           alpha1=0.5,
                           alpha2=0.5,
                           lty1='solid',
                           lty2='solid',
                           xlab="Geographic Distance (km)",
                           ylab='Bee Community Similarity', 
                           add.points=TRUE, 
                           add.legend=TRUE,
                           legend.title='Host Type',
                           log.dist = FALSE)


bombus_model <- get_bee_community_dissim(spec.net, host.type="Bombus", sim.type="similarities")
apis_model <- get_bee_community_dissim(spec.net, host.type="Apis", sim.type="similarities")


plot_decay_ggplot_combined(bombus_model,
                           apis_model,
                           mod1color='purple',
                           group1label = 'Bombus',
                           mod2color='coral',
                           group2label = 'Apis',
                           alpha1=0.5,
                           alpha2=0.5,
                           lty1='solid',
                           lty2='solid',
                           xlab="Geographic Distance (km)",
                           ylab='Bee Community Similarity', 
                           add.points=TRUE, 
                           add.legend=TRUE,
                           legend.title='Host Genus',
                           log.dist = FALSE)
##
# what we need
# site to site comparison of host turnover
# columns should be Species
# rows should be sites
# values should be abundance
abund <- data %>%
  filter(WeightsMicrobe == 1) %>%
  select(Site, GenusSpecies) %>%
  group_by(Site, GenusSpecies) %>%
  arrange(Site) %>%
  filter(GenusSpecies != "") %>%
  summarize(HostAbund = n()) %>%
  pivot_wider(names_from = GenusSpecies, values_from = HostAbund) %>%
  mutate_all(~replace(., is.na(.), 0))
sitenames <- abund$Site
abund$Site <- NULL
rownames(abund) <- sitenames

#distance matrix of sites
geo <- data %>%
  filter(WeightsMicrobe == 1) %>%
  group_by(Site) %>%
  select(Site, Long, Lat) %>%
  distinct() %>%
  arrange(Site)
rownames(geo) <- sitenames
geo$Site <- NULL



#abundance data frame - bray curtis dissimilarity
dist.abund <- vegdist(abund, method = "bray")
dist.sim <- 1 - dist.abund

#geographic data frame - haversine distance in m (takes a df with lat and long and calculates dist)
d.geo <- distm(geo, fun = distHaversine)
dist.geo <- as.dist(d.geo)/1000

#abundance vs geographic mantel test
abund_geo  = mantel(dist.sim, dist.geo, method = "spearman", permutations = 999, na.rm = TRUE)
print(abund_geo)

dist_decay_model <- betapart::decay.model(dist.sim,
                                          dist.geo,
                                          y.type="similarities",
                                          model.type = "exponential",
                                          perm=999)
dist_decay_model

bee.decay <- plot_decay_ggplot_single(dist_decay_model,
                                      xlab="Geographic Distance (km)",
                                      ylab="Bee Community Similarity",
                                      mod1color='navy',
                                      col = "black",
                                      lty = "solid",
                                      lwd = 1.5,
                                      cex = 1)
bee.decay
