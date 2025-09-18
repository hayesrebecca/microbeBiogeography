
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

load("../../microbeBiogeographyData.Rdata")

microbe_data <- prepMicrobeWeights(spec.net)

these_sites <- microbe_data %>%
  select(WeightsMicrobe, SiteYearSr) %>%
  filter(WeightsMicrobe == 1) %>%
  select(!WeightsMicrobe) %>%
  distinct()

site_names <- these_sites$SiteYearSr

## calculate bray curtis dissim for veg community at each site
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
