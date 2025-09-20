
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

## prep microbe weights
spec.net <- prepMicrobeWeights(spec.net)


## model social and solitary dissimilarity
social_model <- get_bee_community_dissim(spec.net, host.type="Social", sim.type="dissimilarities")
solitary_model <- get_bee_community_dissim(spec.net, host.type="Solitary", sim.type="dissimilarities")

## plot distance decay for social versus solitary
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
                           ylab='Bee Community Dissimilarity', 
                           add.points=TRUE, 
                           add.legend=TRUE,
                           legend.title='Host Type',
                           log.dist = FALSE)


## and again for just social, split by genus
bombus_model <- get_bee_community_dissim(spec.net, host.type="Bombus", sim.type="dissimilarities")
apis_model <- get_bee_community_dissim(spec.net, host.type="Apis", sim.type="dissimilarities")


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
                           ylab='Bee Community Dissimilarity', 
                           add.points=TRUE, 
                           add.legend=TRUE,
                           legend.title='Host Genus',
                           log.dist = FALSE)

