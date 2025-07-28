
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
load("../../../skyIslands/data/networks/microNets.RData")

## set hosts="All" to run models for full host dataset with the full list of solitary and social strong HAMs
## set hosts="Social" to run mods for social host dataset with social strong HAMS
## set hosts="Solitary" to run mods for solitary host dataset with solitary strong HAMS
## set hosts="AllPathogens" to run models for full host dataset with the pathogenic microbes

hosts="Social"

## **********************************************************
## Prep obligate and transient networks
## **********************************************************

full_obligate_list <- c("Lactobacillaceae",
                        "Bifidobacteriaceae",
                        "Neisseriaceae",
                        "Orbaceae",
                        "Bartonella",
                        "Acetobacteraceae",
                        "Bacillaceae",
                        "Burkholderiaceae",
                        "Clostridiaceae",
                        "Comamonadaceae",
                        "Enterobacteriaceae",
                        "Lachnospiraceae",
                        "Methylobacteriaceae",
                        "Moraxellaceae",
                        "Sphingomonadaceae", 
                        "Oxalobacteraceae"
)

social_obligate_list <- c("Lactobacillaceae",
                          "Bifidobacteriaceae",
                          "Neisseriaceae",
                          "Orbaceae",
                          "Bartonella",
                          "Acetobacteraceae",
                          "Hafnia",
                          "Wolbachia",
                          "Erwinia"
)

solitary_obligate_list <- c("Acetobacteraceae",
                            "Bacillaceae",
                            "Burkholderiaceae",
                            "Clostridiaceae",
                            "Comamonadaceae",
                            "Enterobacteriaceae",
                            "Lachnospiraceae",
                            "Lactobacillaceae",
                            "Methylobacteriaceae",
                            "Moraxellaceae",
                            "Sphingomonadaceae", 
                            "Oxalobacteraceae",
                            "Hafnia",
                            "Wolbachia",
                            "Erwinia"
)

pathogens <- c("Wolbachia",
               "Erwinia",
               "Hafnia")



## Social host community, social obligates
if(hosts=="Social"){
  social_obligate_network <- prep_obligate_network(raw_network=spNet_micro,
                                                   these_obligates = social_obligate_list,
                                                   genera_to_keep=c("Bombus", "Apis"))
  
  social_transient_network <- prep_transient_network(raw_network=spNet_micro,
                                                     these_obligates = social_obligate_list,
                                                     genera_to_keep=c("Bombus", "Apis"))
}


## **********************************************************
## Run network betalinkr function and prep output table
## **********************************************************

#source("src/betalinkrPrep.R")


if(hosts=='Social'){
  
  ## SOCIAL HOSTS SOCIAL OBLIGATE NETWORKS
  find_sites_for_betalinkr(social_obligate_network)
  
  ## enter the site matrices printed above 
  obligate_social_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, SM, SC, RP),
                                              partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)
  
  obligate_social_betalink_clean <- fix_betalinkr_output(obligate_social_betalink)
  
  ## SOCIAL HOSTS SOCIAL transient NETWORKS
  find_sites_for_betalinkr(social_transient_network)
  
  ## enter the site matrices printed above 
  transient_social_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, SM, SC, RP),
                                               partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)
  
  transient_social_betalink_clean <- fix_betalinkr_output(transient_social_betalink)
  
  
}

hosts="Solitary"
## Solitary host community, solitary obligates
if(hosts=="Solitary"){
  solitary_obligate_network <- prep_obligate_network(raw_network=spNet_micro,
                                                     these_obligates = solitary_obligate_list,
                                                     genera_to_keep=c("Melissodes", "Megachile", "Anthophora", "Andrena"))
  
  solitary_transient_network <- prep_transient_network(raw_network=spNet_micro,
                                                       these_obligates = solitary_obligate_list,
                                                       genera_to_keep=c("Melissodes", "Megachile", "Anthophora", "Andrena"))
}

if(hosts=='Solitary'){
  
  ## SOLITARY HOSTS SOLITARY OBLIGATE NETWORKS
  find_sites_for_betalinkr(solitary_obligate_network)
  
  ## enter the site matrices printed above 
  obligate_solitary_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, SM, SC, RP),
                                                partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)
  
  obligate_solitary_betalink_clean <- fix_betalinkr_output(obligate_solitary_betalink)
  
  ## SOCIAL HOSTS SOCIAL transient NETWORKS
  find_sites_for_betalinkr(solitary_transient_network)
  
  ## enter the site matrices printed above 
  transient_solitary_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, SM, SC, RP),
                                                 partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)
  
  transient_solitary_betalink_clean <- fix_betalinkr_output(transient_solitary_betalink)
}



## **********************************************************
## Pairwise bray curtis dissimilarity calculation and plots
## **********************************************************

run.decay.mictype.mods=FALSE

## prep microbe weights
spec.net <- prepMicrobeWeights(spec.net)

hosts="Social"
 if (hosts=="Social"){
   if (run.decay.mictype.mods == TRUE){
     #load("../../../skyIslands/data/spec_RBCL_16s.Rdata")
     meta_cols <- c('UniqueID', 'Genus', 'Species', 'GenusSpecies', 'Site', 'Lat', 'Long',"WeightsObligateSocialMicrobe", "WeightsTransientSocialMicrobe")

     spec16s <- spec.net %>%
       filter(Apidae == 1) %>%
       select(all_of(meta_cols), starts_with('16s')) %>%
       na.omit()

     ob_model <- microbe_type_decay_model(spec16s, 'ObligateSocial', model.type = 'power', decay.type = "sim") ## power if plan to log transform x axis
     trans_model <- microbe_type_decay_model(spec16s, 'TransientSocial', model.type='power', decay.type = "sim")
     ## save out models
     save(ob_model,
          trans_model,
          file="../../../skyIslands/analysis/microbiome/saved/decay_mictype_mods_social_power.Rdata") ## TODO update filepaths
   } else {
     load("../../../skyIslands/analysis/microbiome/saved/decay_mictype_mods_social_power.Rdata") ## TODO update filepaths
   }
 }

#load("C:/Users/rah10/University of Oregon Dropbox/Rebecca Hayes/skyIslands/analysis/microbiome/saved/decay_mictype_mods_social.Rdata")

## microbe type comparison
social_bray <- plot_decay_ggplot_combined(ob_model,
                                          trans_model,
                                          mod1color='darkgreen',
                                          mod2color='darkorange',
                                          alpha1=0.003,
                                          alpha2=0.005,
                                          lty1='solid',
                                          lty2='solid',
                                          xlab="log Geographic Distance (km)",
                                          ylab='Social Associate Similarity',
                                          add.points=TRUE,
                                          log.dist=TRUE)

#plot_decay_ggplot_combined(ob_model, trans_model)

social_bray <- social_bray + labs(tag="A")
social_bray

hosts="Solitary"
run.decay.mictype.mods=TRUE
dist.decay.type="power"
if (hosts=="Solitary"){
  if (run.decay.mictype.mods == TRUE){
    #load("../../../skyIslands/data/spec_RBCL_16s.Rdata")
    meta_cols <- c('UniqueID', 'Genus', 'Species', 'GenusSpecies', 'Site', 'Lat', 'Long',"WeightsObligateSolitaryMicrobe", "WeightsTransientSolitaryMicrobe")
    
    spec16s <- spec.net %>%
      filter(Apidae == 1) %>%
      select(all_of(meta_cols), starts_with('16s')) %>%
      na.omit()
    
    ob_model <- microbe_type_decay_model(spec16s, 'ObligateSolitary', model.type = dist.decay.type, decay.type = "sim")
    trans_model <- microbe_type_decay_model(spec16s, 'TransientSolitary', model.type=dist.decay.type, decay.type = "sim")
    #browser()
    ## save out models
    save(ob_model,
         trans_model,
         file="../../../skyIslands/analysis/microbiome/saved/decay_mictype_mods_solitary_power.Rdata") ## TODO update filepaths
  } else {
    load("../../../skyIslands/analysis/microbiome/saved/decay_mictype_mods_solitary_power.Rdata") ## TODO update filepaths
  }
}
solitary_bray <- plot_decay_ggplot_combined(ob_model,
                                          trans_model,
                                          mod1color='darkgreen',
                                          mod2color='darkorange',
                                          alpha1=0.003,
                                          alpha2=0.005,
                                          lty1='solid',
                                          lty2='solid',
                                          xlab="log Geographic Distance (km)",
                                          ylab='Solitary Associate Similarity', add.points=TRUE, log.dist = TRUE)

solitary_bray <- solitary_bray + labs(tag="B")
solitary_bray
## make panels for interaction turnover
load("C:/Users/rah10/University of Oregon Dropbox/Rebecca Hayes/skyIslands/analysis/microbiome/saved/turnover_mods_social.Rdata")

## B. Interaction turnover
social.int.plot <- plot_network_turnover_mod_compare(mod1=int.obligate.mod,
                                                     mod2=int.transient.mod,
                                                     this.network1=obligate_social_betalink_clean,
                                                     this.network2=transient_social_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="WholeNetworkLinks",
                                                     label="Total Network Dissimilarity")

social.int.plot[[1]]

social_panelC <- social.int.plot[[1]] + labs(tag="C")



load("C:/Users/rah10/University of Oregon Dropbox/Rebecca Hayes/skyIslands/analysis/microbiome/saved/turnover_mods_solitary.Rdata")

## B. Interaction turnover
solitary.int.plot <- plot_network_turnover_mod_compare(mod1=int.obligate.mod,
                                                     mod2=int.transient.mod,
                                                     this.network1=obligate_solitary_betalink_clean,
                                                     this.network2=transient_solitary_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="WholeNetworkLinks",
                                                     label="Total Network Dissimilarity")
solitary.int.plot[[1]]

solitary_panelD <- solitary.int.plot[[1]] + labs(tag="D")


bray_plots=FALSE
if(bray_plots==TRUE){
  # Arrange all panels in the PDF output
  pdf("figures/bray_combined.pdf", width = 7, height = 7)  
  grid.arrange(
    social_bray,
    solitary_bray,
    social_panelC,
    solitary_panelD,
    ncol = 2
  )
  dev.off()
}


### temp
pdf("figures/bray_similarity_logdist.pdf", width = 7, height = 3.5)  
grid.arrange(
  social_bray,
  solitary_bray,
  #social_panelC,
  #solitary_panelD,
  ncol = 2
)
dev.off()
