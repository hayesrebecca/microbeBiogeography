
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

hosts="AllPathogens"

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
                          "Acetobacteraceae"
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
                            "Oxalobacteraceae"
                            )

pathogens <- c("Wolbachia",
               "Erwinia",
               "Hafnia")


## Full host community, full set of strong host associates
if(hosts=="All"){
  only_obligate_network <- prep_obligate_network(raw_network=spNet_micro, 
                                                 these_obligates=full_obligate_list
                                                 )
  
  only_transient_network <- prep_transient_network(raw_network=spNet_micro,
                                                   these_obligates=full_obligate_list
                                                   )
}
## Social host community, social obligates
if(hosts=="Social"){
  social_obligate_network <- prep_obligate_network(raw_network=spNet_micro,
                                                   these_obligates = social_obligate_list,
                                                   genera_to_keep=c("Bombus", "Apis"))
  
  social_transient_network <- prep_transient_network(raw_network=spNet_micro,
                                                     these_obligates = social_obligate_list,
                                                     genera_to_keep=c("Bombus", "Apis"))
}
## Solitary host community, solitary obligates
if(hosts=="Solitary"){
  solitary_obligate_network <- prep_obligate_network(raw_network=spNet_micro,
                                                   these_obligates = solitary_obligate_list,
                                                   genera_to_keep=c("Melissodes", "Megachile", "Anthophora", "Andrena"))
  
  solitary_transient_network <- prep_transient_network(raw_network=spNet_micro,
                                                     these_obligates = solitary_obligate_list,
                                                     genera_to_keep=c("Melissodes", "Megachile", "Anthophora", "Andrena"))
}
## Full host community, all pathogens
if(hosts=="AllPathogens"){
  pathogen_obligate_network <- prep_obligate_network(raw_network=spNet_micro, 
                                                 these_obligates=pathogens
  )
  
  pathogen_transient_network <- prep_transient_network(raw_network=spNet_micro,
                                                   these_obligates=c(pathogens, full_obligate_list)
  )
}

## **********************************************************
## Run network betalinkr function and prep output table
## **********************************************************

#source("src/betalinkrPrep.R")

if(hosts=='All'){
  ## ALL HOSTS ALL STRONG ASSOCIATE NETWORKS
  find_sites_for_betalinkr(only_obligate_network)
  
  ## enter the site matrices printed above 
  obligate_poll_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, SM, SC, RP),
                                            partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)
  
  obligate_poll_betalink_clean <- fix_betalinkr_output(obligate_poll_betalink)

  ## ALL HOSTS ALL WEAK ASSOCIATE NETWORKS
  find_sites_for_betalinkr(only_transient_network)
  
  ## enter the site matrices printed above 
  transient_poll_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, SM, SC, RP),
                                            partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)
  
  transient_poll_betalink_clean <- fix_betalinkr_output(transient_poll_betalink)
}

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

if(hosts=='AllPathogens'){
  ## ALL HOSTS ALL STRONG ASSOCIATE NETWORKS
  find_sites_for_betalinkr(pathogen_obligate_network)
  
  ## enter the site matrices printed above 
  obligate_pathogen_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, SM, SC, RP),
                                            partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)
  
  obligate_pathogen_betalink_clean <- fix_betalinkr_output(obligate_pathogen_betalink)
  
  ## ALL HOSTS ALL WEAK ASSOCIATE NETWORKS
  find_sites_for_betalinkr(pathogen_transient_network)
  
  ## enter the site matrices printed above 
  transient_pathogen_betalink <- betalinkr_multi(webarray = webs2array(CH, HM, JC, MM, PL, SM, SC, RP),
                                             partitioning="commondenom", binary=FALSE, distofempty='zero', partition.st=TRUE, partition.rr=FALSE)
  
  transient_pathogen_betalink_clean <- fix_betalinkr_output(transient_pathogen_betalink)
}

## **********************************************************
## Run or load turnover by geo distance models
## **********************************************************

## All hosts, both definitions of strong associated microbes included
## only need to run models once, otherwise will load models

if(hosts=="All"){
  run_all_turnover_mods(run.mods=FALSE, # TRUE if never ran model before, false if you want to load models
                        ob.net=obligate_poll_betalink_clean, # Null by default, if run.mods==TRUE input obligate network here
                        trans.net=transient_poll_betalink_clean, # Null by default, if run.mods==TRUE input transient network here
                        filepath="../../../skyIslands/analysis/microbiome/saved/turnover_mods_allhosts.Rdata" # if run.mods=TRUE, input desired save filepath, otherwise input the filepath to load model results
  )
}

## Social hosts, social strong associated microbes included
if(hosts=="Social"){
  run_all_turnover_mods(run.mods=FALSE, # TRUE if never ran model before, false if you want to load models
                        ob.net=obligate_social_betalink_clean, # Null by default, if run.mods==TRUE input obligate network here
                        trans.net=transient_social_betalink_clean, # Null by default, if run.mods==TRUE input transient network here
                        filepath="../../../skyIslands/analysis/microbiome/saved/turnover_mods_social.Rdata" # if run.mods=TRUE, input desired save filepath, otherwise input the filepath to load model results
  )
}

## Social hosts, social strong associated microbes included
if(hosts=="Solitary"){
  run_all_turnover_mods(run.mods=FALSE, # TRUE if never ran model before, false if you want to load models
                        ob.net=obligate_solitary_betalink_clean, # Null by default, if run.mods==TRUE input obligate network here
                        trans.net=transient_solitary_betalink_clean, # Null by default, if run.mods==TRUE input transient network here
                        filepath="../../../skyIslands/analysis/microbiome/saved/turnover_mods_solitary.Rdata" # if run.mods=TRUE, input desired save filepath, otherwise input the filepath to load model results
  )
}

if(hosts=="AllPathogens"){
  run_all_turnover_mods(run.mods=TRUE, # TRUE if never ran model before, false if you want to load models
                        ob.net=obligate_pathogen_betalink_clean, # Null by default, if run.mods==TRUE input obligate network here
                        trans.net=NULL, # Null by default, if run.mods==TRUE input transient network here
                        filepath="../../../skyIslands/analysis/microbiome/saved/turnover_mods_allpathogens.Rdata" # if run.mods=TRUE, input desired save filepath, otherwise input the filepath to load model results
  )
}

## Pairwise bray curtis distance decay models

## prep microbe weights
spec.net <- prepMicrobeWeights(spec.net)

run.decay.mictype.mods=FALSE

if (hosts=="All"){
  if (run.decay.mictype.mods == TRUE){
    #load("../../../skyIslands/data/spec_RBCL_16s.Rdata")
    meta_cols <- c('UniqueID', 'Genus', 'Species', 'GenusSpecies', 'Site', 'Lat', 'Long',"WeightsObligateMicrobe", "WeightsTransientMicrobe")
    
    spec16s <- spec.net %>%
      filter(Apidae == 1) %>%
      select(all_of(meta_cols), starts_with('16s')) %>%
      na.omit()
    
    ob_model <- microbe_type_decay_model(spec16s, 'ObligateAll', model.type = 'exp')
    trans_model <- microbe_type_decay_model(spec16s, 'TransientAll', model.type='exp')
    ## save out models
    save(ob_model,
         trans_model,
         file="../../../skyIslands/analysis/microbiome/saved/decay_mictype_mods_all.Rdata") ## TODO update filepaths
  } else {
    load("../../../skyIslands/analysis/microbiome/saved/decay_mictype_mods_all.Rdata") ## TODO update filepaths
  }
}

if (hosts=="Social"){
  if (run.decay.mictype.mods == TRUE){
    #load("../../../skyIslands/data/spec_RBCL_16s.Rdata")
    meta_cols <- c('UniqueID', 'Genus', 'Species', 'GenusSpecies', 'Site', 'Lat', 'Long',"WeightsObligateSocialMicrobe", "WeightsTransientSocialMicrobe")
    
    spec16s <- spec.net %>%
      filter(Apidae == 1) %>%
      select(all_of(meta_cols), starts_with('16s')) %>%
      na.omit()
    
    ob_model <- microbe_type_decay_model(spec16s, 'ObligateSocial', model.type = 'exp')
    trans_model <- microbe_type_decay_model(spec16s, 'TransientSocial', model.type='exp')
    ## save out models
    save(ob_model,
         trans_model,
         file="../../../skyIslands/analysis/microbiome/saved/decay_mictype_mods_social.Rdata") ## TODO update filepaths
  } else {
    load("../../../skyIslands/analysis/microbiome/saved/decay_mictype_mods_social.Rdata") ## TODO update filepaths
  }
}

if (hosts=="Solitary"){
  if (run.decay.mictype.mods == TRUE){
    #load("../../../skyIslands/data/spec_RBCL_16s.Rdata")
    meta_cols <- c('UniqueID', 'Genus', 'Species', 'GenusSpecies', 'Site', 'Lat', 'Long',"WeightsObligateSolitaryMicrobe", "WeightsTransientSolitaryMicrobe")
    
    spec16s <- spec.net %>%
      filter(Apidae == 1) %>%
      select(all_of(meta_cols), starts_with('16s')) %>%
      na.omit()
    
    ob_model <- microbe_type_decay_model(spec16s, 'ObligateSolitary', model.type = 'exp')
    trans_model <- microbe_type_decay_model(spec16s, 'TransientSolitary', model.type='exp')
    ## save out models
    save(ob_model,
         trans_model,
         file="../../../skyIslands/analysis/microbiome/saved/decay_mictype_mods_solitary.Rdata") ## TODO update filepaths
  } else {
    load("../../../skyIslands/analysis/microbiome/saved/decay_mictype_mods_solitary.Rdata") ## TODO update filepaths
  }
}

## **********************************************************
## Make combined plots for model results for obligate vs
##  transient networks
## **********************************************************

if (hosts=="All") {
  # microbe type comparison
  altpanelA <- plot_decay_ggplot_combined(ob_model,
                                       trans_model,
                                       mod1color='darkgreen',
                                       mod2color='darkorange',
                                       alpha1=0.01,
                                       alpha2=0.01,
                                       lty1='solid',
                                       lty2='solid',
                                       xlab="Geographic Distance (km)",
                                       ylab='Bray-Curtis Dissimilarity', add.points=TRUE)
  
  altpanelA <- altpanelA + labs(tag="A.")
  
  ## B. Interaction turnover
  int.plot <- plot_network_turnover_mod_compare(mod1=int.obligate.mod,
                                                     mod2=int.transient.mod,
                                                     this.network1=obligate_poll_betalink_clean,
                                                     this.network2=transient_poll_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="WholeNetworkLinks",
                                                     label="Total Interaction Turnover")
  int.plot[[1]]
  
  panelB <- int.plot[[1]] + labs(tag="B.")
  int.table <- int.plot[[2]]
  
  
  ## C. rewiring
  
  rewiring.plot <- plot_network_turnover_mod_compare(mod1=rewiring.obligate.mod,
                                                mod2=rewiring.transient.mod,
                                                this.network1=obligate_poll_betalink_clean,
                                                this.network2=transient_poll_betalink_clean,
                                                network_type1='Obligate',
                                                network_type2='Transient',
                                                this.effect="GeoDist",
                                                this.resp="OnlySharedLinks",
                                                label="Rewiring")
  rewiring.plot[[1]]
  
  panelC <- rewiring.plot[[1]] + labs(tag="C.")
  rewiring.table <- rewiring.plot[[2]]
  
  ## D. Host-driven turnover
  
  host.driven.plot <- plot_network_turnover_mod_compare(mod1=host.driven.obligate.mod,
                                                     mod2=host.driven.transient.mod,
                                                     this.network1=obligate_poll_betalink_clean,
                                                     this.network2=transient_poll_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="TurnoverAbsencePollinators",
                                                     label="Host-Driven Turnover")
  host.driven.plot[[1]]
  panelD <- host.driven.plot[[1]] + labs(tag="D.")
  host.table <- host.driven.plot[[2]]
  
  ## E. Microbe-driven turnover
  
  microbe.driven.plot <- plot_network_turnover_mod_compare(mod1=microbe.driven.obligate.mod,
                                                        mod2=microbe.driven.transient.mod,
                                                        this.network1=obligate_poll_betalink_clean,
                                                        this.network2=transient_poll_betalink_clean,
                                                        network_type1='Obligate',
                                                        network_type2='Transient',
                                                        this.effect="GeoDist",
                                                        this.resp="TurnoverAbsenceMicrobes",
                                                        label="Microbe-Driven Turnover")
  microbe.driven.plot[[1]]
  panelE <- microbe.driven.plot[[1]] + labs(tag="E.")
  microbe.table <- microbe.driven.plot[[2]]
  
  ## F. Complete turnover
  
  complete.plot <- plot_network_turnover_mod_compare(mod1=complete.obligate.mod,
                                                           mod2=complete.transient.mod,
                                                           this.network1=obligate_poll_betalink_clean,
                                                           this.network2=transient_poll_betalink_clean,
                                                           network_type1='Obligate',
                                                           network_type2='Transient',
                                                           this.effect="GeoDist",
                                                           this.resp="TurnoverAbsenceBoth",
                                                           label="Complete Turnover")
  
  complete.plot[[1]]
  panelF <- complete.plot[[1]] + labs(tag="F.")
  complete.table <- complete.plot[[2]]
  
  
  # Arrange all panels in the PDF output
  pdf("figures/turnover_combined_all.pdf", width = 8.5, height = 11)  
  grid.arrange(
      altpanelA,
      panelB,
      panelC,
      panelD,
      panelE,
      panelF,
      ncol = 2
  )
  dev.off()
  
  ## Combine results tables and save out
  combined.table <- bind_rows(int.table,
                              rewiring.table,
                              host.table,
                              microbe.table,
                              complete.table)
  
  write.csv(combined.table,
            file=sprintf("saved/tables/turnover_all.csv")) 
}


if (hosts=="Social") {
  # microbe type comparison
  altpanelA <- plot_decay_ggplot_combined(ob_model,
                                          trans_model,
                                          mod1color='darkgreen',
                                          mod2color='darkorange',
                                          alpha1=0.01,
                                          alpha2=0.01,
                                          lty1='solid',
                                          lty2='solid',
                                          xlab="Geographic Distance (km)",
                                          ylab='Bray-Curtis Dissimilarity', add.points=TRUE)
  
  altpanelA <- altpanelA + labs(tag="A.")
  
  ## B. Interaction turnover
  int.plot <- plot_network_turnover_mod_compare(mod1=int.obligate.mod,
                                                mod2=int.transient.mod,
                                                this.network1=obligate_social_betalink_clean,
                                                this.network2=transient_social_betalink_clean,
                                                network_type1='Obligate',
                                                network_type2='Transient',
                                                this.effect="GeoDist",
                                                this.resp="WholeNetworkLinks",
                                                label="Total Interaction Turnover")
  int.plot[[1]]
  
  panelB <- int.plot[[1]] + labs(tag="B.")
  int.table <- int.plot[[2]]
  
  
  ## C. rewiring
  
  rewiring.plot <- plot_network_turnover_mod_compare(mod1=rewiring.obligate.mod,
                                                     mod2=rewiring.transient.mod,
                                                     this.network1=obligate_social_betalink_clean,
                                                     this.network2=transient_social_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="OnlySharedLinks",
                                                     label="Rewiring")
  rewiring.plot[[1]]
  
  panelC <- rewiring.plot[[1]] + labs(tag="C.")
  rewiring.table <- rewiring.plot[[2]]
  
  ## D. Host-driven turnover
  
  host.driven.plot <- plot_network_turnover_mod_compare(mod1=host.driven.obligate.mod,
                                                        mod2=host.driven.transient.mod,
                                                        this.network1=obligate_social_betalink_clean,
                                                        this.network2=transient_social_betalink_clean,
                                                        network_type1='Obligate',
                                                        network_type2='Transient',
                                                        this.effect="GeoDist",
                                                        this.resp="TurnoverAbsencePollinators",
                                                        label="Host-Driven Turnover")
  host.driven.plot[[1]]
  panelD <- host.driven.plot[[1]] + labs(tag="D.")
  host.table <- host.driven.plot[[2]]
  
  ## E. Microbe-driven turnover
  
  microbe.driven.plot <- plot_network_turnover_mod_compare(mod1=microbe.driven.obligate.mod,
                                                           mod2=microbe.driven.transient.mod,
                                                           this.network1=obligate_social_betalink_clean,
                                                           this.network2=transient_social_betalink_clean,
                                                           network_type1='Obligate',
                                                           network_type2='Transient',
                                                           this.effect="GeoDist",
                                                           this.resp="TurnoverAbsenceMicrobes",
                                                           label="Microbe-Driven Turnover")
  microbe.driven.plot[[1]]
  panelE <- microbe.driven.plot[[1]] + labs(tag="E.")
  microbe.table <- microbe.driven.plot[[2]]
  
  ## F. Complete turnover
  
  complete.plot <- plot_network_turnover_mod_compare(mod1=complete.obligate.mod,
                                                     mod2=complete.transient.mod,
                                                     this.network1=obligate_social_betalink_clean,
                                                     this.network2=transient_social_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="TurnoverAbsenceBoth",
                                                     label="Complete Turnover")
  
  complete.plot[[1]]
  panelF <- complete.plot[[1]] + labs(tag="F.")
  complete.table <- complete.plot[[2]]
  
  
  # Arrange all panels in the PDF output
  pdf("figures/turnover_combined_social.pdf", width = 8.5, height = 11)  
  grid.arrange(
    altpanelA,
    panelB,
    panelC,
    panelD,
    panelE,
    panelF,
    ncol = 2
  )
  dev.off()
  
  ## Combine results tables and save out
  combined.table <- bind_rows(int.table,
                              rewiring.table,
                              host.table,
                              microbe.table,
                              complete.table)
  
  write.csv(combined.table,
            file=sprintf("saved/tables/turnover_social.csv")) 
}


if (hosts=="Solitary") {
  # microbe type comparison
  altpanelA <- plot_decay_ggplot_combined(ob_model,
                                          trans_model,
                                          mod1color='darkgreen',
                                          mod2color='darkorange',
                                          alpha1=0.01,
                                          alpha2=0.01,
                                          lty1='solid',
                                          lty2='solid',
                                          xlab="Geographic Distance (km)",
                                          ylab='Bray-Curtis Dissimilarity', add.points=TRUE)
  
  altpanelA <- altpanelA + labs(tag="A.")
  
  ## B. Interaction turnover
  int.plot <- plot_network_turnover_mod_compare(mod1=int.obligate.mod,
                                                mod2=int.transient.mod,
                                                this.network1=obligate_solitary_betalink_clean,
                                                this.network2=transient_solitary_betalink_clean,
                                                network_type1='Obligate',
                                                network_type2='Transient',
                                                this.effect="GeoDist",
                                                this.resp="WholeNetworkLinks",
                                                label="Total Interaction Turnover")
  int.plot[[1]]
  
  panelB <- int.plot[[1]] + labs(tag="B.")
  int.table <- int.plot[[2]]
  
  
  ## C. rewiring
  
  rewiring.plot <- plot_network_turnover_mod_compare(mod1=rewiring.obligate.mod,
                                                     mod2=rewiring.transient.mod,
                                                     this.network1=obligate_solitary_betalink_clean,
                                                     this.network2=transient_solitary_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="OnlySharedLinks",
                                                     label="Rewiring")
  rewiring.plot[[1]]
  
  panelC <- rewiring.plot[[1]] + labs(tag="C.")
  rewiring.table <- rewiring.plot[[2]]
  
  ## D. Host-driven turnover
  
  host.driven.plot <- plot_network_turnover_mod_compare(mod1=host.driven.obligate.mod,
                                                        mod2=host.driven.transient.mod,
                                                        this.network1=obligate_solitary_betalink_clean,
                                                        this.network2=transient_solitary_betalink_clean,
                                                        network_type1='Obligate',
                                                        network_type2='Transient',
                                                        this.effect="GeoDist",
                                                        this.resp="TurnoverAbsencePollinators",
                                                        label="Host-Driven Turnover")
  host.driven.plot[[1]]
  panelD <- host.driven.plot[[1]] + labs(tag="D.")
  host.table <- host.driven.plot[[2]]
  
  ## E. Microbe-driven turnover
  
  microbe.driven.plot <- plot_network_turnover_mod_compare(mod1=microbe.driven.obligate.mod,
                                                           mod2=microbe.driven.transient.mod,
                                                           this.network1=obligate_solitary_betalink_clean,
                                                           this.network2=transient_solitary_betalink_clean,
                                                           network_type1='Obligate',
                                                           network_type2='Transient',
                                                           this.effect="GeoDist",
                                                           this.resp="TurnoverAbsenceMicrobes",
                                                           label="Microbe-Driven Turnover")
  microbe.driven.plot[[1]]
  panelE <- microbe.driven.plot[[1]] + labs(tag="E.")
  microbe.table <- microbe.driven.plot[[2]]
  
  ## F. Complete turnover
  
  complete.plot <- plot_network_turnover_mod_compare(mod1=complete.obligate.mod,
                                                     mod2=complete.transient.mod,
                                                     this.network1=obligate_solitary_betalink_clean,
                                                     this.network2=transient_solitary_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="TurnoverAbsenceBoth",
                                                     label="Complete Turnover")
  
  complete.plot[[1]]
  panelF <- complete.plot[[1]] + labs(tag="F.")
  complete.table <- complete.plot[[2]]
  
  
  # Arrange all panels in the PDF output
  pdf("figures/turnover_combined_solitary.pdf", width = 8.5, height = 11)  
  grid.arrange(
    altpanelA,
    panelB,
    panelC,
    panelD,
    panelE,
    panelF,
    ncol = 2
  )
  dev.off()
  
  ## Combine results tables and save out
  combined.table <- bind_rows(int.table,
                              rewiring.table,
                              host.table,
                              microbe.table,
                              complete.table)
  
  write.csv(combined.table,
            file=sprintf("saved/tables/turnover_solitary.csv")) 
}


