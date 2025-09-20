
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

#load("../../microbeBiogeographyData.Rdata") 

## set hosts="All" to run models for full host dataset with the full list of solitary and social strong HAMs
## set hosts="Social" to run mods for social host dataset with social strong HAMS
## set hosts="Solitary" to run mods for solitary host dataset with solitary strong HAMS
## set hosts="AllPathogens" to run models for full host dataset with the pathogenic microbes
## set hosts="SocialPathogens" to run models for social host dataset with the pathogenic microbes
## set hosts="SolitaryPathogens" to run models for social host dataset with the pathogenic microbes

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

## Full host community, all pathogens
if(hosts=="SocialPathogens"){
  pathogen_obligate_network <- prep_obligate_network(raw_network=spNet_micro, 
                                                     these_obligates=pathogens,
                                                     genera_to_keep=c("Bombus", "Apis")
  )
  
  pathogen_transient_network <- prep_transient_network(raw_network=spNet_micro,
                                                       these_obligates=pathogens,
                                                       genera_to_keep=c("Bombus", "Apis")
  )
}

## Full host community, all pathogens
if(hosts=="SolitaryPathogens"){
  pathogen_obligate_network <- prep_obligate_network(raw_network=spNet_micro, 
                                                     these_obligates=pathogens,
                                                     genera_to_keep=c("Melissodes", "Megachile", "Anthophora", "Andrena")
  )
  
  pathogen_transient_network <- prep_transient_network(raw_network=spNet_micro,
                                                       these_obligates=pathogens,
                                                       genera_to_keep=c("Melissodes", "Megachile", "Anthophora", "Andrena")
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

if(hosts=='SocialPathogens'){
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

if(hosts=='SolitaryPathogens'){
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
                        filepath="C:/Users/rah10/University of Oregon Dropbox/Rebecca Hayes/skyIslands/analysis/microbiome/saved/turnover_mods_social.Rdata" # if run.mods=TRUE, input desired save filepath, otherwise input the filepath to load model results
  )
}

## Social hosts, social strong associated microbes included
if(hosts=="Solitary"){
  run_all_turnover_mods(run.mods=FALSE, # TRUE if never ran model before, false if you want to load models
                        ob.net=obligate_solitary_betalink_clean, # Null by default, if run.mods==TRUE input obligate network here
                        trans.net=transient_solitary_betalink_clean, # Null by default, if run.mods==TRUE input transient network here
                        filepath="C:/Users/rah10/University of Oregon Dropbox/Rebecca Hayes/skyIslands/analysis/microbiome/saved/turnover_mods_solitary.Rdata" # if run.mods=TRUE, input desired save filepath, otherwise input the filepath to load model results
  )
}

if(hosts=="AllPathogens"){
  run_all_turnover_mods(run.mods=FALSE, # TRUE if never ran model before, false if you want to load models
                        ob.net=obligate_pathogen_betalink_clean, # Null by default, if run.mods==TRUE input obligate network here
                        trans.net=transient_pathogen_betalink_clean, # Null by default, if run.mods==TRUE input transient network here
                        filepath="C:/Users/rah10/University of Oregon Dropbox/Rebecca Hayes/skyIslands/analysis/microbiome/saved/turnover_mods_allpathogens.Rdata" # if run.mods=TRUE, input desired save filepath, otherwise input the filepath to load model results
  )
}

if(hosts=="SocialPathogens"){
  run_all_turnover_mods(run.mods=FALSE, # TRUE if never ran model before, false if you want to load models
                        ob.net=obligate_pathogen_betalink_clean, # Null by default, if run.mods==TRUE input obligate network here
                        trans.net=transient_pathogen_betalink_clean, # Null by default, if run.mods==TRUE input transient network here
                        filepath="C:/Users/rah10/University of Oregon Dropbox/Rebecca Hayes/skyIslands/analysis/microbiome/saved/turnover_mods_social_pathogens.Rdata" # if run.mods=TRUE, input desired save filepath, otherwise input the filepath to load model results
  )
}

if(hosts=="SolitaryPathogens"){
  run_all_turnover_mods(run.mods=FALSE, # TRUE if never ran model before, false if you want to load models
                        ob.net=obligate_pathogen_betalink_clean, # Null by default, if run.mods==TRUE input obligate network here
                        trans.net=transient_pathogen_betalink_clean, # Null by default, if run.mods==TRUE input transient network here
                        filepath="C:/Users/rah10/University of Oregon Dropbox/Rebecca Hayes/skyIslands/analysis/microbiome/saved/turnover_mods_solitary_pathogens.Rdata" # if run.mods=TRUE, input desired save filepath, otherwise input the filepath to load model results
  )
}

## **********************************************************
## Make combined plots for model results for obligate vs
##  transient networks
## **********************************************************

if (hosts=="All") {
  
  ## A. Species turnover
  speccomp.plot <- plot_network_turnover_mod_compare(mod1=speccomp.obligate.mod,
                                                mod2=speccomp.transient.mod,
                                                this.network1=obligate_poll_betalink_clean,
                                                this.network2=transient_poll_betalink_clean,
                                                network_type1='Obligate',
                                                network_type2='Transient',
                                                this.effect="GeoDist",
                                                this.resp="DissimilaritySpeciesComposition",
                                                label="Total Composition Turnover")
  speccomp.plot[[1]]
  
  panelA <- speccomp.plot[[1]] + labs(tag="A.")
  speccomp.table <- speccomp.plot[[2]]
  
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
      panelA,
      panelB,
      panelC,
      panelD,
      panelE,
      panelF,
      ncol = 2
  )
  dev.off()
  
  ## Combine results tables and save out
  combined.table <- bind_rows(speccomp.table,
                              int.table,
                              rewiring.table,
                              host.table,
                              microbe.table,
                              complete.table)
  
  write.csv(combined.table,
            file=sprintf("saved/tables/turnover_all.csv")) 
}

if (hosts=="Social") {
 
  # A. Species turnover
  speccomp.plot <- plot_network_turnover_mod_compare(mod1=speccomp.obligate.mod,
                                                     mod2=speccomp.transient.mod,
                                                     this.network1=obligate_social_betalink_clean,
                                                     this.network2=transient_social_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="DissimilaritySpeciesComposition",
                                                     label="Total Composition Turnover")
  speccomp.plot[[1]]

  #panelA <- speccomp.plot[[1]] + labs(tag="A.")
  speccomp.table <- speccomp.plot[[2]]
  
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
  
  #social_panelC <- int.plot[[1]] + labs(tag="C.")
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
  
  panelA <- rewiring.plot[[1]] + labs(tag="A")
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
  panelC <- host.driven.plot[[1]] + labs(tag="C")
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
  panelD <- microbe.driven.plot[[1]] + labs(tag="D")
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
  panelB <- complete.plot[[1]] + labs(tag="B")
  complete.table <- complete.plot[[2]]
  
  
  # Arrange all panels in the PDF output
  pdf("figures/turnover_combined_social.pdf", width = 7, height = 7)  
  grid.arrange(
    #panelA,
    #panelB,
    panelA,
    panelB,
    panelC,
    panelD,
    ncol = 2
  )
  dev.off()
  
  ## Combine results tables and save out
  combined.table <- bind_rows(speccomp.table,
                              int.table,
                              rewiring.table,
                              host.table,
                              microbe.table,
                              complete.table)
  
  write.csv(combined.table,
            file=sprintf("saved/tables/turnover_social_similarity.csv")) 
}

if (hosts=="Solitary") {

  
  ## A. Species turnover
  speccomp.plot <- plot_network_turnover_mod_compare(mod1=speccomp.obligate.mod,
                                                     mod2=speccomp.transient.mod,
                                                     this.network1=obligate_solitary_betalink_clean,
                                                     this.network2=transient_solitary_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="DissimilaritySpeciesComposition",
                                                     label="Total Composition Turnover")
  speccomp.plot[[1]]
  
  #panelA <- speccomp.plot[[1]] + labs(tag="A.")
  speccomp.table <- speccomp.plot[[2]]
  
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
  
  #solitary_panelD <- int.plot[[1]] + labs(tag="D.")
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
  
  panelA <- rewiring.plot[[1]] + labs(tag="A")
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
  panelC <- host.driven.plot[[1]] + labs(tag="C")
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
  panelD <- microbe.driven.plot[[1]] + labs(tag="D")
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
  panelB <- complete.plot[[1]] + labs(tag="B")
  complete.table <- complete.plot[[2]]
  
  
  # Arrange all panels in the PDF output
  pdf("figures/turnover_combined_solitary.pdf", width = 7, height = 7)  
  grid.arrange(
    #panelA,
    #panelB,
    panelA,
    panelB,
    panelC,
    panelD,
    ncol = 2
  )
  dev.off()
  
  ## Combine results tables and save out
  combined.table <- bind_rows(speccomp.table,
                              int.table,
                              rewiring.table,
                              host.table,
                              microbe.table,
                              complete.table)
  
  write.csv(combined.table,
            file=sprintf("saved/tables/turnover_solitary_similarity.csv")) 
}

if (hosts=="AllPathogens") {

  ## B. Interaction turnover
  int.plot <- plot_network_turnover_mod_compare(mod1=int.obligate.mod,
                                                mod2=int.transient.mod,
                                                this.network1=obligate_pathogen_betalink_clean,
                                                this.network2=transient_pathogen_betalink_clean,
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
                                                     this.network1=obligate_pathogen_betalink_clean,
                                                     this.network2=transient_pathogen_betalink_clean,
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
                                                        this.network1=obligate_pathogen_betalink_clean,
                                                        this.network2=transient_pathogen_betalink_clean,
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
                                                           this.network1=obligate_pathogen_betalink_clean,
                                                           this.network2=transient_pathogen_betalink_clean,
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
                                                     this.network1=obligate_pathogen_betalink_clean,
                                                     this.network2=transient_pathogen_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="TurnoverAbsenceBoth",
                                                     label="Complete Turnover")
  
  complete.plot[[1]]
  panelF <- complete.plot[[1]] + labs(tag="F.")
  complete.table <- complete.plot[[2]]
  
  
  # Arrange all panels in the PDF output
  pdf("figures/turnover_combined_pathogens.pdf", width = 8.5, height = 11)  
  grid.arrange(
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
            file=sprintf("saved/tables/turnover_pathogens.csv")) 
}

if (hosts=="SolitaryPathogens") {

  
  ## A. Species turnover
  speccomp.plot <- plot_network_turnover_mod_compare(mod1=speccomp.obligate.mod,
                                                     mod2=speccomp.transient.mod,
                                                     this.network1=obligate_pathogen_betalink_clean,
                                                     this.network2=transient_pathogen_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="DissimilaritySpeciesComposition",
                                                     label="Total Composition Turnover")
  speccomp.plot[[1]]
  
  #panelA <- speccomp.plot[[1]] + labs(tag="A.")
  speccomp.table <- speccomp.plot[[2]]
  
  ## B. Interaction turnover
  int.plot <- plot_network_turnover_mod_compare(mod1=int.obligate.mod,
                                                mod2=int.transient.mod,
                                                this.network1=obligate_pathogen_betalink_clean,
                                                this.network2=transient_pathogen_betalink_clean,
                                                network_type1='Obligate',
                                                network_type2='Transient',
                                                this.effect="GeoDist",
                                                this.resp="WholeNetworkLinks",
                                                label="Total Interaction Turnover")
  int.plot[[1]]
  
  #solitary_panelD <- int.plot[[1]] + labs(tag="D.")
  int.table <- int.plot[[2]]
  
  
  ## C. rewiring
  
  rewiring.plot <- plot_network_turnover_mod_compare(mod1=rewiring.obligate.mod,
                                                     mod2=rewiring.transient.mod,
                                                     this.network1=obligate_pathogen_betalink_clean,
                                                     this.network2=transient_pathogen_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="OnlySharedLinks",
                                                     label="Rewiring")
  rewiring.plot[[1]]
  
  panelA <- rewiring.plot[[1]] + labs(tag="A")
  rewiring.table <- rewiring.plot[[2]]
  
  ## D. Host-driven turnover
  
  host.driven.plot <- plot_network_turnover_mod_compare(mod1=host.driven.obligate.mod,
                                                        mod2=host.driven.transient.mod,
                                                        this.network1=obligate_pathogen_betalink_clean,
                                                        this.network2=transient_pathogen_betalink_clean,
                                                        network_type1='Obligate',
                                                        network_type2='Transient',
                                                        this.effect="GeoDist",
                                                        this.resp="TurnoverAbsencePollinators",
                                                        label="Host-Driven Turnover")
  host.driven.plot[[1]]
  panelC <- host.driven.plot[[1]] + labs(tag="C")
  host.table <- host.driven.plot[[2]]
  
  ## E. Microbe-driven turnover
  
  microbe.driven.plot <- plot_network_turnover_mod_compare(mod1=microbe.driven.obligate.mod,
                                                           mod2=microbe.driven.transient.mod,
                                                           this.network1=obligate_pathogen_betalink_clean,
                                                           this.network2=transient_pathogen_betalink_clean,
                                                           network_type1='Obligate',
                                                           network_type2='Transient',
                                                           this.effect="GeoDist",
                                                           this.resp="TurnoverAbsenceMicrobes",
                                                           label="Microbe-Driven Turnover")
  microbe.driven.plot[[1]]
  panelD <- microbe.driven.plot[[1]] + labs(tag="D")
  microbe.table <- microbe.driven.plot[[2]]
  
  ## F. Complete turnover
  
  complete.plot <- plot_network_turnover_mod_compare(mod1=complete.obligate.mod,
                                                     mod2=complete.transient.mod,
                                                     this.network1=obligate_pathogen_betalink_clean,
                                                     this.network2=transient_pathogen_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="TurnoverAbsenceBoth",
                                                     label="Complete Turnover")
  
  complete.plot[[1]]
  panelB <- complete.plot[[1]] + labs(tag="B")
  complete.table <- complete.plot[[2]]
  
  
  # Arrange all panels in the PDF output
  pdf("figures/turnover_combined_solitary_pathogens.pdf", width = 7, height = 7)  
  grid.arrange(
    #panelA,
    #panelB,
    panelA,
    panelB,
    panelC,
    panelD,
    ncol = 2
  )
  dev.off()
  
  ## Combine results tables and save out
  combined.table <- bind_rows(speccomp.table,
                              int.table,
                              rewiring.table,
                              host.table,
                              microbe.table,
                              complete.table)
  
  write.csv(combined.table,
            file=sprintf("saved/tables/turnover_solitary_pathogens.csv")) 
}

if (hosts=="SocialPathogens") {
  
  # A. Species turnover
  speccomp.plot <- plot_network_turnover_mod_compare(mod1=speccomp.obligate.mod,
                                                     mod2=speccomp.transient.mod,
                                                     this.network1=obligate_pathogen_betalink_clean,
                                                     this.network2=transient_pathogen_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="DissimilaritySpeciesComposition",
                                                     label="Total Composition Turnover")
  speccomp.plot[[1]]
  
  #panelA <- speccomp.plot[[1]] + labs(tag="A.")
  speccomp.table <- speccomp.plot[[2]]
  
  ## B. Interaction turnover
  int.plot <- plot_network_turnover_mod_compare(mod1=int.obligate.mod,
                                                mod2=int.transient.mod,
                                                this.network1=obligate_pathogen_betalink_clean,
                                                this.network2=transient_pathogen_betalink_clean,
                                                network_type1='Obligate',
                                                network_type2='Transient',
                                                this.effect="GeoDist",
                                                this.resp="WholeNetworkLinks",
                                                label="Total Interaction Turnover")
  int.plot[[1]]
  
  #social_panelC <- int.plot[[1]] + labs(tag="C.")
  int.table <- int.plot[[2]]
  
  
  ## C. rewiring
  
  rewiring.plot <- plot_network_turnover_mod_compare(mod1=rewiring.obligate.mod,
                                                     mod2=rewiring.transient.mod,
                                                     this.network1=obligate_pathogen_betalink_clean,
                                                     this.network2=transient_pathogen_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="OnlySharedLinks",
                                                     label="Rewiring")
  rewiring.plot[[1]]
  
  panelA <- rewiring.plot[[1]] + labs(tag="A")
  rewiring.table <- rewiring.plot[[2]]
  
  ## D. Host-driven turnover
  
  host.driven.plot <- plot_network_turnover_mod_compare(mod1=host.driven.obligate.mod,
                                                        mod2=host.driven.transient.mod,
                                                        this.network1=obligate_pathogen_betalink_clean,
                                                        this.network2=transient_pathogen_betalink_clean,
                                                        network_type1='Obligate',
                                                        network_type2='Transient',
                                                        this.effect="GeoDist",
                                                        this.resp="TurnoverAbsencePollinators",
                                                        label="Host-Driven Turnover")
  host.driven.plot[[1]]
  panelC <- host.driven.plot[[1]] + labs(tag="C")
  host.table <- host.driven.plot[[2]]
  
  ## E. Microbe-driven turnover
  
  microbe.driven.plot <- plot_network_turnover_mod_compare(mod1=microbe.driven.obligate.mod,
                                                           mod2=microbe.driven.transient.mod,
                                                           this.network1=obligate_pathogen_betalink_clean,
                                                           this.network2=transient_pathogen_betalink_clean,
                                                           network_type1='Obligate',
                                                           network_type2='Transient',
                                                           this.effect="GeoDist",
                                                           this.resp="TurnoverAbsenceMicrobes",
                                                           label="Microbe-Driven Turnover")
  microbe.driven.plot[[1]]
  panelD <- microbe.driven.plot[[1]] + labs(tag="D")
  microbe.table <- microbe.driven.plot[[2]]
  
  ## F. Complete turnover
  
  complete.plot <- plot_network_turnover_mod_compare(mod1=complete.obligate.mod,
                                                     mod2=complete.transient.mod,
                                                     this.network1=obligate_pathogen_betalink_clean,
                                                     this.network2=transient_pathogen_betalink_clean,
                                                     network_type1='Obligate',
                                                     network_type2='Transient',
                                                     this.effect="GeoDist",
                                                     this.resp="TurnoverAbsenceBoth",
                                                     label="Complete Turnover")
  
  complete.plot[[1]]
  panelB <- complete.plot[[1]] + labs(tag="B")
  complete.table <- complete.plot[[2]]
  
  
  # Arrange all panels in the PDF output
  pdf("figures/turnover_combined_solitary_pathogens.pdf", width = 7, height = 7)  
  grid.arrange(
    #panelA,
    #panelB,
    panelA,
    panelB,
    panelC,
    panelD,
    ncol = 2
  )
  dev.off()
  
  ## Combine results tables and save out
  combined.table <- bind_rows(speccomp.table,
                              int.table,
                              rewiring.table,
                              host.table,
                              microbe.table,
                              complete.table)
  
  write.csv(combined.table,
            file=sprintf("saved/tables/turnover_solitary_pathogens.csv")) 
}

