
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
## **********************************************************
## Prep obligate and transient networks
## **********************************************************

only_obligate_network <- prep_obligate_network(raw_network=spNet_micro)

only_obligate_network_BM <- prep_obligate_network(raw_network=spNet_micro,  
                                                      genera_to_keep=c("Melissodes", "Bombus", "Apis"))


only_transient_network <- prep_transient_network(raw_network=spNet_micro)
only_transient_network_BM <- prep_transient_network(raw_network=spNet_micro,genera_to_keep=c("Melissodes", "Bombus", "Apis"))

## **********************************************************
## Run network betalinkr function and prep output table
## **********************************************************

source("src/betalinkrPrep.R")

## **********************************************************
## Run or load turnover by geo distance models
## **********************************************************

## only need to run models once, otherwise will load models

run.mods=FALSE

if (run.mods==TRUE){
  
## Interaction turnover
int.obligate.mod <- run_network_turnover_mod(this_component="WholeNetworkLinks",
                                                    this_network=obligate_poll_betalink)
  
int.transient.mod <- run_network_turnover_mod(this_component="WholeNetworkLinks",
                                                     this_network=transient_poll_betalink)

## Turnover species composition
speccomp.obligate.mod <- run_network_turnover_mod(this_component="DissimilaritySpeciesComposition",
                                                  this_network=obligate_poll_betalink)

speccomp.transient.mod <- run_network_turnover_mod(this_component="DissimilaritySpeciesComposition",
                                                   this_network=transient_poll_betalink)
  
  
## Rewiring
rewiring.obligate.mod <- run_network_turnover_mod(this_component="OnlySharedLinks",
                                                  this_network=obligate_poll_betalink)

rewiring.transient.mod <- run_network_turnover_mod(this_component="OnlySharedLinks",
                                                   this_network=transient_poll_betalink)

## Host-driven turnover
host.driven.obligate.mod <- run_network_turnover_mod(this_component="TurnoverAbsencePollinators",
                                                     this_network=obligate_poll_betalink)

host.driven.transient.mod <- run_network_turnover_mod(this_component="TurnoverAbsencePollinators",
                                                      this_network=transient_poll_betalink)

## Microbe-driven turnover
microbe.driven.obligate.mod <- run_network_turnover_mod(this_component="TurnoverAbsenceMicrobes",
                                                        this_network=obligate_poll_betalink)

microbe.driven.transient.mod <- run_network_turnover_mod(this_component="TurnoverAbsenceMicrobes",
                                                         this_network=transient_poll_betalink)

## Complete turnover
complete.obligate.mod <- run_network_turnover_mod(this_component="TurnoverAbsenceBoth",
                                                  this_network=obligate_poll_betalink)

complete.transient.mod <- run_network_turnover_mod(this_component="TurnoverAbsenceBoth",
                                                   this_network=transient_poll_betalink)

## save out models
save(int.obligate.mod,
     int.transient.mod,
     speccomp.obligate.mod,
     speccomp.transient.mod,
     rewiring.obligate.mod,
       rewiring.transient.mod,
       host.driven.obligate.mod,
       host.driven.transient.mod,
       microbe.driven.obligate.mod,
       microbe.driven.transient.mod,
       complete.obligate.mod,
       complete.transient.mod,
       file="../../../skyIslands/analysis/microbiome/saved/turnover_mods_BM.Rdata") ## TODO update with correct filepath
} else { 
  load("../../../skyIslands/analysis/microbiome/saved/turnover_mods_BM.Rdata") ## TODO update with correct filepath
}


## Pairwise bray curtis distance decay models
source("src/distDecay.R")

run.decay.genus.mods=FALSE

if (run.decay.genus.mods == TRUE){
  bombus_model <- genusspecies_decay_model(spec16s, 'Bombus', type='Genus', model.type = 'exp')
  melissodes_model <- genusspecies_decay_model(spec16s, 'Melissodes', type='Genus', model.type='exp')
  ## save out models
  save(bombus_model,
       melissodes_model,
       file="../../../skyIslands/analysis/microbiome/saved/decay_genus_mods.Rdata") ## TODO update corrected filepaths
} else {
    load("../../../skyIslands/analysis/microbiome/saved/decay_genus_mods.Rdata") ## TODO update corrected filepaths
}

run.decay.mictype.mods=FALSE

if (run.decay.mictype.mods == TRUE){
  #load("../../../skyIslands/data/spec_RBCL_16s.Rdata")
  meta_cols <- c('UniqueID', 'Genus', 'Species', 'GenusSpecies', 'Site', 'Lat', 'Long',"WeightsObligateMicrobe", "WeightsTransientMicrobe")
  
  spec16s <- spec.net %>%
    filter(Apidae == 1) %>%
    select(all_of(meta_cols), starts_with('16s')) %>%
    na.omit()
  ob_model <- microbe_type_decay_model(spec16s, 'Obligate', model.type = 'exp')
  trans_model <- microbe_type_decay_model(spec16s, 'Facultative', model.type='exp')
  ## save out models
  save(ob_model,
       trans_model,
       file="../../../skyIslands/analysis/microbiome/saved/decay_mictype_mods_BM.Rdata") ## TODO update filepaths
} else {
  load("../../../skyIslands/analysis/microbiome/saved/decay_mictype_mods_BM.Rdata") ## TODO update filepaths
}


## **********************************************************
## Make combined plots for model results for obligate vs
##  transient networks
## **********************************************************


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
                                                   this.network1=obligate_poll_betalink,
                                                   this.network2=transient_poll_betalink,
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
                                              this.network1=obligate_poll_betalink,
                                              this.network2=transient_poll_betalink,
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
                                                   this.network1=obligate_poll_betalink,
                                                   this.network2=transient_poll_betalink,
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
                                                      this.network1=obligate_poll_betalink,
                                                      this.network2=transient_poll_betalink,
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
                                                         this.network1=obligate_poll_betalink,
                                                         this.network2=transient_poll_betalink,
                                                         network_type1='Obligate',
                                                         network_type2='Transient',
                                                         this.effect="GeoDist",
                                                         this.resp="TurnoverAbsenceBoth",
                                                         label="Complete Turnover")

complete.plot[[1]]
panelF <- complete.plot[[1]] + labs(tag="F.")
complete.table <- complete.plot[[2]]


# Arrange all panels in the PDF output
pdf("figures/turnover_combined_BM.pdf", width = 8.5, height = 11)  
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
          file=sprintf("saved/tables/turnover_BM.csv")) 
