library(ggplot2)
library(brms)
library(dplyr)
library(tidyr)
library(tidybayes)
library(ggthemes)
library(bayestestR)
library(gtools)

## plotting
library(viridis)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(rstantools)
library(performance)
library(bayestestR)

save.dir <- "saved/tables"
if(!dir.exists(save.dir)) {
  dir.create(save.dir, showWarnings = FALSE)
}

fig.dir <- "figures"
if(!dir.exists(fig.dir)) {
  dir.create(save.dir, showWarnings = FALSE)
}



load('../../microbeBiogeography_data.Rdata')
load("../../../skyIslands/data/trees.Rdata") #TODO update once ask LP about how to handle this additional data
site.sum <- read.csv("../../../skyIslands/data/sitestats.csv") #TODO update once ask LP about how to handle this additional data


spec.net <- spec.net[!is.na(spec.net$GenusSpecies),]


parasites <- c(#"AspergillusSpp", ## problematic parasite!
  "AscosphaeraSpp",
  "ApicystisSpp",
  "CrithidiaExpoeki",
  "CrithidiaBombi",
  "CrithidiaSpp",
  "NosemaBombi",
  "NosemaCeranae")

# site.sum <- site.sum[,c("Site",
#                         "Year",
#                         "SampleRound")]
# 
# spec.net <- merge(spec.net, site.sum, all.x=TRUE)
# 
# traits <-
#     read.csv("../../../skyIslands_saved/data/raw/bee_traits.csv") #TODO update once ask LP about how to handle this additional data
# traits$GenusSpecies <- fix.white.space(traits$GenusSpecies)
# traits <- traits[, c("GenusSpecies", "Sociality", "Lecty", "MeanITD"),]
# 
# net.traits <- read.csv("../../../skyIslands/data/networks_traits.csv") #TODO update once ask LP about how to handle this additional data
# net.traits <- net.traits[, c("GenusSpecies", "r.degree"),]
# 
# traits <- merge(traits, net.traits, by="GenusSpecies", all.x=TRUE)
# 
# spec.net <- merge(spec.net, traits, all.x=TRUE, by="GenusSpecies")

dir.create(path="saved", showWarnings = FALSE)
dir.create(path="saved/tables", showWarnings = FALSE)

spec.net <- spec.net[order(spec.net$Site),]

## fixing blank spaces with noID
spec.net$GenusSpecies <- if_else(spec.net$GenusSpecies=='', 'NoID', spec.net$GenusSpecies)

## raw, non standardized data for plotting
spec.orig <- prepDataSEM(spec.net, variables.to.log, 
                         standardize=FALSE)

## Make SEM weights and standardize data.
spec.net <- prepDataSEM(spec.net, variables.to.log, 
                        vars_yearsr = vars_yearsr, vars_sp = vars_sp, 
                        vars_yearsrsp = vars_yearsrsp,
                        vars_site=vars_site)



##load tree from :
##Henriquez Piskulich, Patricia Andrea; Hugall, Andrew F.; Stuart-Fox, Devi (2023). A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila) [Dataset]. Dryad. https://doi.org/10.5061/dryad.80gb5mkw1

phylo <- ape::read.tree("../../../skyIslands/data/BEE_mat7_fulltree.nwk") #TODO update once ask LP about how to handle this additional data


##clean up unwanted portion of labels
pattern <- "(_n\\d+m\\d+_[A-Za-z0-9]+)?$"
phylo$tip.label <- gsub(pattern, "", phylo$tip.label)

## replace underscore with space
phylo$tip.label <- gsub("_", " ", phylo$tip.label)

## Species that are not in the phylogeny are not used. brms is not allowing an incomplete
## phylogeny, to avoid the error we changed the species not present to one that is in the phylogeny. 
## We chose a species for which we did not do parasite screening and should not influence results.
not_in_phylo <- unique(spec.net$GenusSpecies[!spec.net$GenusSpecies %in% phylo$tip.label])
spec.net$GenusSpecies[spec.net$GenusSpecies %in% not_in_phylo]<- "Agapostemon angelicus"



phylo_matrix <- ape::vcv.phylo(phylo)

## check which individuals don't have microbe data
drop.PD.NA <- unique(spec.net$UniqueID[spec.net$WeightsMicrobe == 1 &
                                         is.na(spec.net$PD)])

## drop individuals that had parasite screen but not microbe
## filling in zeros for NAs in PD
spec.net <- spec.net[!(spec.net$UniqueID %in% drop.PD.NA),] %>%
  mutate(PD = ifelse(!is.na(PD), PD, 0))

# ##adding abundance weights column
# abund_csv <- data.frame(read.csv("../../../skyIslands/data/sp_year_site_round.csv")) #TODO update once ask LP about how to handle this additional data
# 
# #join abundance csv
# spec.net <- merge(spec.net, abund_csv)

## what I want:
## for subset -- should be ones and zeroes ones for bombus that had microbe screening
## for weights -- log abundance weights for obligate and transient, but only for that genus 

## make genus specific weights 0s and  bombus
spec.net$BombusWeights = ifelse(spec.net$Genus=='Bombus'&spec.net$WeightsMicrobe==1, 1, 0)

## make genus specific weights 0s and 1s apis
spec.net$ApisWeights = ifelse(spec.net$Genus=='Apis'&spec.net$WeightsMicrobe==1, 1, 0)

## make genus specific weights 0s and 1s melissodes
spec.net$MelissodesWeights = ifelse(spec.net$Genus=='Melissodes'&spec.net$WeightsMicrobe==1, 1, 0)



spec.net$WeightsObligateBombus = spec.net$WeightsObligateMicrobe*spec.net$BombusWeights


spec.net$WeightsTransientBombus = spec.net$WeightsTransientMicrobe*spec.net$BombusWeights

spec.net$WeightsObligateMelissodes = spec.net$WeightsObligateMicrobe*spec.net$MelissodesWeights

spec.net$WeightsTransientMelissodes = spec.net$WeightsTransientMicrobe*spec.net$MelissodesWeights

spec.net$PD.obligate.log <- log(spec.net$PD.obligate + 1)
spec.net$PD.obligate.log <- ifelse(is.na(spec.net$PD.obligate.log), 0, spec.net$PD.obligate.log)

spec.net$PD.transient.log <- log(spec.net$PD.transient + 1)
spec.net$PD.transient.log <- ifelse(is.na(spec.net$PD.transient.log), 0, spec.net$PD.transient.log)


save(spec.net, file="../../spec_microbes.Rdata")
