#rm(list=ls())

## packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtreeExtra")
library(ggtreeExtra)
BiocManager::install("phyloseq")
library(phyloseq)

BiocManager::install("SparseArray")
library(SparseArray)



library(tidyr)
library(dplyr)
library(bipartite)
library(phyloseq)
library(TreeTools)
library(devtools)
library(ape)


BiocManager::install("TreeSummarizedExperiment")
library(TreeSummarizedExperiment)

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)


library(ggtree)
library(tidyverse)
library(ggtreeExtra)
library(ggstar)
library(ggplot2)
library(treeio)
library(ggnewscale)
library(tibble)
library(pals)
library(viridis)
library(phyloseq)
library(randomcoloR)
library(phytools)


meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species',
               'GenusSpecies', 'Sex', 'Site', 'Meadow')


## can probs figure out a way to do this automatically but was getting
## frustrated so hard coded this
bombus_sites <- c('JC', 'SM', 'SC', 'MM', 'HM', 'PL', 'CH')
melissodes_sites <- c('JC', 'SC', 'MM', 'HM', 'PL', 'CH', 'RP')

meta_cols <- c('UniqueID', 'Family', 'Genus', 'Species', 'GenusSpecies', 'Sex', 'Site', 'Meadow')

meta <- spec.net %>%
  filter(Apidae == 1) %>%
  select(all_of(meta_cols), Apidae, starts_with('16s')) %>%
  na.omit() %>%
  select(!starts_with('16s')) 

comm_presabs <- as.data.frame(indiv.comm.16sR0) #load in the pres/abs table
comm_presabs[comm_presabs > 0] <- 1 #change all rel abund to 1
comm_presabs <- tibble::rownames_to_column(comm_presabs, "UniqueID") #make rownames (UniqueID) into column

finalASV <- as.data.frame(finalASVtable)
finalASV[finalASV > 0] <- 1 #change all rel abund to 1
finalASV <- tibble::rownames_to_column(finalASV, "UniqueID") #make rownames (UniqueID) into column
