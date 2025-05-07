
## **********************************************************
## Load libraries and source files
## **********************************************************

## based on the following tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html

rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd('microbeBiogeography/figures')

library(pals)
library(cowplot)
library(ggpubr)
library(gridExtra)
library(grid)
library(ggplot2)
library(dplyr)
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

load("../microbeBiogeographyData.Rdata")

source('src/trees_init.R')
source('src/tree_functions.R')

meta_cols <- c('UniqueID', 'Genus', 'Species', 'GenusSpecies', 'Site', 'Meadow')

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


## **********************************************************
## Create microbe phylo trees by genus
## **********************************************************

## Bombus tree
bombus_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta=meta, "Bombus", genus.or.spp='Genus', finalASV, bombus_sites, do_collapse = TRUE, add_tip_labs = FALSE)
panelA <- bombus_tree[[1]] + labs(tag="A. Bombus (n=283)")
bombus_meta <- bombus_tree[[2]]
panelA

# Open a PDF device to save the plot
pdf("bombus.pdf",
    height=24, width=33)  # Adjust width and height as needed

# Plot the final combined figure
plot(panelA)

# Close the PDF device to complete saving
dev.off()

## Melissodes tree
melissodes_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes", genus.or.spp='Genus', finalASV, melissodes_sites, do_collapse = TRUE, add_tip_labs = FALSE)
panelD <- melissodes_tree[[1]] + labs(tag="D. Melissodes (n=46)")
melissodes_meta <- melissodes_tree[[2]]
panelD

# Open a PDF device to save the plot
pdf("melissodes.pdf",
    height=24, width=33)  # Adjust width and height as needed

# Plot the final combined figure
plot(panelD)

# Close the PDF device to complete saving
dev.off()


## Apis tree
apis_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis", genus.or.spp='Genus', finalASV, apis_sites, do_collapse = TRUE, add_tip_labs = FALSE)
panelB <- apis_tree[[1]] + labs(tag="B. Apis (n=220)")
apis_meta <- apis_tree[[2]]
panelB

# Open a PDF device to save the plot
pdf("apis.pdf",
    height=24, width=33)  # Adjust width and height as needed

# Plot the final combined figure
plot(panelB)

# Close the PDF device to complete saving
dev.off()

## Megachile tree
megachile_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Megachile", genus.or.spp='Genus', finalASV, megachile_sites, do_collapse = TRUE, add_tip_labs = FALSE)
panelF <- megachile_tree[[1]] + labs(tag="F. Megachile (n=41)")
megachile_meta <- megachile_tree[[2]]
panelF

# Open a PDF device to save the plot
pdf("megachile.pdf",
    height=24, width=33)  # Adjust width and height as needed

# Plot the final combined figure
plot(panelF)

# Close the PDF device to complete saving
dev.off()

## Anthophora tree
anthophora_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Anthophora", genus.or.spp='Genus', finalASV, anthophora_sites, do_collapse = TRUE, add_tip_labs = FALSE)
panelC <- anthophora_tree[[1]] + labs(tag="C. Anthophora (n=34)")
anthophora_meta <- anthophora_tree[[2]]
panelC

# Open a PDF device to save the plot
pdf("anthophora.pdf",
    height=24, width=33)  # Adjust width and height as needed

# Plot the final combined figure
plot(panelC)

# Close the PDF device to complete saving
dev.off()

## Andrena tree
andrena_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Andrena", genus.or.spp='Genus', finalASV, andrena_sites, do_collapse = TRUE, add_tip_labs=FALSE)
panelE <- andrena_tree[[1]] + labs(tag="E. Andrena (n=80)")
andrena_meta <- andrena_tree[[2]]
panelE

## save out andrena tree
# Open a PDF device to save the plot
pdf("andrena.pdf",
    height=24, width=33)  # Adjust width and height as needed

# Plot the final combined figure
plot(panelE)

# Close the PDF device to complete saving
dev.off()

## **********************************************************
## Create custom legend
## **********************************************************
library(ggplot2) 
library(cowplot)  # for get_legend()

# All 15 taxa
these_obligates <- c(
  "Acetobacteraceae", # BothObligate
  "Lactobacillaceae", # BothObligate
  "Bartonella",       # SocialObligate
  "Bifidobacteriaceae",# SocialObligate
  "Neisseriaceae",    # SocialObligate
  "Orbaceae",         # SocialObligate
  "Bacillaceae",      # SolitaryObligate
  "Burkholderiaceae", # SolitaryObligate
  "Clostridiaceae",   # SolitaryObligate
  "Comamonadaceae",   # SolitaryObligate
  "Enterobacteriaceae", # SolitaryObligate
  "Methylobacteriaceae", # SolitaryObligate
  "Moraxellaceae",    # SolitaryObligate
  "Oxalobacteraceae", # SolitaryObligate
  "Sphingomonadaceae" # SolitaryObligate
)

# Matching colors
leg_col_order <- c(
  "#F6A600", # Acetobacteraceae
  "#882D17", # Lactobacillaceae
  "#7C6BD0", # Bartonella
  "#B3446C", # Bifidobacteriaceae
  "#DCD300", # Neisseriaceae
  "#8DB600", # Orbaceae
  "#5AC8FA", # Bacillaceae
  "#E66100", # Burkholderiaceae
  "#3CB371", # Clostridiaceae
  "#5D8AA8", # Comamonadaceae
  "#C83737", # Enterobacteriaceae
  "#A85C90", # Methylobacteriaceae
  "#1F78B4", # Moraxellaceae
  "#D2691E", # Oxalobacteraceae
  "#20B2AA"  # Sphingomonadaceae
)

# Matching shapes
leg_shape_order <- c(
  21, # Acetobacteraceae
  21, # Lactobacillaceae
  22, # Bartonella
  22, # Bifidobacteriaceae
  22, # Neisseriaceae
  22, # Orbaceae
  24, # Bacillaceae
  24, # Burkholderiaceae
  24, # Clostridiaceae
  24, # Comamonadaceae
  24, # Enterobacteriaceae
  24, # Methylobacteriaceae
  24, # Moraxellaceae
  24, # Oxalobacteraceae
  24  # Sphingomonadaceae
)

# Create a DataFrame for legend
data.leg <- data.frame(
  Xdata = rnorm(length(these_obligates)),
  Ydata = rnorm(length(these_obligates)),
  Family = these_obligates,
  leg_color = leg_col_order,
  leg_shape = leg_shape_order
)

# Order the taxa within each shape group alphabetically
data.leg$Family <- factor(data.leg$Family, levels = unique(c(
  # First, order taxa with shape 21 alphabetically
  sort(data.leg$Family[data.leg$leg_shape == 21]), 
  # Then, order taxa with shape 22 alphabetically
  sort(data.leg$Family[data.leg$leg_shape == 22]),
  # Then, order taxa with shape 24 alphabetically
  sort(data.leg$Family[data.leg$leg_shape == 24])
)))

# Create plot for legend
gplot <- ggplot(data.leg, aes(Xdata, Ydata, fill = Family, shape = Family)) +
  geom_point(size = 7, color = "black") +  # Black outline for clarity
  scale_fill_manual(values = setNames(leg_col_order, these_obligates)) +
  scale_shape_manual(values = setNames(leg_shape_order, these_obligates)) +
  theme(
    legend.position = 'bottom',
    legend.box = "vertical",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key = element_blank()
  ) +
  labs(fill = 'Bacteria Family', shape = 'Bacteria Family') +
  guides(
    fill = guide_legend(nrow = 3, ncol = 5),  # Adjust the number of rows and columns
    shape = guide_legend(nrow = 3, ncol = 5)  # Adjust the number of rows and columns
  )

# Draw only the legend
panelG <- get_legend(gplot)
plot(panelG)




## **********************************************************
## Create full paneled figure with legend
## **********************************************************

#Set up the layout matrix so that the bottom legend spans both columns
layout <- rbind(c(1, 2, 3, 4, 5, 6),
                c(7, 7, 7, 7, 7, 7)) # The legend will span both columns

# Create the final layout
final_plot <- arrangeGrob(panelA, panelB, panelC, panelD, panelE, panelF, panelG,
                          layout_matrix = layout,
                          heights = c(9, 1)) # Adjust the heights as needed

# Center the legend panelD properly within its grid
panelG_centered <- arrangeGrob(panelG, 
                               ncol = 1, 
                               padding = unit(1, "lines"))  # Add padding if necessary

# Open a PDF device to save the plot
pdf("grid_trees_all.pdf",
    height=24, width=33)  # Adjust width and height as needed

# Plot the final combined figure
grid.arrange(final_plot)

# Close the PDF device to complete saving
dev.off()
