
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
panelA <- bombus_tree[[1]] + labs(tag="A Bombus")
bombus_meta <- bombus_tree[[2]]
panelA

# Open a PDF device to save the plot
pdf("figures/bombus.pdf",
    height=24, width=33)  # Adjust width and height as needed

# Plot the final combined figure
plot(panelA)

# Close the PDF device to complete saving
dev.off()

## Melissodes tree
melissodes_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes", genus.or.spp='Genus', finalASV, melissodes_sites, do_collapse = TRUE, add_tip_labs = FALSE)
panelD <- melissodes_tree[[1]] + labs(tag="D Melissodes")
melissodes_meta <- melissodes_tree[[2]]
panelD

# Open a PDF device to save the plot
pdf("figures/melissodes.pdf",
    height=24, width=33)  # Adjust width and height as needed

# Plot the final combined figure
plot(panelD)

# Close the PDF device to complete saving
dev.off()


## Apis tree
apis_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis", genus.or.spp='Genus', finalASV, apis_sites, do_collapse = TRUE, add_tip_labs = FALSE)
panelB <- apis_tree[[1]] + labs(tag="B Apis")
apis_meta <- apis_tree[[2]]
panelB

# Open a PDF device to save the plot
pdf("figures/apis.pdf",
    height=24, width=33)  # Adjust width and height as needed

# Plot the final combined figure
plot(panelB)

# Close the PDF device to complete saving
dev.off()

## Megachile tree
megachile_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Megachile", genus.or.spp='Genus', finalASV, megachile_sites, do_collapse = TRUE, add_tip_labs = FALSE)
panelF <- megachile_tree[[1]] + labs(tag="F Megachile")
megachile_meta <- megachile_tree[[2]]
panelF

# Open a PDF device to save the plot
pdf("figures/megachile.pdf",
    height=24, width=33)  # Adjust width and height as needed

# Plot the final combined figure
plot(panelF)

# Close the PDF device to complete saving
dev.off()

## Anthophora tree
anthophora_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Anthophora", genus.or.spp='Genus', finalASV, anthophora_sites, do_collapse = TRUE, add_tip_labs = FALSE)
panelC <- anthophora_tree[[1]] + labs(tag="C Anthophora")
anthophora_meta <- anthophora_tree[[2]]
panelC

# Open a PDF device to save the plot
pdf("figures/anthophora.pdf",
    height=24, width=33)  # Adjust width and height as needed

# Plot the final combined figure
plot(panelC)

# Close the PDF device to complete saving
dev.off()

## Andrena tree
andrena_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Andrena", genus.or.spp='Genus', finalASV, andrena_sites, do_collapse = TRUE, add_tip_labs=FALSE)
panelE <- andrena_tree[[1]] + labs(tag="E Andrena")
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
library(cowplot) # for get_legend()

# Define all families and their colors
legend_families <- c(
  # BothObligate
  "Acetobacteraceae", "Lactobacillaceae",
  
  # SocialObligate
  "Bartonellaceae", "Bifidobacteriaceae", "Neisseriaceae", "Orbaceae",
  
  # SolitaryObligate
  "Bacillaceae", "Burkholderiaceae", "Clostridiaceae", "Comamonadaceae",
  "Enterobacteriaceae", "Methylobacteriaceae", "Moraxellaceae", "Oxalobacteraceae", "Sphingomonadaceae",
  
  # Pathogens
  "Wolbachia", "Erwinia", "Hafnia"
)

legend_colors <- c(
  # BothObligate
  "#F6A600", "#882D17",
  
  # SocialObligate
  "#7C6BD0", "#B3446C", "#DCD300", "#8DB600",
  
  # SolitaryObligate
  "#5AC8FA", "#E66100", "#3CB371", "#5D8AA8",
  "#C83737", "#A85C90", "#1F78B4", "#D2691E", "#20B2AA",
  
  # Pathogens
  "#E41A1C", "#377EB8", "#4DAF4A"
)

legend_shapes <- c(
  # BothObligate
  21, 21,
  
  # SocialObligate
  22, 22, 22, 22,
  
  # SolitaryObligate
  24, 24, 24, 24, 24, 24, 24, 24, 24,
  
  # Pathogens
  23, 23, 23
)

# Create a DataFrame
legend_df <- data.frame(
  Xdata = rnorm(length(legend_families)),
  Ydata = rnorm(length(legend_families)),
  Family = legend_families,
  Color = legend_colors,
  Shape = legend_shapes
)

# Make a scatter plot just for the legend
gplot <- ggplot(legend_df, aes(x = Xdata, y = Ydata, fill = Family, shape = Family)) +
  geom_point(size = 7, color = "black") +  # 'color' sets outline, 'fill' sets inside
  scale_fill_manual(values = setNames(legend_colors, legend_families)) +
  scale_shape_manual(values = setNames(legend_shapes, legend_families)) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  ) +
  guides(
    fill = guide_legend(ncol = 1),  # stacked vertically
    shape = guide_legend(ncol = 1)  # stacked vertically
  ) +
  labs(fill = "Bacteria Family", shape = "Bacteria Family")

# Extract legend only
panelG <- get_legend(gplot)

# Plot just the legend
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
pdf("figures/grid_trees_all.pdf",
    height=20, width=24)  # Adjust width and height as needed

# Plot the final combined figure
grid.arrange(final_plot)

# Close the PDF device to complete saving
dev.off()
