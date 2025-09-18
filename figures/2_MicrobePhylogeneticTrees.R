
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

load("C:/Users/rah10/Dropbox (University of Oregon)/skyIslands/data/spec_RBCL_16s.Rdata")
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
panelA <- bombus_tree[[1]] + labs(tag="A") + theme(plot.tag = element_text(size = 30))
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
panelD <- melissodes_tree[[1]] + labs(tag="D") + theme(plot.tag = element_text(size = 30))
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
panelB <- apis_tree[[1]] + labs(tag="B") + theme(plot.tag = element_text(size = 30))
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
panelF <- megachile_tree[[1]] + labs(tag="F") + theme(plot.tag = element_text(size = 30))
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
panelC <- anthophora_tree[[1]] + labs(tag="C") + theme(plot.tag = element_text(size = 30))
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
panelE <- andrena_tree[[1]] + labs(tag="E") + theme(plot.tag = element_text(size = 30))
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
library(cowplot)

# Define all families, colors, shapes
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

# Named vector mapping family names to hex colors
taxa_colors <- c(
  Acetobacteraceae = '#e6194B',
  Lactobacillaceae = '#3cb44b',
  
  Bartonella = '#ffe119',
  Bifidobacteriaceae = '#4363d8',
  Neisseriaceae = '#f58231',
  Orbaceae = '#911eb4',
  
  Bacillaceae = '#42d4f4',
  Burkholderiaceae = '#f032e6',
  Clostridiaceae = '#bfef45',
  Comamonadaceae = '#fabed4',
  Enterobacteriaceae = '#469990',
  Methylobacteriaceae = '#dcbeff',
  Moraxellaceae = '#9A6324',
  Oxalobacteraceae = '#800000',
  Sphingomonadaceae = '#aaffc3',
  
  Wolbachia = '#808000',
  Erwinia = '#ffd8b1',
  Hafnia = '#000075'
)


taxa_shapes <- c(
  Acetobacteraceae = 21, Lactobacillaceae = 21,              # circles
  Bartonella = 22, Bifidobacteriaceae = 22, Neisseriaceae = 22, Orbaceae = 22, # squares
  Bacillaceae = 24, Burkholderiaceae = 24, Clostridiaceae = 24, Comamonadaceae = 24,
  Enterobacteriaceae = 24, Methylobacteriaceae = 24, Moraxellaceae = 24,
  Oxalobacteraceae = 24, Sphingomonadaceae = 24,             # triangles
  Wolbachia = 23, Erwinia = 23, Hafnia = 23                  # diamonds
)

legend_families <- names(taxa_colors)

legend_df <- data.frame(
  Family = legend_families,
  Color = unname(taxa_colors[legend_families]),
  Shape = unname(taxa_shapes[legend_families]),
  stringsAsFactors = FALSE
)

# Sort by shape and then alphabetically
legend_df <- legend_df[order(legend_df$Shape, legend_df$Family), ]

# Factorize Family based on sorted order
legend_df$Family <- factor(legend_df$Family, levels = legend_df$Family)

# Add dummy X and Y values
legend_df$Xdata <- rnorm(nrow(legend_df))
legend_df$Ydata <- rnorm(nrow(legend_df))


library(ggplot2)
library(cowplot)

gplot <- ggplot(legend_df, aes(x = Xdata, y = Ydata, fill = Family, shape = Family)) +
  geom_point(size = 7, color = "black") +
  scale_fill_manual(values = setNames(legend_df$Color, legend_df$Family)) +
  scale_shape_manual(values = setNames(legend_df$Shape, legend_df$Family)) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.key = element_blank()
  ) +
  labs(fill = "Bacteria\n Family", shape = "Bacteria\n Family") +
  guides(
    fill = guide_legend(nrow = 3, ncol = 6, byrow = FALSE),
    shape = guide_legend(nrow = 3, ncol = 6, byrow = FALSE)
  )

panelG <- get_legend(gplot)

plot(panelG)


df_continuous <- data.frame(
  x = 1:100,
  y = 1,
  percent_sites = 0:99  # full range of values
)


# Build the plot with dummy point instead of tile
gradient_legend_plot <- ggplot(df_continuous, aes(x = x, y = y, fill = percent_sites)) +
  geom_point(shape = 21, size = 5) +
  scale_fill_gradient(
    low = "lightgrey", high = "black",
    name = "%\n Sites\n Present",
    limits = c(0, 100)
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.key.width = unit(1, "cm")
    #legend.key.height = unit(0.6, "cm")
  )

# Extract just the legend
gradient_legend <- get_legend(gradient_legend_plot)

# Plot just the legend (for visual test)
plot(gradient_legend)

## **********************************************************
## Create full paneled figure with legend
## **********************************************************

#Set up the layout matrix so that the bottom legend spans both columns
layout <- rbind(c(1, 2, 3, 4, 5, 6),
                c(7, 7, 7, 7, 7, 8)) # The legend will span both columns

# Create the final layout
final_plot <- arrangeGrob(panelA, panelB, panelC, panelD, panelE, panelF, panelG, gradient_legend,
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
