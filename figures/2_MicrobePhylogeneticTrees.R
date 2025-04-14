
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
bombus_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta=meta, "Bombus", genus.or.spp='Genus', finalASV, bombus_sites, do_collapse = TRUE, add_tip_labs = TRUE)
panelA <- bombus_tree[[1]] + labs(tag="A. Bombus (n=444)")
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
melissodes_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Melissodes", genus.or.spp='Genus', finalASV, melissodes_sites, do_collapse = TRUE, add_tip_labs = TRUE)
panelD <- melissodes_tree[[1]] + labs(tag="D. Melissodes (n=51)")
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
apis_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Apis", genus.or.spp='Genus', finalASV, apis_sites, do_collapse = TRUE, add_tip_labs = TRUE)
panelB <- apis_tree[[1]] + labs(tag="B. Apis (n=245)")
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
megachile_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Megachile", genus.or.spp='Genus', finalASV, megachile_sites, do_collapse = TRUE, add_tip_labs = TRUE)
panelF <- megachile_tree[[1]] + labs(tag="F. Megachile (n=43)")
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
anthophora_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Anthophora", genus.or.spp='Genus', finalASV, anthophora_sites, do_collapse = TRUE, add_tip_labs = TRUE)
panelC <- anthophora_tree[[1]] + labs(tag="C. Anthophora (n=38)")
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
andrena_tree <- phylotree_heatmap_byGenus(physeq16sR0, meta, "Andrena", genus.or.spp='Genus', finalASV, andrena_sites, do_collapse = TRUE, add_tip_labs=TRUE)
panelE <- andrena_tree[[1]] + labs(tag="E. Andrena (n=115)")
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

## Obligate symbionts
these_obligates <- c("Acetobacteraceae",
                     "Bartonellaceae",
                     "Bifidobacteriaceae",
                     "Lactobacillaceae",
                     "Neisseriaceae",
                     "Orbaceae")


leg_col_order <- c("#F6A600", #Acetobacteriaceae
                   "#604E97", #Bartonellaceae
                   "#B3446C", #Bifidobacteriaceae
                   "#882D17", #Lactobacillaceae
                   "#DCD300", #Neisseriaceae
                   "#8DB600") #Orbaceae

# Create a DataFrame to store legend data
data.leg <- data.frame( 
  Xdata = rnorm(6),
  Ydata = rnorm(6), 
  Family = these_obligates,
  leg_color = leg_col_order)

# Create a Scatter Plot to generate proper legend
gplot <- ggplot(data.leg, aes(Xdata, Ydata, color = Family)) +    
  geom_point(size = 7) +
  scale_color_manual(values=data.leg$leg_color) +
  theme(legend.position='bottom') +
  labs(color='Bacteria Family') +
  guides(colour = guide_legend(nrow = 1)) + theme(legend.key=element_blank(),
                                                   legend.text=element_text(size=12),
                                                  legend.title=element_text(size=12))

## Draw only legend without plot 
panelG  <- get_legend(gplot)                     
plot(get_legend(gplot) )

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
