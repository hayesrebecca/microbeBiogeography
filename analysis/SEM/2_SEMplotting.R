
## **********************************************************
## Load libraries
## **********************************************************

rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd("microbeBiogeography/analysis/SEM/")

library(plyr)
library(lme4)
library(R2admb)
#library(shinystan)
library(picante)
library(bayesplot)
library(pscl)
library(brms)
library(performance)
# install.packages("glmmADMB", 
#                  repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                          getOption("repos")),
#                  type="source")

## **********************************************************
## Prep model variables and load data files
## **********************************************************

## all of the variables that are explanatory variables and thus need
## to be centered

variables.to.log <- c("rare.degree",
                      "BeeAbundance"
)

vars_yearsr <- c("MeanFloralDiversity",
                 "BeeAbundance",
                 "BeeDiversity"
)

vars_yearsrsp <- "rare.degree"
vars_sp <- "MeanITD"
vars_site <- "Lat"

source("src/misc_microbe.R")

load("../../microbeBiogeographyData.Rdata")

## **********************************************************
## Set up model weights
## **********************************************************

# site level weights
spec.uni <- spec.orig[spec.orig$Weights ==1,]

## bee abund
## TODO unlog transform the data for axes?
labs.bee.abund <- (pretty(c(spec.uni$BeeAbundance), n=8))
axis.bee.abund <-  standardize.axis(labs.bee.abund,
                                  spec.uni$BeeAbundance)

## bee diversity
labs.bee.div <- (pretty(c(spec.uni$BeeDiversity), n=8))
axis.bee.div <-  standardize.axis(labs.bee.div,
                                  spec.uni$BeeDiversity)

## species level weights
spec.sp <-  spec.orig[spec.orig$WeightsObligateMicrobe == 1 & spec.orig$Genus == "Bombus",]

## rare.degree
labs.degree <- (pretty(spec.sp$rare.degree, n=10))
axis.degree <-  standardize.axis(labs.degree,
                                 spec.sp$rare.degree)

spec.sp <-  spec.orig[spec.orig$WeightsTransientMicrobe == 1 & spec.orig$Genus == "Bombus",]
## meanITD
labs.itd <- (pretty(spec.sp$MeanITD, n=8))
axis.itd <-  standardize.axis(labs.itd,
                              spec.sp$MeanITD)


## **********************************************************
## Load models
## **********************************************************

## load model output data

load(file="../../../skyIslands/analysis/microbiome/saved/fullMicrobeBombusFit_2_log.Rdata") #TODO update with filepath for this fork once working
load(file="../../../skyIslands/analysis/microbiome/saved/fullMicrobeMelissodesFit_2.Rdata") #TODO update with filepath for this fork once working


#load(file="saved/fullMicrobeBombusFit_2_log.Rdata") #TODO update with filepath for this fork once working
#load(file="saved/fullMicrobeMelissodesFit_2.Rdata") #TODO update with filepath for this fork once working


## **********************************************************
## Filter data to genus level and correct colnames for
##  plotting functions
## **********************************************************

spec.bombus <- spec.net[spec.net$Genus == "Bombus",]
bombus.obligate <- spec.bombus[spec.bombus$WeightsObligateMicrobe == 1,]
bombus.obligate$PDobligate <- bombus.obligate$PD.obligate
bombus.obligate$PDobligatelog <- bombus.obligate$PD.obligate.log
bombus.transient <- spec.bombus[spec.bombus$WeightsTransientMicrobe == 1,]
bombus.transient$PDtransient <- bombus.transient$PD.transient
bombus.transient$PDtransientlog <- bombus.transient$PD.transient.log

spec.melissodes <- spec.net[spec.net$Genus == "Melissodes",]
melissodes.obligate <- spec.melissodes[spec.melissodes$WeightsObligateMicrobe == 1,]
melissodes.obligate$PDobligate <- melissodes.obligate$PD.obligate
melissodes.obligate$PDobligatelog <- melissodes.obligate$PD.obligate.log
melissodes.transient <- spec.melissodes[spec.melissodes$WeightsTransientMicrobe == 1,]
melissodes.transient$PDtransient <- melissodes.transient$PD.transient
melissodes.transient$PDtransientlog <- melissodes.transient$PD.transient.log

## **********************************************************
## Plot conditional effects model results
## **********************************************************

source("src/plotting_functions.R")

## obligate microbe PD ~ Bee Diversity
obligate.bee.div.plot <- plot_model_condeff_compare(model.a=fit.microbe.bombus,
                                       model.b=fit.microbe.melissodes,
                                       this.effect='BeeDiversity',
                                       this.resp.a='PDobligatelog', 
                                       this.resp.b='PDobligatelog',
                                       point.data.a=bombus.obligate,
                                       point.data.b=melissodes.obligate,
                                       axis.breaks=axis.bee.div,
                                       axis.labs=labs.bee.div,
                                       xlabel="Bee Diversity",
                                       ylabel="Strong Associate PD (logged)",
                                       mod1color='navy',
                                       mod2color='gold',
                                       fill.a=TRUE,
                                       fill.b=FALSE
                                       )
obligate.bee.div.plot

panelA <- obligate.bee.div.plot + labs(tag="A.")

## Obligate microbe PD ~ Rarefied degree
  
obligate.rare.degree.plot <- plot_model_condeff_compare(model.a=fit.microbe.bombus,
                                         model.b=fit.microbe.melissodes,
                                         this.effect='rare.degree',
                                         this.resp.a='PDobligatelog', ## TODO: potentially update to both be PD obligate log?
                                         this.resp.b='PDobligatelog',
                                         point.data.a=bombus.obligate,
                                         point.data.b=melissodes.obligate,
                                         axis.breaks=axis.degree,
                                         axis.labs=labs.degree,
                                         xlabel="Diet Breadth (logged)",
                                         ylabel="Strong Associate PD (logged)",
                                         mod1color='navy',
                                         mod2color='gold',
                                         fill.a=TRUE,
                                         fill.b=FALSE
                                         )
obligate.rare.degree.plot 
panelB <- obligate.rare.degree.plot + labs(tag="B.")

## transient microbe pd ~ mean itd
## can only make for Bombus since melissodes doesn't have mean ITD

transient.itd.plot  <- plot_model_condeff_single(model=fit.microbe.bombus,
                                      this.effect='MeanITD.x',
                                      this.resp='PDtransientlog', ## TODO: potentially update to both be PD obligate log?
                                      point.data=bombus.transient,
                                      axis.breaks=axis.itd,
                                      axis.labs=labs.itd,
                                      xlabel="Body Size (mm)",
                                      ylabel="Weak Associate PD (logged)",
                                      mod1color='navy')

transient.itd.plot
panelC <- transient.itd.plot + labs(tag="C.")


## obligate microbe PD ~ Bee abundance
obligate.abund.plot <- plot_model_condeff_compare(model.a=fit.microbe.bombus,
                                                model.b=fit.microbe.melissodes,
                                                this.effect='BeeAbundance',
                                                this.resp.a="PDobligatelog",
                                                this.resp.b="PDobligatelog",
                                                point.data.a=bombus.obligate,
                                                point.data.b=melissodes.obligate,
                                                axis.breaks=axis.bee.abund,
                                                axis.labs=labs.bee.abund,
                                                xlabel="Bee Abundance (logged)",
                                                ylabel="Strong Associate PD (logged)",
                                                mod1color='navy',
                                                mod2color='gold',
                                                fill.a=FALSE,
                                                fill.b=TRUE
)
obligate.abund.plot
panelD <- obligate.abund.plot + labs(tag="D.")

## **********************************************************
## Save out panel figure
## **********************************************************

fig.dir <- "figures"
if(!dir.exists(fig.dir)) {
  dir.create(fig.dir, showWarnings = FALSE)
}

## Make panel figs and save out
pdf("figures/SEM_panel_figure.pdf", width = 8.5, height = 8.5) # Open a new pdf file
grid.arrange(panelA,
             panelB,
             panelC,
             panelD,
             ncol=2) 
dev.off()

