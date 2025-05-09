## **********************************************************
## Load libraries
## **********************************************************

rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path
setwd(local.path)
setwd("microbeBiogeography/analysis/SEM/")

## Toggle on or off to run individual genus models
run.bombus = TRUE
run.melissodes = TRUE

library(picante)
library(plyr)
library(bayesplot)
library(pscl)
library(brms)
library(performance)
library(R2admb)
library(shinystan)
## If not already installed
# install.packages("glmmADMB", 
#                  repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#                          getOption("repos")),
#                  type="source")


## **********************************************************
## Standardize, center, and transform data
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

## **********************************************************
## Load Data
## **********************************************************


## TODO need to fix colnames here and leading X / double .. in microbe colnames
load("../../microbeBiogeographyData.Rdata")

ncores <- 1


## **********************************************************
## Flower diversity
## **********************************************************

flower.div.vars <- c("Lat",
                     "(1|Site)")

flower.div.x <- paste(flower.div.vars, collapse="+")
flower.div.y <- "MeanFloralDiversity | subset(Weights)"
formula.flower.div <- as.formula(paste(flower.div.y, "~",flower.div.x))

## **********************************************************
## Bee abundance
## **********************************************************

tot.bee.abund.vars <- c("(1|Site)")

tot.bee.abund.x <- paste(tot.bee.abund.vars, collapse="+")
tot.bee.abund.y <- "BeeAbundance | subset(Weights)"
formula.tot.bee.abund <- as.formula(paste(tot.bee.abund.y, "~",tot.bee.abund.x))

## **********************************************************
## Bee diversity
## **********************************************************

tot.bee.div.vars <- c("MeanFloralDiversity",
                      "Lat",
                      "(1|Site)")

tot.bee.div.x <- paste(tot.bee.div.vars, collapse="+")
tot.bee.div.y <- "BeeDiversity | subset(Weights)"
formula.tot.bee.div <- as.formula(paste(tot.bee.div.y, "~",tot.bee.div.x))

## **********************************************************
## convert formulas to brms format
## **********************************************************
bf.fdiv <- bf(formula.flower.div)
bf.tot.babund <- bf(formula.tot.bee.abund)
bf.tot.bdiv <- bf(formula.tot.bee.div)

## **********************************************************
## change NAs to 0 to prevent brms dropping missing data
## **********************************************************

spec.net[is.na(spec.net)] <- 0

## **********************************************************
## Bombus microbe PD model
## **********************************************************

## Obligate PD model
ob.microbe.bombus.vars <- c("BeeAbundance",
                         "BeeDiversity",
                         "Lat", 
                         "MeanFloralDiversity",
                         "MeanITD",
                         "rare.degree",
                         "(1|Site)",
                         "(1|gr(GenusSpecies, cov = phylo_matrix))") 


ob.microbe.bombus.x <- paste(ob.microbe.bombus.vars, collapse="+")
ob.microbe.bombus.y <- "PD.obligate.log | subset(WeightsObligateBombus)"
formula.ob.microbe.bombus <- as.formula(paste(ob.microbe.bombus.y, "~",
                                           ob.microbe.bombus.x))

bf.ob.microbe.bombus.skew <- bf(formula.ob.microbe.bombus, family=skew_normal())

## Facultative model
non.ob.microbe.bombus.vars <- c("BeeAbundance",
                            "BeeDiversity",
                            "Lat", 
                            "MeanFloralDiversity",
                            "MeanITD",
                            "rare.degree",
                            "(1|Site)",
                            "(1|gr(GenusSpecies, cov = phylo_matrix))") 

non.ob.microbe.bombus.x <- paste(non.ob.microbe.bombus.vars, collapse="+")
non.ob.microbe.bombus.y <- "PD.transient.log | subset(WeightsTransientBombus)"
formula.non.ob.microbe.bombus <- as.formula(paste(non.ob.microbe.bombus.y, "~",
                                              non.ob.microbe.bombus.x))

bf.non.ob.microbe.bombus.student <- bf(formula.non.ob.microbe.bombus, family=student())

## combined model
bform.bombus <- bf.fdiv +
  bf.tot.bdiv +
  bf.tot.babund +
  bf.ob.microbe.bombus.skew +
  bf.non.ob.microbe.bombus.student +
  set_rescor(FALSE)


if(run.bombus){
  fit.microbe.bombus <- brm(bform.bombus, spec.net,
                          cores=ncores,
                          iter = 10000,
                          chains =1,
                          thin=1,
                          init=0,
                          open_progress = FALSE,
                          control = list(adapt_delta = 0.99,
                                         max_treedepth = 15),
                          save_pars = save_pars(all = TRUE),
                          data2 = list(phylo_matrix=phylo_matrix))
  
  write.ms.table(fit.microbe.bombus, "bombus_microbe_2")
  r2loo.bombus <- loo_R2(fit.microbe.bombus)
  r2.bombus <- rstantools::bayes_R2(fit.microbe.bombus)
  save(fit.microbe.bombus, spec.net, r2.bombus, r2loo.bombus,
       file="saved/fullMicrobeBombusFit_2.Rdata")
}


## **********************************************************
## Melissodes microbe PD model
## **********************************************************
## obligate PD model
ob.microbe.melissodes.vars <- c("BeeAbundance",
                          "BeeDiversity",
                          "Lat", 
                          "MeanFloralDiversity",
                          "rare.degree",
                          "(1|Site)") 


ob.microbe.melissodes.x <- paste(ob.microbe.melissodes.vars, collapse="+")
ob.microbe.melissodes.y <- "PD.obligate.log | subset(WeightsObligateMelissodes)"
formula.ob.microbe.melissodes <- as.formula(paste(ob.microbe.melissodes.y, "~",
                                            ob.microbe.melissodes.x))

## model formula
bf.ob.microbe.melissodes.skew <- bf(formula.ob.microbe.melissodes, family=skew_normal())


## facultative PD model
non.ob.microbe.melissodes.vars <- c("BeeAbundance",
                              "BeeDiversity",
                              "Lat", 
                              "MeanFloralDiversity",
                              "rare.degree",
                              "(1|Site)") 


non.ob.microbe.melissodes.x <- paste(non.ob.microbe.melissodes.vars, collapse="+")
non.ob.microbe.melissodes.y <- "PD.transient.log | subset(WeightsTransientMelissodes)"
formula.non.ob.microbe.melissodes <- as.formula(paste(non.ob.microbe.melissodes.y, "~",
                                                non.ob.microbe.melissodes.x))


## model formula
bf.non.ob.microbe.melissodes.skew <- bf(formula.non.ob.microbe.melissodes, family=skew_normal())


##combine forms
bform.melissodes <- bf.fdiv +
  bf.tot.bdiv +
  bf.tot.babund +
  bf.ob.microbe.melissodes.skew +
  bf.non.ob.microbe.melissodes.skew +
  set_rescor(FALSE)

if(run.melissodes){
fit.microbe.melissodes <- brm(bform.melissodes , spec.net,
                           cores=ncores,
                           iter = 10000,
                           chains =1,
                           thin=1,
                           init=0,
                           open_progress = FALSE,
                           control = list(adapt_delta = 0.9999,
                                          max_treedepth = 15
                                          ),
                           save_pars = save_pars(all = TRUE))

write.ms.table(fit.microbe.melissodes, "melissodes_microbe_2")
r2loo.melissodes <- loo_R2(fit.microbe.melissodes)
r2.melissodes <- rstantools::bayes_R2(fit.microbe.melissodes)
save(fit.microbe.melissodes, spec.net, r2.melissodes, r2loo.melissodes,
      file="saved/fullMicrobeMelissodesFit_2.Rdata")
}
