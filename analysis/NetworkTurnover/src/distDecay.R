# Module: Distance Decay Model Analysis and Visualization
 #
 # Description:
 # This module provides functions for analyzing and visualizing decay models in ecological data. Specifically, it includes:
 # 1. **genusspecies_decay_model**: Performs a decay model analysis for genus-species data, calculating dissimilarity matrices for microbial communities 
 #    and geographic distance, and performing a Mantel test to evaluate their correlation.
 # 2. **microbe_type_decay_model**: Computes decay models for different types of microbes (Obligate or Facultative) in ecological datasets, using 
 #    Bray-Curtis dissimilarity for abundance data and geographic distance, and returns a decay model and a Mantel test.
 # 3. **plot_decay_ggplot_single**: Creates a ggplot visualization of a single decay model, including raw data points and the fitted decay curve. 
 #    It provides options for adjusting plot appearance, such as color, transparency, and line style.
 # 4. **plot_decay_ggplot_combined**: Generates a combined ggplot visualization for two decay models, allowing for the inclusion of both 
 #    raw data points and fitted curves for two datasets with customizable appearance.
 #
 # The functions in this module are designed to analyze the relationship between ecological dissimilarities (e.g., microbial community 
 # dissimilarities) and geographic distance, using both statistical testing (Mantel test) and visualizations (ggplot2). The decay models 
 # produced by these functions help to explore the distance-decay relationship, a key concept in ecology.
 #
 # Functions:
 # - genusspecies_decay_model: Performs decay model analysis for genus-species data, with a Mantel test and a decay model.
 # - microbe_type_decay_model: Computes decay models for Obligate or Facultative microbes using Bray-Curtis dissimilarity and geographic distance.
 # - plot_decay_ggplot_single: Visualizes a single decay model, including raw data points and fitted curve, with customizable appearance.
 # - plot_decay_ggplot_combined: Creates a combined plot for two decay models, with raw data points and fitted curves for each model, allowing 
 #   comparison between them.
 #
 # Notes:
 # - The decay models rely on calculating dissimilarities (using Bray-Curtis) for community data and geographic distances (using Haversine 
 #   or Euclidean methods).
 # - Mantel tests assess the correlation between ecological dissimilarity and geographic distance to evaluate the distance-decay relationship.
 # - The visualization functions use `ggplot2` to produce customizable plots, enabling users to visualize both the data points and the 
 #   fitted decay curves, with options for adjusting the plot aesthetics.
 #
 # Dependencies:
 # - `vegan`: for calculating Bray-Curtis dissimilarities (`vegdist`).
 # - `geosphere`: for calculating geographic distances using Haversine (`distm`).
 # - `betapart`: for performing decay model analysis (`decay.model`).
 # - `ggplot2`: for creating visualizations of the decay models and data.
 #
 # Example Usage:
 # 1. For genus-species decay modeling:
 #    decay_model <- genusspecies_decay_model(data, model_type = "exp")
 #
 # 2. For microbe type decay modeling:
 #    microbe_decay_model <- microbe_type_decay_model(data, type = "Obligate", model_type = "exp")
 #
 # 3. To visualize a single decay model:
 #    plot_decay_ggplot_single(decay_model, xlab = "Geographic Distance", ylab = "Dissimilarity")
 #
 # 4. To visualize combined decay models:
 #    plot_decay_ggplot_combined(decay_model1, decay_model2, xlab = "Geographic Distance", ylab = "Dissimilarity")
 
 
 set.seed(777)
 library(vegan)
 library(geosphere)
 library(betapart)
 library(tidyverse)
 
 
 # Function: genusspecies_decay_model
 #
 # Description:
 # This function performs a decay model analysis to explore the relationship between ecological dissimilarity (based on Bray-Curtis 
 # dissimilarity of 16S abundance data) and geographic distance using the Haversine method. The function filters data based on a 
 # specified taxonomic level (either "Genus" or "GenusSpecies"), calculates the Bray-Curtis dissimilarity matrix for species 
 # abundance data, and computes the geographic distance matrix. It then runs a Mantel test to assess the correlation between 
 # abundance dissimilarity and geographic distance, followed by fitting a decay model using the `betapart` package.
 #
 # Parameters:
 # - data: A data frame containing the ecological and geographic data, including columns for taxonomic identification (Genus, 
 #   GenusSpecies), 16S abundance data, and site coordinates (Longitude, Latitude).
 # - type: A string indicating the taxonomic level to analyze. Acceptable values are "Genus" and "GenusSpecies".
 # - which: A string specifying the particular genus or genus-species combination to filter for in the data.
 # - model.type: A string specifying the type of model to fit in the `betapart::decay.model` function. Acceptable values are 
 #   typically "exp".
 #
 # Output:
 # - Prints the result of the Mantel test comparing Bray-Curtis dissimilarity and geographic distance (`abund_geo`).
 # - Returns a fitted decay model object from the `betapart::decay.model` function, representing the relationship between 
 #   abundance dissimilarity and geographic distance.
 #
 # Notes:
 # - The Bray-Curtis dissimilarity matrix is calculated based on the abundance data of 16S sequences.
 # - The geographic distance is computed using the Haversine formula, and the distances are converted to kilometers.
 # - The Mantel test is used to assess the statistical correlation between ecological dissimilarity and geographic distance 
 #   using the Spearman method.
 # - The function requires the `vegdist` function from the `vegan` package to calculate the Bray-Curtis dissimilarity and the 
 #   `distm` function from the `geosphere` package for geographic distance calculation.
 # - The decay model is fit using the `decay.model` function from the `betapart` package.
 #
 # The function will output the result of the Mantel test and return the fitted decay model object.
 #
 genusspecies_decay_model <- function(data, type, which, model.type){
   #bray curtis dissimilarity matrix of 16s
 
   if(type == 'Genus'){
     abund <- data %>%
       filter(Genus == which) %>%
       select(UniqueID, starts_with('16s')) %>%
       select(-UniqueID)
 
     #distance matrix of sites
     geo <- data %>%
       filter(Genus == which) %>%
       select(UniqueID, Long, Lat) %>%
       select(-UniqueID) %>%
       mutate()
   } else if (type == 'GenusSpecies'){
 
     abund <- data %>%
       filter(GenusSpecies == which) %>%
       select(UniqueID, starts_with('16s')) %>%
       select(-UniqueID)
     #distance matrix of sites
     geo <- data %>%
       filter(GenusSpecies == which) %>%
       select(UniqueID, Long, Lat) %>%
       select(-UniqueID) %>%
       mutate()
   }
 
   #abundance data frame - bray curtis dissimilarity
   dist.abund <- vegdist(abund, method = "bray")
 
   #geographic data frame - haversine distance in m (takes a df with lat and long and calculates dist)
   d.geo <- distm(geo, fun = distHaversine)
   dist.geo <- as.dist(d.geo)/1000
 
   #abundance vs geographic mantel test
   abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 999, na.rm = TRUE)
   print(abund_geo)
 
   dist_decay_model <- betapart::decay.model(dist.abund,
                                             dist.geo,
                                             y.type='dissim',
                                             model.type = model.type,
                                             perm=999)
 
   dist_decay_model
 
 }
 
 #' Fit Decay Model for Microbe Type
 #'
 #' This function fits a decay model to test the relationship between Bray-Curtis dissimilarity
 #' in microbial community composition and geographic distance, using species filtered by microbial type.
 #'
 #' @param data A data frame containing 16S abundance data, microbial type weights, and geographic coordinates.
 #' @param type A character string specifying the microbial type to analyze. 
 #'        Options include: "ObligateAll", "TransientAll", 
 #'        "ObligateSocial", "TransientSocial", "ObligateSolitary", "TransientSolitary".
 #' @param model.type A character string indicating the model type to be used in `betapart::decay.model` (e.g., "exp", "power").
 #' @param decay.type A character string for `y.type` in the `decay.model` function (e.g., "dissimilarities" or "similarities").
 #'
 #' @return A fitted decay model object from `betapart::decay.model`. Also prints results from a Mantel test
 #'         comparing Bray-Curtis dissimilarity with geographic distance.
 #'
 #' @details
 #' This function:
 #' \itemize{
 #'   \item Filters the input data by a microbial type.
 #'   \item Calculates Bray-Curtis dissimilarity using `vegan::vegdist`.
 #'   \item Computes pairwise geographic distances using `geosphere::distm`.
 #'   \item Performs a Mantel test for correlation.
 #'   \item Fits a decay model using `betapart::decay.model`.
 #' }
 #'
 #' @importFrom vegan vegdist mantel
 #' @importFrom geosphere distm distHaversine
 #' @importFrom dplyr filter select mutate
 #' @export
 microbe_type_decay_model <- function(data, 
                                      type,
                                      model.type,
                                      decay.type) {
    
    # Select microbial abundance and coordinates based on type
    if (type == 'ObligateAll') {
       abund <- data %>%
          filter(WeightsObligateMicrobe == 1) %>%
          select(UniqueID, starts_with('16s')) %>%
          select(-UniqueID)
       
       geo <- data %>%
          filter(WeightsObligateMicrobe == 1) %>%
          select(Long, Lat)
       
    } else if (type == 'TransientAll') {
       abund <- data %>%
          filter(WeightsTransientMicrobe == 1) %>%
          select(UniqueID, starts_with('16s')) %>%
          select(-UniqueID)
       
       geo <- data %>%
          filter(WeightsTransientMicrobe == 1) %>%
          select(Long, Lat)
       
    } else if (type == 'ObligateSocial') {
       abund <- data %>%
          filter(WeightsObligateSocialMicrobe == 1) %>%
          select(UniqueID, starts_with('16s')) %>%
          select(-UniqueID)
       
       geo <- data %>%
          filter(WeightsObligateSocialMicrobe == 1) %>%
          select(Long, Lat)
       
    } else if (type == 'TransientSocial') {
       abund <- data %>%
          filter(WeightsTransientSocialMicrobe == 1) %>%
          select(UniqueID, starts_with('16s')) %>%
          select(-UniqueID)
       
       geo <- data %>%
          filter(WeightsTransientSocialMicrobe == 1) %>%
          select(Long, Lat)
       
    } else if (type == 'ObligateSolitary') {
       abund <- data %>%
          filter(WeightsObligateSolitaryMicrobe == 1) %>%
          select(UniqueID, starts_with('16s')) %>%
          select(-UniqueID)
       
       geo <- data %>%
          filter(WeightsObligateSolitaryMicrobe == 1) %>%
          select(Long, Lat)
       
    } else if (type == 'TransientSolitary') {
       abund <- data %>%
          filter(WeightsTransientSolitaryMicrobe == 1) %>%
          select(UniqueID, starts_with('16s')) %>%
          select(-UniqueID)
       
       geo <- data %>%
          filter(WeightsTransientSolitaryMicrobe == 1) %>%
          select(Long, Lat)
       
    } else {
       stop("Unknown microbe type specified.")
    }
    
    # Bray-Curtis dissimilarity matrix 
    dist.abund <- vegan::vegdist(abund, method = "bray")
    
    # Geographic distance matrix (in km)
    d.geo <- geosphere::distm(geo, fun = distHaversine)
    dist.geo <- as.dist(d.geo) / 1000  # convert meters to kilometers
    
    # Mantel test: abundance vs geographic distance
    abund_geo <- vegan::mantel(dist.abund, dist.geo, method = "spearman", permutations = 999, na.rm = TRUE)
    print(abund_geo)
    
    # Transform abundances into similarities by subtracting the dissimilarity from 1 if desired
    if (decay.type == 'similarities'){
      dist.abund <- 1 - dist.abund
    }
    # Adjust for zero distances in power models to avoid -Inf
    if (model.type == "power") {
       dist.geo[dist.geo == 0] <- 0.0001
    }
    
    # Fit the decay model
    dist_decay_model <- betapart::decay.model(
       dist.abund,
       dist.geo,
       y.type = decay.type,
       model.type = model.type,
       perm = 999
    )
    
    return(dist_decay_model)
 }
 
 
 #' Plot a Single Decay Model Using ggplot2
 #'
 #' This function creates a ggplot visualization of a single distance decay model fit. It plots 
 #' the raw data points and overlays a fitted decay curve based on model predictions, typically 
 #' from a `betapart::decay.model` object.
 #'
 #' @param x A list or model object containing:
 #'   \describe{
 #'     \item{data.x}{A numeric vector of independent variable values (e.g., geographic distances).}
 #'     \item{data.y}{A numeric vector of dependent variable values (e.g., dissimilarities).}
 #'     \item{model}{A fitted decay model (e.g., from `betapart::decay.model`).}
 #'   }
 #' @param xlab A character string for the x-axis label.
 #' @param ylab A character string for the y-axis label.
 #' @param mod1color Color for the fitted decay curve (default is `"navy"`).
 #' @param col Color for the raw data points (default is `"black"`).
 #' @param lty Line type for the fitted curve (default is `"solid"`; currently unused).
 #' @param lwd Line width for the fitted decay curve (default is `1.5`).
 #' @param cex Size of the raw data points (default is `1`).
 #'
 #' @return A `ggplot` object showing the decay relationship with raw data and model fit.
 #'
 #' @examples
 #' plot_decay_ggplot_single(my_decay_model, xlab = "Distance (km)", ylab = "Dissimilarity")
 #'
 #' @export
 plot_decay_ggplot_single <- function(x,
                                      xlab,
                                      ylab,
                                      mod1color = 'navy',
                                      col = "black",
                                      lty = "solid",
                                      lwd = 1.5,
                                      cex = 1) {
    
    # Create a dataframe of input data
    data <- data.frame(
       x = x$data.x,
       y = x$data.y
    )
    
    # Sort the data by x to ensure smooth fitted line plotting
    sorted_data <- data[order(data$x), ]
    
    # Extract model and compute fitted values ordered by x
    model <- x$model
    fitted_values <- fitted(model)[order(data$x)]
    
    # Combine sorted x with fitted values
    fit_df <- data.frame(
       x = sorted_data$x,
       y = fitted_values
    )
    
    # Build the ggplot
    p <- ggplot(data, aes(x = x, y = y)) +
       # Add raw data points with transparency and custom size
       geom_point(color = "black", fill = mod1color, alpha = 0.1, shape = 21, size = cex * 3) +
       
       # Add fitted decay curve (first as a thick black line for outline)
       geom_line(data = fit_df, aes(x = x, y = 1-y), color = "black", linewidth = lwd + 1) +
       
       # Overlay colored fitted line for emphasis
       geom_line(data = fit_df, aes(x = x, y = 1-y), color = mod1color, linewidth = lwd) +
       
       # Customize axes labels
       labs(x = xlab, y = ylab) +
       
       # Clean theme with enhanced text sizing
       theme_classic() +
       theme(
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          text = element_text(size = 16)
       )
    
    return(p)
 }
 
 
 #' Calculate community dissimilarity decay with geographic distance for bee hosts
 #'
 #' This function computes pairwise Bray-Curtis dissimilarities (or similarities) between sites
 #' based on host-associated microbe-weighted bee communities. It then fits a distance decay model
 #' using exponential decay and optionally returns a decay model for "All", "Social", or "Solitary" hosts.
 #'
 #' @param data A data frame containing bee and site information. Must include columns: `Site`, `GenusSpecies`, `Genus`, `Lat`, `Long`, and `WeightsMicrobe`.
 #' @param host.type Character string specifying host group to use. Options: `"All"`, `"Social"`, or `"Solitary"`.
 #' @param sim.type Character string specifying output type. Options: `"similarities"` or `"dissimilarities"`.
 #'
 #' @return An object of class `decay.model` from `betapart::decay.model`, representing the distance decay model fitted to the data.
 #' Also prints results of a Mantel test comparing Bray-Curtis dissimilarity to geographic distance.
 #' @export
 #'
 #' @examples
 #' get_bee_community_dissim(mydata, host.type = "Social", sim.type = "dissimilarities")
 
 get_bee_community_dissim <- function(data, host.type, sim.type) {
    
    # Filter and aggregate abundance data depending on host type
    if (host.type == "All") {
       abund <- data %>%
          filter(WeightsMicrobe == 1) %>%
          select(Site, GenusSpecies) %>%
          group_by(Site, GenusSpecies) %>%
          arrange(Site) %>%
          filter(GenusSpecies != "") %>%
          summarize(HostAbund = n()) %>%
          pivot_wider(names_from = GenusSpecies, values_from = HostAbund) %>%
          mutate_all(~replace(., is.na(.), 0))
       
       sitenames <- abund$Site
       abund$Site <- NULL
       rownames(abund) <- sitenames
       
       # Get site coordinates
       geo <- data %>%
          filter(WeightsMicrobe == 1) %>%
          group_by(Site) %>%
          select(Site, Long, Lat) %>%
          distinct() %>%
          arrange(Site)
       
       rownames(geo) <- sitenames
       geo$Site <- NULL
    }
    
    if (host.type == "Social") {
       abund <- data %>%
          filter(WeightsMicrobe == 1) %>%
          filter(Genus %in% c("Apis", "Bombus")) %>%
          select(Site, GenusSpecies) %>%
          group_by(Site, GenusSpecies) %>%
          arrange(Site) %>%
          filter(GenusSpecies != "") %>%
          summarize(HostAbund = n()) %>%
          pivot_wider(names_from = GenusSpecies, values_from = HostAbund) %>%
          mutate_all(~replace(., is.na(.), 0))
       
       sitenames <- abund$Site
       abund$Site <- NULL
       rownames(abund) <- sitenames
       
       geo <- data %>%
          filter(WeightsMicrobe == 1) %>%
          filter(Genus %in% c("Apis", "Bombus")) %>%
          group_by(Site) %>%
          select(Site, Long, Lat) %>%
          distinct() %>%
          arrange(Site)
       
       rownames(geo) <- sitenames
       geo$Site <- NULL
    }
    
    if (host.type == "Solitary") {
       abund <- data %>%
          filter(WeightsMicrobe == 1) %>%
          filter(Genus %in% c("Anthophora", "Andrena", "Melissodes", "Megachile")) %>%
          select(Site, GenusSpecies) %>%
          group_by(Site, GenusSpecies) %>%
          arrange(Site) %>%
          filter(GenusSpecies != "") %>%
          summarize(HostAbund = n()) %>%
          pivot_wider(names_from = GenusSpecies, values_from = HostAbund) %>%
          mutate_all(~replace(., is.na(.), 0))
       
       sitenames <- abund$Site
       abund$Site <- NULL
       rownames(abund) <- sitenames
       
       geo <- data %>%
          filter(WeightsMicrobe == 1) %>%
          filter(Genus %in% c("Anthophora", "Andrena", "Melissodes", "Megachile")) %>%
          group_by(Site) %>%
          select(Site, Long, Lat) %>%
          distinct() %>%
          arrange(Site)
       
       rownames(geo) <- sitenames
       geo$Site <- NULL
    }
   if (host.type == "Bombus") {
     abund <- data %>%
       filter(WeightsMicrobe == 1) %>%
       filter(Genus == "Bombus") %>%
       select(Site, GenusSpecies) %>%
       group_by(Site, GenusSpecies) %>%
       arrange(Site) %>%
       filter(GenusSpecies != "") %>%
       summarize(HostAbund = n()) %>%
       pivot_wider(names_from = GenusSpecies, values_from = HostAbund) %>%
       mutate_all(~replace(., is.na(.), 0))
     
     sitenames <- abund$Site
     abund$Site <- NULL
     rownames(abund) <- sitenames
     
     geo <- data %>%
       filter(WeightsMicrobe == 1) %>%
       filter(Genus == "Bombus") %>%
       group_by(Site) %>%
       select(Site, Long, Lat) %>%
       distinct() %>%
       arrange(Site)
     
     rownames(geo) <- sitenames
     geo$Site <- NULL
     browser()
   }
   if (host.type == "Apis") {
     abund <- data %>%
       filter(WeightsMicrobe == 1) %>%
       filter(Genus == "Apis") %>%
       select(Site, GenusSpecies) %>%
       group_by(Site, GenusSpecies) %>%
       arrange(Site) %>%
       filter(GenusSpecies != "") %>%
       summarize(HostAbund = n()) %>%
       pivot_wider(names_from = GenusSpecies, values_from = HostAbund) %>%
       mutate_all(~replace(., is.na(.), 0))
     
     sitenames <- abund$Site
     abund$Site <- NULL
     rownames(abund) <- sitenames
     
     geo <- data %>%
       filter(WeightsMicrobe == 1) %>%
       filter(Genus == "Apis") %>%
       group_by(Site) %>%
       select(Site, Long, Lat) %>%
       distinct() %>%
       arrange(Site)
     
     rownames(geo) <- sitenames
     geo$Site <- NULL
   }
    
    # Calculate Bray-Curtis dissimilarity matrix from abundance data
    dist.abund <- vegdist(abund, method = "bray")
    
    # Convert to similarity matrix if requested
    if (sim.type == "similarities") {
       dist.sim <- 1 - dist.abund
    } else {
       dist.sim <- dist.abund
    }
    
    # Compute pairwise geographic distances (in kilometers)
    d.geo <- distm(geo, fun = distHaversine)
    dist.geo <- as.dist(d.geo) / 1000  # convert meters to km
    
    # Mantel test between Bray-Curtis dissimilarity and geographic distance
    set.seed(1298)
    abund_geo <- mantel(dist.abund, dist.geo, method = "spearman", permutations = 999, na.rm = TRUE)
    print(abund_geo)
    
    # Fit decay model (from betapart)
    set.seed(1298)
    dist_decay_model <- betapart::decay.model(
       dist.sim,
       dist.geo,
       y.type = sim.type,
       model.type = "exponential",
       perm = 999
    )
    
    return(dist_decay_model)
 }
 
 
 # Function: plot_decay_ggplot_combined
 #
 # Description:
 # This function creates a combined ggplot visualization for two decay models. The plot displays both the raw data points and the fitted 
 # decay curves for two datasets, each corresponding to a different model. The points are jittered to prevent overplotting, and the fitted 
 # curves represent the decay relationship for each model. This function can either show the data points or just the fitted curves, 
 # depending on the `add.points` parameter.
 #
 # Parameters:
 # - x: A list or model object containing the first decay model output and its data. This object must include:
 #   - x$data.x: A numeric vector or column representing the independent variable (e.g., geographic distance or other explanatory variable).
 #   - x$data.y: A numeric vector or column representing the dependent variable (e.g., Bray-Curtis dissimilarity or ecological distance).
 #   - x$model: A fitted model object (e.g., decay model) containing the model's predictions.
 # - z: A list or model object containing the second decay model output and its data. The structure is similar to `x`.
 # - xlab: A string representing the label for the x-axis (e.g., "Distance").
 # - ylab: A string representing the label for the y-axis (e.g., "Dissimilarity").
 # - mod1color: A string specifying the color of the first model's fitted curve (default is 'navy').
 # - mod2color: A string specifying the color of the second model's fitted curve (default is 'gold').
 # - alpha1: A numeric value controlling the transparency of the first dataset points (default is 0.1).
 # - alpha2: A numeric value controlling the transparency of the second dataset points (default is 0.5).
 # - col: A string specifying the color of the points (default is 'black'). This color is used for the points in both datasets.
 # - lty1: A string specifying the line type for the first fitted curve (e.g., "solid", "dashed").
 # - lty2: A string specifying the line type for the second fitted curve (e.g., "solid", "dashed").
 # - lwd: A numeric value specifying the line width for the fitted curves (default is 1.5).
 # - cex: A numeric value for the size of the scatter points (default is 1).
 # - add.points: A logical value (TRUE or FALSE) that determines whether to include the raw data points in the plot. If FALSE, only the fitted curves are shown (default is TRUE).
 #
 # Output:
 # - Returns a ggplot object representing the combined decay plot. The plot consists of:
 #   - Jittered points for both datasets (if `add.points` is TRUE).
 #   - Fitted decay curves for both models, each drawn with a unique color and line type.
 #   - X and Y axis labels as defined by the `xlab` and `ylab` parameters.
 #
 # Notes:
 # - The function assumes that `x` and `z` contain fitted decay models and their corresponding data.
 # - The fitted model's predictions are used to create the fitted lines, with the data ordered by the independent variable for both models.
 # - The jittered points use `geom_point` with random noise (using `position_jitter`) to avoid overplotting.
 # - If `add.points` is FALSE, only the fitted curves are shown in the plot.
 # - The plot is created using `ggplot2`, with the appearance customized by `theme_classic` and adjustments to axis labels, text size, and line appearance.
 #
 # plot_decay_ggplot_combined <- function(x,
 #                                        z,
 #                                        xlab,
 #                                        ylab,
 #                                        mod1color='navy',
 #                                        mod2color='gold',
 #                                        alpha1 = 0.1,
 #                                        alpha2 = 0.5,
 #                                        col = "black",
 #                                        lty1,
 #                                        lty2,
 #                                        lwd = 1.5,
 #                                        cex = 1,
 #                                        add.points=TRUE,
 #                                        log.dist=FALSE) {
 # 
 #   # Extract data and fitted values
 #   data1 <- data.frame(x$data.x, x$data.y)
 #   model1 <- x$model
 #   fitted_values1 <- data.frame(fitted(model1)[order(data1$x.data.x)])
 #   sorted_data1 <- data1[order(data1$x.data.x), ]
 # 
 #   # Extract data and fitted values
 #   data2 <- data.frame(z$data.x, z$data.y)
 #   model2 <- z$model
 #   fitted_values2 <- data.frame(fitted(model2)[order(data2$z.data.x)])
 #   sorted_data2 <- data2[order(data2$z.data.x), ]
 # 
 #   if(add.points==TRUE){
 #   # Create the ggplot object
 #   p <- ggplot(data1, aes(x = x.data.x, y = x.data.y)) +
 #     geom_point(data=data1, aes(x = x.data.x, y = x.data.y), fill = mod1color, alpha = alpha1 , color="black", pch=21, cex=3, position = position_jitter(w = 10, h = 0)) +
 #     geom_point(data=data2, aes(x = z.data.x, y = z.data.y), fill = mod2color, alpha = alpha2 , color="black", pch=21, cex=3, position = position_jitter(w = 10, h = 0)) +
 #     geom_line(aes(x = sorted_data1$x.data.x, y = fitted_values1$fitted.model1..order.data1.x.data.x..),
 #               color = 'black', linewidth=2.5,) +
 #     geom_line(aes(x = sorted_data1$x.data.x, y = fitted_values1$fitted.model1..order.data1.x.data.x..), color = mod1color, linewidth=2, linetype=lty1) +
 #     geom_line(data=data2, aes(x = sorted_data2$z.data.x, y = fitted_values2$fitted.model2..order.data2.z.data.x..),
 #               color = 'black', linewidth=2.5,) +
 #     geom_line(data=data2, aes(x = sorted_data2$z.data.x, y = fitted_values2$fitted.model2..order.data2.z.data.x..), color = mod2color, linewidth=2, linetype=lty2) +
 #     labs(x=xlab, y=ylab) +
 #     theme_classic() +
 #     theme(axis.title.x = element_text(size=16),
 #           axis.title.y = element_text(size=16),
 #           text = element_text(size=16))
 #   }
 #   
 #   if(log.dist==TRUE){
 #      # Create the ggplot object
 #      p <- ggplot(data1, aes(x = x.data.x, y = x.data.y)) +
 #         geom_point(data=data1, aes(x = x.data.x, y = x.data.y), fill = mod1color, alpha = alpha1 , color="black", pch=21, cex=3, position = position_jitter(w = 0.1, h = 0)) +
 #         geom_point(data=data2, aes(x = z.data.x, y = z.data.y), fill = mod2color, alpha = alpha2 , color="black", pch=21, cex=3, position = position_jitter(w = 0.1, h = 0)) +
 #         geom_line(aes(x = sorted_data1$x.data.x, y = fitted_values1$fitted.model1..order.data1.x.data.x..),
 #                   color = 'black', linewidth=2.5,) +
 #         geom_line(aes(x = sorted_data1$x.data.x, y = fitted_values1$fitted.model1..order.data1.x.data.x..), color = mod1color, linewidth=2, linetype=lty1) +
 #         geom_line(data=data2, aes(x = sorted_data2$z.data.x, y = fitted_values2$fitted.model2..order.data2.z.data.x..),
 #                   color = 'black', linewidth=2.5,) +
 #         geom_line(data=data2, aes(x = sorted_data2$z.data.x, y = fitted_values2$fitted.model2..order.data2.z.data.x..), color = mod2color, linewidth=2, linetype=lty2) +
 #         labs(x=xlab, y=ylab) +
 #         scale_x_log10()+
 #            #limits = c(60, 350)) +  # Adjust upper limit based on your data
 #         #   breaks = c(1, 10, 100, 1000),
 #         #   labels = scales::label_number()) +
 #         theme_classic() +
 #         theme(axis.title.x = element_text(size=16),
 #               axis.title.y = element_text(size=16),
 #               text = element_text(size=16)) 
 #   }
 # 
 #   if(add.points==FALSE){
 #     # Create the ggplot object
 #     p <- ggplot(data1, aes(x = x.data.x, y = x.data.y)) + 
 #       geom_line(aes(x = sorted_data1$x.data.x, y = (1 - fitted_values1$fitted.model1..order.data1.x.data.x..)),
 #                 color = 'black', linewidth=2.5,) +
 #       geom_line(aes(x = sorted_data1$x.data.x, y = (1 - fitted_values1$fitted.model1..order.data1.x.data.x..)), color = mod1color, linewidth=2, linetype=lty1) +
 #       geom_line(data=data2, aes(x = sorted_data2$z.data.x, y = (1 - fitted_values2$fitted.model2..order.data2.z.data.x..)),
 #                 color = 'black', linewidth=2.5,) +
 #       geom_line(data=data2, aes(x = sorted_data2$z.data.x, y = (1 - fitted_values2$fitted.model2..order.data2.z.data.x..)), color = mod2color, linewidth=2, linetype=lty2) +
 #       labs(x=xlab, y=ylab) +
 #       theme_classic() +
 #       theme(axis.title.x = element_text(size=16),
 #             axis.title.y = element_text(size=16),
 #             text = element_text(size=16)) + ylim(0,1)
 #   }
 # 
 #   return(p)
 # 
 # }
 

 
 plot_decay_ggplot_combined <- function(x,
                                        z,
                                        xlab,
                                        ylab,
                                        mod1color='navy',
                                        group1label='Social Hosts',
                                        mod2color='gold',
                                        group2label='Solitary Hosts',
                                        alpha1 = 0.1,
                                        alpha2 = 0.5,
                                        col = "black",
                                        lty1,
                                        lty2,
                                        lwd = 1.5,
                                        cex = 1,
                                        add.points=TRUE,
                                        log.dist=FALSE,
                                        add.legend=FALSE,
                                        legend.title=NULL) {
    # Prepare data with group for legend
    data1 <- data.frame(x = x$data.x, y = x$data.y, group = group1label)
    data2 <- data.frame(x = z$data.x, y = z$data.y, group = group2label)
    model1 <- x$model
    model2 <- z$model
    
    sorted_data1 <- data1[order(data1$x), ]
    sorted_data2 <- data2[order(data2$x), ]
    sorted_data1$fitted <- fitted(model1)[order(data1$x)]
    sorted_data2$fitted <- fitted(model2)[order(data2$x)]
    
    p <- ggplot()
    
    if(add.points){
       p <- p +
          geom_point(data = data1, aes(x = x, y = y, fill = group),
                     alpha = alpha1, color = "black", shape = 21, size = 3,
                     position = position_jitter(w = ifelse(log.dist, 0.1, 10), h = 0)) +
          geom_point(data = data2, aes(x = x, y = y, fill = group),
                     alpha = alpha2, color = "black", shape = 21, size = 3,
                     position = position_jitter(w = ifelse(log.dist, 0.1, 10), h = 0))
    }
    
    # Add lines regardless of add.points
    p <- p +
       geom_line(data = sorted_data1, aes(x = x, y = fitted, linetype = group, color = group), linewidth = 2) +
       geom_line(data = sorted_data2, aes(x = x, y = fitted, linetype = group, color = group), linewidth = 2)
    
    # Add scales (always include to support group mapping)
    p <- p +
       scale_fill_manual(name = legend.title, values = setNames(c(mod1color, mod2color), c(group1label, group2label))) +
       scale_color_manual(name = legend.title, values = setNames(c(mod1color, mod2color), c(group1label, group2label))) +
       scale_linetype_manual(name = legend.title, values = setNames(c(lty1, lty2), c(group1label, group2label)))
    
    # Axes, theme, etc.
    p <- p + labs(x = xlab, y = ylab) +
       theme_classic(base_size = 16)
    
    if(log.dist){
       p <- p + scale_x_log10()
    }
    
    if(!add.legend){
       p <- p + theme(legend.position = "none")
    } else {
       p <- p + theme(
          legend.position = c(0.95, 0.95),
          legend.justification = c("right", "top"),
          legend.key.size = unit(0.5, "lines"),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 11)
       )
    }
    
    return(p)
 }
 
 