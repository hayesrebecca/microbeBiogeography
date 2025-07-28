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
 
 # Function: microbe_type_decay_model
 #
 # Description:
 # This function performs a decay model analysis to explore the relationship between ecological dissimilarity (based on Bray-Curtis 
 # dissimilarity of 16S abundance data) and geographic distance using the Haversine method. The function filters data based on a 
 # specified microbe type (either "Obligate" or "Facultative"), calculates the Bray-Curtis dissimilarity matrix for species 
 # abundance data, and computes the geographic distance matrix. It then runs a Mantel test to assess the correlation between 
 # abundance dissimilarity and geographic distance, followed by fitting a decay model using the `betapart` package.
 #
 # Parameters:
 # - data: A data frame containing the ecological and geographic data, including columns for microbe weights (for obligate or 
 #   facultative microbes), 16S abundance data, and site coordinates (Longitude, Latitude).
 # - type: A string indicating the type of microbe to analyze. Acceptable values are "Obligate" and "Facultative".
 # - model.type: A string specifying the type of model to fit in the `betapart::decay.model` function. Acceptable values are 
 #   typically "poisson" or "binomial" depending on the model choice for fitting.
 #
 # Output:
 # - Prints the result of the Mantel test comparing Bray-Curtis dissimilarity and geographic distance (`abund_geo`).
 # - Returns a fitted decay model object from the `betapart::decay.model` function, representing the relationship between 
 #   abundance dissimilarity and geographic distance for the specified microbe type.
 #
 # Notes:
 # - The Bray-Curtis dissimilarity matrix is calculated based on the abundance data of 16S sequences.
 # - The geographic distance is computed using the Haversine formula, and the distances are converted to kilometers.
 # - The Mantel test is used to assess the statistical correlation between ecological dissimilarity and geographic distance 
 #   using the Spearman method.
 # - The function requires the `vegdist` function from the `vegan` package to calculate the Bray-Curtis dissimilarity and the 
 #   `distm` function from the `geosphere` package for geographic distance calculation.
 # - The decay model is fit using the `decay.model` function from the `betapart` package.
 # - The function filters the data based on the `WeightsObligateMicrobe` or `WeightsTransientMicrobe` column to select either 
 #   obligate or facultative microbes.
 #
 microbe_type_decay_model <- function(data, type, model.type, decay.type, log.dist=FALSE){
   #bray curtis dissimilarity matrix of 16s
   if(type == 'ObligateAll'){
     abund <- data %>%
       filter(WeightsObligateMicrobe == 1) %>%
       select(UniqueID, starts_with('16s')) %>%
       select(-UniqueID)
 
     #distance matrix of sites
     geo <- data %>%
       filter(WeightsObligateMicrobe == 1) %>%
       select(UniqueID, Long, Lat) %>%
       select(-UniqueID) %>%
       mutate()
   } else if (type == 'Facultative'){
   } else if (type == 'TransientAll'){
     abund <- data %>%
       filter(WeightsTransientMicrobe == 1) %>%
       select(UniqueID, starts_with('16s')) %>%
       select(-UniqueID)
 
     #distance matrix of sites
     geo <- data %>%
       filter(WeightsTransientMicrobe == 1) %>%
       select(UniqueID, Long, Lat) %>%
       select(-UniqueID) %>%
       mutate()
   } else if (type == 'ObligateSocial'){
     abund <- data %>%
       filter(WeightsObligateSocialMicrobe == 1) %>%
       select(UniqueID, starts_with('16s')) %>%
       select(-UniqueID)
     
     #distance matrix of sites
     geo <- data %>%
       filter(WeightsObligateSocialMicrobe == 1) %>%
       select(UniqueID, Long, Lat) %>%
       select(-UniqueID) %>%
       mutate()
   } else if (type == 'TransientSocial'){
     abund <- data %>%
       filter(WeightsTransientSocialMicrobe == 1) %>%
       select(UniqueID, starts_with('16s')) %>%
       select(-UniqueID)
     
     #distance matrix of sites
     geo <- data %>%
       filter(WeightsTransientSocialMicrobe == 1) %>%
       select(UniqueID, Long, Lat) %>%
       select(-UniqueID) %>%
       mutate()
   } else if (type == 'ObligateSolitary'){
     abund <- data %>%
       filter(WeightsObligateSolitaryMicrobe == 1) %>%
       select(UniqueID, starts_with('16s')) %>%
       select(-UniqueID)
     
     #distance matrix of sites
     geo <- data %>%
       filter(WeightsObligateSolitaryMicrobe == 1) %>%
       select(UniqueID, Long, Lat) %>%
       select(-UniqueID) %>%
       mutate()
   } else if (type == 'TransientSolitary'){
     abund <- data %>%
       filter(WeightsTransientSolitaryMicrobe == 1) %>%
       select(UniqueID, starts_with('16s')) %>%
       select(-UniqueID)
     
     #distance matrix of sites
     geo <- data %>%
       filter(WeightsTransientSolitaryMicrobe == 1) %>%
       select(UniqueID, Long, Lat) %>%
       select(-UniqueID) %>%
       mutate()
     }
   
 
   #abundance data frame - bray curtis dissimilarity
   dist.abund <- vegdist(abund, method = "bray")
   dist.abund <- 1 - dist.abund
 
   #geographic data frame - haversine distance in m (takes a df with lat and long and calculates dist)
   d.geo <- distm(geo, fun = distHaversine)
   dist.geo <- as.dist(d.geo)/1000
 
   #abundance vs geographic mantel test
   abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 999, na.rm = TRUE)
   print(abund_geo)
   #browser()
   if(model.type=="power"){
      dist.geo[dist.geo == 0] <- 0.0001
   }   
   dist_decay_model <- betapart::decay.model(dist.abund,
                                             dist.geo,
                                             y.type=decay.type,
                                             model.type = model.type,
                                             perm=999)
   dist_decay_model
 
 }
 
 # Function: plot_decay_ggplot_single
 #
 # Description:
 # This function creates a single ggplot visualization of a decay model fit. The plot displays the relationship between the 
 # independent and dependent variables along with the fitted decay curve. The scatter points represent the raw data, and the 
 # fitted curve is drawn to show the decay relationship. This is typically used to visualize the results of a decay model, 
 # such as those produced by the `betapart::decay.model` function.
 #
 # Parameters:
 # - x: A list or model object that contains the decay model output and data. This object must include:
 #   - x$data.x: A numeric vector or column representing the independent variable (e.g., geographic distance or other explanatory variable).
 #   - x$data.y: A numeric vector or column representing the dependent variable (e.g., Bray-Curtis dissimilarity or ecological distance).
 #   - x$model: A fitted model object (e.g., decay model) containing the model's predictions.
 # - xlab: A string representing the label for the x-axis (e.g., "Distance").
 # - ylab: A string representing the label for the y-axis (e.g., "Dissimilarity").
 # - mod1color: A string specifying the color of the fitted line (default is 'navy'). This color is used for the second line on the plot.
 # - col: A string specifying the color of the points (default is 'black'). This is used for the jittered points.
 # - lty: A string specifying the line type for the fitted curve (default is "solid"). This argument is not currently used but can be modified to support different line types.
 # - lwd: A numeric value specifying the line width of the fitted curve (default is 1.5).
 # - cex: A numeric value for the size of the scatter points (default is 1).
 #
 # Output:
 # - Returns a ggplot object representing the decay plot. The plot consists of:
 #   - Jittered points to represent the raw data (with customizable point color, size, and transparency).
 #   - A fitted decay curve (two lines: one black and one in the specified color, representing the model fit).
 #   - X and Y axis labels as defined by the `xlab` and `ylab` parameters.
 # 
 # Notes:
 # - The function assumes that `x` contains a fitted decay model and its corresponding data.
 # - The fitted model's predictions are used to create the fitted line, with the data ordered by the independent variable.
 # - The jittered points use the `geom_jitter` function for a visually appealing display of data points with a slight random noise to prevent overplotting.
 # - The plot is created using `ggplot2`, with the appearance customized by `theme_classic` and adjustments to axis labels and text size.
 #
 plot_decay_ggplot_single <- function(x,
                                      xlab,
                                      ylab,
                                      mod1color='navy',
                                      col = "black",
                                      lty = "solid",
                                      lwd = 1.5,
                                      cex = 1) {
 
   # Extract data and fitted values
   data <- data.frame(x$data.x, x$data.y)
   model <- x$model
   fitted_values <- data.frame(fitted(model)[order(data$x.data.x)])
   sorted_data <- data[order(data$x.data.x), ]
 
   # Create the ggplot object
   p <- ggplot(data, aes(x = x.data.x, y = x.data.y)) +
     geom_point(fill = mod1color, alpha = 0.1 , color="black", pch=21, cex=3) +
     geom_line(aes(x = sorted_data$x.data.x, y = fitted_values$fitted.model..order.data.x.data.x..),
               color = 'black', linewidth=2.5,) +
     geom_line(aes(x = sorted_data$x.data.x, y = fitted_values$fitted.model..order.data.x.data.x..), color = mod1color, linewidth=2) +
     labs(x=xlab, y=ylab) +
     theme_classic() +
     theme(axis.title.x = element_text(size=16),
           axis.title.y = element_text(size=16),
           text = element_text(size=16))
 
   return(p)
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
 plot_decay_ggplot_combined <- function(x,
                                        z,
                                        xlab,
                                        ylab,
                                        mod1color='navy',
                                        mod2color='gold',
                                        alpha1 = 0.1,
                                        alpha2 = 0.5,
                                        col = "black",
                                        lty1,
                                        lty2,
                                        lwd = 1.5,
                                        cex = 1,
                                        add.points=TRUE,
                                        log.dist=FALSE) {
 
   # Extract data and fitted values
   data1 <- data.frame(x$data.x, x$data.y)
   model1 <- x$model
   fitted_values1 <- data.frame(fitted(model1)[order(data1$x.data.x)])
   sorted_data1 <- data1[order(data1$x.data.x), ]
 
   # Extract data and fitted values
   data2 <- data.frame(z$data.x, z$data.y)
   model2 <- z$model
   fitted_values2 <- data.frame(fitted(model2)[order(data2$z.data.x)])
   sorted_data2 <- data2[order(data2$z.data.x), ]
 
   if(add.points==TRUE){
   # Create the ggplot object
   p <- ggplot(data1, aes(x = log(x.data.x), y = x.data.y)) +
     geom_point(data=data1, aes(x = log(x.data.x), y = (1-x.data.y)), fill = mod1color, alpha = alpha1 , color="black", pch=21, cex=3, position = position_jitter(w = 10, h = 0)) +
     geom_point(data=data2, aes(x = log(z.data.x), y = (1-z.data.y)), fill = mod2color, alpha = alpha2 , color="black", pch=21, cex=3, position = position_jitter(w = 10, h = 0)) +
     geom_line(aes(x = log(sorted_data1$x.data.x), y = (1 - fitted_values1$fitted.model1..order.data1.x.data.x..)),
               color = 'black', linewidth=2.5,) +
     geom_line(aes(x = log(sorted_data1$x.data.x), y = (1 - fitted_values1$fitted.model1..order.data1.x.data.x..)), color = mod1color, linewidth=2, linetype=lty1) +
     geom_line(data=data2, aes(x = log(sorted_data2$z.data.x), y = (1 - fitted_values2$fitted.model2..order.data2.z.data.x..)),
               color = 'black', linewidth=2.5,) +
     geom_line(data=data2, aes(x = log(sorted_data2$z.data.x), y = (1 - fitted_values2$fitted.model2..order.data2.z.data.x..)), color = mod2color, linewidth=2, linetype=lty2) +
     labs(x=xlab, y=ylab) +
     theme_classic() +
     theme(axis.title.x = element_text(size=16),
           axis.title.y = element_text(size=16),
           text = element_text(size=16))
   }
   
   if(log.dist==TRUE){
      # Create the ggplot object
      p <- ggplot(data1, aes(x = x.data.x, y = x.data.y)) +
         geom_point(data=data1, aes(x = x.data.x, y = x.data.y), fill = mod1color, alpha = alpha1 , color="black", pch=21, cex=3, position = position_jitter(w = 0.01, h = 0)) +
         geom_point(data=data2, aes(x = z.data.x, y = z.data.y), fill = mod2color, alpha = alpha2 , color="black", pch=21, cex=3, position = position_jitter(w = 0.01, h = 0)) +
         geom_line(aes(x = sorted_data1$x.data.x, y = fitted_values1$fitted.model1..order.data1.x.data.x..),
                   color = 'black', linewidth=2.5,) +
         geom_line(aes(x = sorted_data1$x.data.x, y = fitted_values1$fitted.model1..order.data1.x.data.x..), color = mod1color, linewidth=2, linetype=lty1) +
         geom_line(data=data2, aes(x = sorted_data2$z.data.x, y = fitted_values2$fitted.model2..order.data2.z.data.x..),
                   color = 'black', linewidth=2.5,) +
         geom_line(data=data2, aes(x = sorted_data2$z.data.x, y = fitted_values2$fitted.model2..order.data2.z.data.x..), color = mod2color, linewidth=2, linetype=lty2) +
         labs(x=xlab, y=ylab) +
         scale_x_log10()+
            #limits = c(60, 350)) +  # Adjust upper limit based on your data
         #   breaks = c(1, 10, 100, 1000),
         #   labels = scales::label_number()) +
         theme_classic() +
         theme(axis.title.x = element_text(size=16),
               axis.title.y = element_text(size=16),
               text = element_text(size=16)) 
   }
 
   if(add.points==FALSE){
     # Create the ggplot object
     p <- ggplot(data1, aes(x = x.data.x, y = x.data.y)) + 
       geom_line(aes(x = sorted_data1$x.data.x, y = (1 - fitted_values1$fitted.model1..order.data1.x.data.x..)),
                 color = 'black', linewidth=2.5,) +
       geom_line(aes(x = sorted_data1$x.data.x, y = (1 - fitted_values1$fitted.model1..order.data1.x.data.x..)), color = mod1color, linewidth=2, linetype=lty1) +
       geom_line(data=data2, aes(x = sorted_data2$z.data.x, y = (1 - fitted_values2$fitted.model2..order.data2.z.data.x..)),
                 color = 'black', linewidth=2.5,) +
       geom_line(data=data2, aes(x = sorted_data2$z.data.x, y = (1 - fitted_values2$fitted.model2..order.data2.z.data.x..)), color = mod2color, linewidth=2, linetype=lty2) +
       labs(x=xlab, y=ylab) +
       theme_classic() +
       theme(axis.title.x = element_text(size=16),
             axis.title.y = element_text(size=16),
             text = element_text(size=16)) + ylim(0,1)
   }
 
   return(p)
 
 }
 
 