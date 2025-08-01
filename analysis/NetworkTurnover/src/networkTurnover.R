# Module: Network Turnover Analysis
#
# Description:
# This script provides a suite of functions for analyzing and visualizing network turnover in ecological networks. It includes 
# functions for preparing obligate and transient networks, running network turnover models, and plotting the results. 
# Specifically, the script is designed to:
# 1. Prepare the obligate and transient networks using the `prep_obligate_network` and `prep_transient_network` functions.
# 2. Fit network turnover models using the `run_network_turnover_mod` function.
# 3. Generate plots of the model results, both for individual network types using the `plot_network_turnover_mod_single` function, 
#    and for comparisons between two network types using the `plot_network_turnover_mod_compare` function.
# 4. Visualize the conditional effects of the models, including credible intervals (95%, 80%, and 50%), using ribbons and line plots.
#
# Functions:
# 1. **prep_obligate_network**: Prepares the obligate network data, organizing and filtering it for use in turnover models.
#    - Input: Raw network data (e.g., species interactions or ecological variables).
#    - Output: Processed data ready for modeling.
#
# 2. **prep_transient_network**: Prepares the transient network data, organizing and filtering it for use in turnover models.
#    - Input: Raw network data for transient interactions or variables.
#    - Output: Processed data ready for modeling.
#
# 3. **run_network_turnover_mod**: Fits a network turnover model using the `brms` package to analyze turnover or dissimilarity 
#    in ecological networks. The model is fitted for either obligate or transient networks, depending on the input data.
#    - Input: Prepared network data and model specifications (e.g., formula, priors).
#    - Output: A fitted model object.
#
# 4. **plot_network_turnover_mod_single**: Generates a plot of the network turnover model results for a single network type 
#    ("Obligate" or "Transient"). It visualizes the conditional effects of the model, including credible intervals (95%, 80%, 
#    and 50%) using ribbons and line plots. It also customizes the plot appearance based on the significance of the geographic 
#    distance effect in the model.
#    - Input: A fitted model, network data, effect of interest, response variable, and a label for the y-axis.
#    - Output: A ggplot2 plot object with the conditional effects and credible intervals.
#
# 5. **plot_network_turnover_mod_compare**: Generates a comparative plot of the network turnover model results for two network 
#    types ("Obligate" and "Transient"). It visualizes the conditional effects of both models, including credible intervals (95%, 
#    80%, and 50%) using ribbons and line plots. The plot also includes points for original data and customizes the plot 
#    appearance based on the significance of the geographic distance effect in each model.
#    - Input: Two fitted models, network data for each network type, effect of interest, response variable, and a label for the 
#      y-axis.
#    - Output: A ggplot2 plot object comparing the conditional effects of both models, with ribbons for credible intervals and 
#      lines for estimates, along with original data points.
#
# Notes:
# - All plotting functions visualize the conditional effects from the fitted models, focusing on the relationship between 
#   geographic distance and the response variable.
# - The plot functions use ribbons to represent 95%, 80%, and 50% credible intervals and customize their appearance 
#   based on the significance of the geographic distance effect (`GeoDist`).
# - The color of the ribbons and points depends on the network type ("Obligate" or "Transient"), with each type having a 
#   distinct color (e.g., dark green for obligate, dark orange for transient).
# - The script assumes that the input data are preprocessed and formatted correctly for modeling and visualization.
#
# Usage:
# 1. Prepare the data using the `prep_obligate_network` and `prep_transient_network` functions.
# 2. Fit the network turnover models using the `run_network_turnover_mod` function.
# 3. Generate individual plots for each network type using the `plot_network_turnover_mod_single` function.
# 4. Compare the results of the two network types using the `plot_network_turnover_mod_compare` function.


# Function: prep_obligate_network
#
# Description:
# Prepares a network containing only obligate symbionts from the raw network data.
# The function filters out interactions for specific obligate symbionts from a given network 
# and creates a new list of networks containing only those obligate symbionts.
#
# Parameters:
# - raw_network: A list of networks, where each network corresponds to a site and 
#   contains interactions (typically species or bacteria) as rows and columns. Default is `spNet_micro`.
#
# Output:
# - Returns a list of networks where each network contains only interactions involving 
#   specific obligate symbionts (such as "Lactobacillus", "Bifidobacterium", etc.).
#
# Example:
# # Example call:
# obligate_network <- prep_obligate_network(raw_network = my_network_data)
#
# Notes:
# - The function filters for specific symbionts listed in `these_obligates`.
# - Each site in the original `raw_network` is processed individually, and a new network 
#   containing only the obligate symbionts is created.
# - The output is a list of networks where each site contains only the interactions of obligate species.
#
prep_obligate_network <- function(raw_network=spNet_micro, these_obligates, genera_to_keep=NULL){
  
  bee.obligates <- paste(these_obligates, collapse = "|")
  
  site_list <- names(raw_network)
  
  only_obligate_network <- list()
  
  for (x in site_list){
    
    obligates_rows <- rownames(raw_network[[x]])
    
    ob_rows_to_keep <- grep(bee.obligates, obligates_rows)
    
    if (is.null(genera_to_keep)){
      ob_new_net <- raw_network[[x]][ob_rows_to_keep, , drop=FALSE]
      } else {
        
        keep_these_genera <- grep(paste(genera_to_keep, collapse = "|"), colnames(raw_network[[x]]))
        #browser()
        ob_new_net <- raw_network[[x]][ob_rows_to_keep,keep_these_genera, drop=FALSE]
    }
    
    print(dim(ob_new_net))
    new_name <- x
    
    only_obligate_network[[new_name]] <- as.matrix(ob_new_net)
    
  }
  only_obligate_network
  
}

# Function: prep_transient_network
#
# Description:
# Prepares a network containing only transient symbionts from the raw network data.
# The function filters out interactions for specific obligate symbionts and retains the transient symbionts,
# creating a new list of networks with only transient symbionts for each site.
#
# Parameters:
# - raw_network: A list of networks, where each network corresponds to a site and contains interactions
#   (typically species or bacteria) as rows and columns. Default is `spNet_micro`.
#
# Output:
# - Returns a list of networks where each network contains only interactions involving transient symbionts 
#   (i.e., species or bacteria that are not obligate symbionts).
#
# Notes:
# - The function filters out rows corresponding to obligate symbionts and retains the transient ones.
# - The output is a list of networks with transient symbionts for each site.
#
prep_transient_network <- function(raw_network=spNet_micro, these_obligates, genera_to_keep=NULL){
  site_list <- names(raw_network)
  
  ## obligate symbionts
  bee.obligates <- paste(these_obligates, collapse = "|")
  
  only_transient_network <- list()
  
  for (x in site_list){
    
    trans_rows <- rownames(raw_network[[x]])
    
    trans_rows_to_drop <- !grepl(bee.obligates, trans_rows)
    
    trans_new_net <- raw_network[[x]][trans_rows_to_drop, , drop=FALSE]
    
    if (is.null(genera_to_keep)){
      trans_new_net <- raw_network[[x]][trans_rows_to_drop,]
    } else {
      
      keep_these_genera <- grep(paste(genera_to_keep, collapse = "|"), colnames(raw_network[[x]]))
      #browser()
      trans_new_net <- raw_network[[x]][trans_rows_to_drop,keep_these_genera, drop=FALSE]
    }
    
    
    print(dim(trans_new_net))
    new_name <- x
    
    only_transient_network[[new_name]] <- as.matrix(trans_new_net)
  }
  only_transient_network
}

find_sites_for_betalinkr <- function(these_networks){
  microbe_nets <- lapply(these_networks, function(x){
    #net_dims <- dim()
    ifelse(dim(x)[1]==0&dim(x)[2]==0, 0, 1)
  })
  
  sites_with_data <- names(which(unlist(microbe_nets) == 1))
  
  nets_to_include <- these_networks[names(these_networks) %in% sites_with_data]
  
  ## drops unresolved species
  cleaned_obligate_network <- lapply(nets_to_include, function(df) {
    df[, colnames(df) != "", drop = FALSE]
  })
  
  print(names(cleaned_obligate_network))
  
  ## splits network list elements into their own objects
  return(list2env(cleaned_obligate_network, envir = .GlobalEnv))
}



fix_betalinkr_output <- function(betalinkr_output){
  
  lower.order <- "Microbes"
  higher.order <- "Pollinators"
  #View(obligate_poll_betalink)
  
  colnames(betalinkr_output) <- c("Site1",
                                  "Site2",
                                  "DissimilaritySpeciesComposition",
                                  "OnlySharedLinks",
                                  "WholeNetworkLinks",
                                  "SpeciesTurnoverLinks",
                                  paste("TurnoverAbsence",lower.order,sep=""),
                                  paste("TurnoverAbsence",higher.order,sep=""),
                                  "TurnoverAbsenceBoth")
  
  #browser()
  geo <- spec.net %>%
    filter(Site %in% betalinkr_output$Site1 | Site %in% betalinkr_output$Site2) %>%
    select(Site, Long, Lat) %>%
    distinct()
  
  site_coords1 <- geo %>%
    rename(Site1 = Site, Long1 = Long, Lat1 = Lat)
  
  site_coords2 <- geo %>%
    rename(Site2 = Site, Long2 = Long, Lat2 = Lat)
  
  dissim_df <- betalinkr_output %>%
    left_join(site_coords1, by = "Site1") %>%
    left_join(site_coords2, by = "Site2")
  
  
  betalinkr_output_temp <- dissim_df %>%
    mutate(GeoDist = geosphere::distHaversine(cbind(Long2, Lat2), cbind(Long1, Lat1))/1000)
    
  
  return(betalinkr_output_temp)
}


# Function: run_network_turnover_mod
#
# Description:
# Runs a network turnover model using the `brms` package to model the relationship between the turnover 
# of interactions in a network and geographic distance, with site-level random effects.
#
# Parameters:
# - this_component: The name of the variable (column) in the network data that represents the network component
#   to be modeled (e.g., a dissimilarity measure or a turnover value).
# - this_network: The data frame containing the network data, which includes the specified component and other
#   relevant variables (e.g., `GeoDist`, `Site1`, `Site2`).
#
# Output:
# - Returns a `brms` model object after fitting the network turnover model.
#
# Notes:
# - The function requires a column in the network data corresponding to the component of interest (e.g., `ST`).
# - The model predicts this component as a function of geographic distance (`GeoDist`), with random effects for 
#   site pairings (`Site1`, `Site2`).
# - The model is fit using the `brm` function from the `brms` package.
#
run_network_turnover_mod <- function(this_component,
                                     this_network){
  # Assign the value of 'this_component' to a new variable 'y'.
  y <- this_component
  
  this_network$this_component <- this_network[[this_component]]
  
  # Define a formula for brms with this_component the response variable,
  # 'GeoDist' as a fixed effect, and 'Site1' and 'Site2' as random effects.
  forms <- bf(formula(this_component~GeoDist + (1|Site1) + (1|Site2)))
  
  # Fit model
  mod1 <-  brm(forms, this_network,
               cores=1,
               iter = 10^5,
               chains = 3,
               thin=1,
               init=0,
               control = list(adapt_delta = 0.99),
               save_pars = save_pars(all = TRUE))
  
  mod1
}

run_all_turnover_mods <- function(run.mods=FALSE, # TRUE if never ran model before, false if you want to load models
                                  ob.net=NULL, # Null by default, if run.mods==TRUE input obligate network here
                                  trans.net=NULL, # Null by default, if run.mods==TRUE input transient network here
                                  filepath # if run.mods=TRUE, input desired save filepath, otherwise input the filepath to load model results
                                  ){
  if (run.mods==TRUE){
    
    ## Interaction turnover
    int.obligate.mod <- run_network_turnover_mod(this_component="WholeNetworkLinks",
                                                 this_network=ob.net)
    
    int.transient.mod <- run_network_turnover_mod(this_component="WholeNetworkLinks",
                                                  this_network=trans.net)
    
    ## Turnover species composition
    speccomp.obligate.mod <- run_network_turnover_mod(this_component="DissimilaritySpeciesComposition",
                                                      this_network=ob.net)
    
    speccomp.transient.mod <- run_network_turnover_mod(this_component="DissimilaritySpeciesComposition",
                                                       this_network=trans.net)
    
    
    ## Rewiring
    rewiring.obligate.mod <- run_network_turnover_mod(this_component="OnlySharedLinks",
                                                      this_network=ob.net)
    
    rewiring.transient.mod <- run_network_turnover_mod(this_component="OnlySharedLinks",
                                                       this_network=trans.net)
    
    ## Host-driven turnover
    host.driven.obligate.mod <- run_network_turnover_mod(this_component="TurnoverAbsencePollinators",
                                                         this_network=ob.net)
    
    host.driven.transient.mod <- run_network_turnover_mod(this_component="TurnoverAbsencePollinators",
                                                          this_network=trans.net)
    
    ## Microbe-driven turnover
    microbe.driven.obligate.mod <- run_network_turnover_mod(this_component="TurnoverAbsenceMicrobes",
                                                            this_network=ob.net)
    
    microbe.driven.transient.mod <- run_network_turnover_mod(this_component="TurnoverAbsenceMicrobes",
                                                             this_network=trans.net)
    
    ## Complete turnover
    complete.obligate.mod <- run_network_turnover_mod(this_component="TurnoverAbsenceBoth",
                                                      this_network=ob.net)
    
    complete.transient.mod <- run_network_turnover_mod(this_component="TurnoverAbsenceBoth",
                                                       this_network=trans.net)
    
    ## save out models
    save(int.obligate.mod,
         int.transient.mod,
         speccomp.obligate.mod,
         speccomp.transient.mod,
         rewiring.obligate.mod,
         rewiring.transient.mod,
         host.driven.obligate.mod,
         host.driven.transient.mod,
         microbe.driven.obligate.mod,
         microbe.driven.transient.mod,
         complete.obligate.mod,
         complete.transient.mod,
         file=filepath) ## TODO update with correct filepath
  } else { 
    load(filepath) ## TODO update with correct filepath
  }
}



# Function: plot_network_turnover_mod_single
#
# Description:
# This function generates a plot of the network turnover model results, specifically for a single network type (either 
# "Obligate" or "Transient"). It visualizes the conditional effects of the model, including credible intervals (95%, 80%, 
# and 50%) using ribbons and line plots. It also customizes the plot appearance based on the significance of the 
# geographic distance effect in the model.
#
# Parameters:
# - mod1: A fitted `brms` model object containing the network turnover model results.
# - this.network: A data frame containing the network data, which must include the variables used in the model (e.g., 
#   geographic distance and response variable).
# - network_type: A string indicating the type of network to plot. Acceptable values are "Obligate" or "Transient".
# - this.effect: A string specifying the effect of interest to plot from the conditional effects (usually geographic distance).
# - this.resp: A string specifying the response variable from the model (e.g., turnover or dissimilarity).
# - label: A string to label the y-axis in the plot.
#
# Output:
# - A list containing:
#   1. A ggplot2 plot object representing the model's conditional effects, with ribbons for credible intervals and lines 
#      for estimates.
#   2. A summary of the model's effect of geographic distance (`model_geodist`), including its posterior probability 
#      (`Pgt0`).
#
# Notes:
# - The function assigns different colors for the plot based on the network type ("Obligate" or "Transient").
# - Credible intervals are visualized with three ribbons: 95%, 80%, and 50% intervals, each with different transparency 
#   and line styles (dotted, dashed, and solid).
# - The geographic distance effect (`GeoDist`) is checked for significance, and if it is highly significant (Pgt0 >= 0.95 
#   or Pgt0 <= 0.05), the ribbon color is set to match the network type.
# - The plot is customized with labels, axis titles, and classic ggplot2 themes.
#
plot_network_turnover_mod_single <- function(mod1, 
                                      this.network,
                                      network_type,
                                      this.effect,
                                      this.resp,
                                      label
                                      ){
  

  mod_summary <- write.summ.table(mod1)
  model_geodist <- mod_summary[rownames(mod_summary) == "GeoDist",]
  
  if(network_type == "Obligate") {
    point_color <- "darkgreen"
    if(model_geodist$Pgt0 >= 0.95){
      ribbon_color <- "darkgreen"
    } else if (model_geodist$Pgt0 <= 0.05) {
      ribbon_color <- "darkgreen"
    } else {ribbon_color <- NA}
  }
  
  if(network_type == "Transient") {
    point_color <- "darkorange"
    if(model_geodist$Pgt0 >= 0.95){
      ribbon_color <- "darkorange"
    } else if (model_geodist$Pgt0 <= 0.05) {
      ribbon_color <- "darkorange"
    } else {ribbon_color <- NA}
  }
  
  # Extract the data from conditional_effects
  cond_effects_data <- conditional_effects(mod1, effects = this.effect, resp = this.resp, plot = FALSE)
  plot_data <- cond_effects_data[[this.effect]]
  
  # Plot using ggplot2 for credible intervals with geom_ribbon
  plot_obj <- ggplot(plot_data, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, 
                fill = ribbon_color,
                color = point_color, linetype='dotted') +
    geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__),
                    ymax = upper__ - 0.1 * (upper__ - lower__)),
                alpha = 0.3, 
                fill=ribbon_color, 
                color = point_color, linetype='dashed') +
    geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__),
                    ymax = upper__ - 0.25 * (upper__ - lower__)),
                alpha = 0.4, 
                fill=ribbon_color,
                color = point_color, linetype='solid') +
    #Add line for the estimates
    geom_line(data = plot_data, color = 'black', linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data, color = point_color, linewidth=2, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add points for original data
    geom_point(data = this.network, aes(x = .data[[this.effect]], y = .data[[this.resp]]),
               fill = point_color, alpha = 0.9,color="black", pch=21, cex=3) +
    # Labels and theme
    theme_classic()  +
    labs(x = "Geographic Distance (km)", y = label,
         fill = "Credible Interval") +
    theme_classic() +
    ylim(0,1) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none")
  
  return(list(plot_obj, model_geodist))
}

# Function: plot_network_turnover_mod_compare
#
# Description:
# This function generates a comparative plot of the network turnover model results for two network types (either 
# "Obligate" or "Transient"). It visualizes the conditional effects of both models, including credible intervals (95%, 
# 80%, and 50%) using ribbons and line plots. The plot also includes points for original data and customizes the plot 
# appearance based on the significance of the geographic distance effect in each model.
#
# Parameters:
# - mod1: A fitted `brms` model object containing the network turnover model results for the first network type.
# - mod2: A fitted `brms` model object containing the network turnover model results for the second network type.
# - this.network1: A data frame containing the network data for the first network type, which must include the 
#   variables used in the model (e.g., geographic distance and response variable).
# - this.network2: A data frame containing the network data for the second network type, which must include the 
#   variables used in the model.
# - network_type1: A string indicating the type of the first network to plot. Acceptable values are "Obligate" or "Transient".
# - network_type2: A string indicating the type of the second network to plot. Acceptable values are "Obligate" or "Transient".
# - this.effect: A string specifying the effect of interest to plot from the conditional effects (usually geographic distance).
# - this.resp: A string specifying the response variable from the model (e.g., turnover or dissimilarity).
# - label: A string to label the y-axis in the plot.
#
# Output:
# - A list containing:
#   1. A ggplot2 plot object representing the conditional effects of both models, with ribbons for credible intervals 
#      and lines for estimates, along with original data points.
#   2. A summary of the geographic distance effect (`GeoDist`) for both models, including their posterior probabilities 
#      (`Pgt0`), with the results combined into a single data frame.
#
# Notes:
# - The function assigns different colors for the plot based on the network type ("Obligate" or "Transient").
# - Credible intervals are visualized with three ribbons: 95%, 80%, and 50% intervals, each with different transparency 
#   and line styles (dotted, dashed, and solid).
# - The geographic distance effect (`GeoDist`) is checked for significance, and if it is highly significant (Pgt0 >= 0.95 
#   or Pgt0 <= 0.05), the ribbon color is set to match the network type.
# - The plot includes the original data points for each network type, with different point colors depending on the network.
# - The plot is customized with labels, axis titles, and classic ggplot2 themes.
plot_network_turnover_mod_compare <- function(mod1,
                                              mod2,
                                              this.network1,
                                              this.network2,
                                              network_type1,
                                              network_type2,
                                              this.effect,
                                              this.resp,
                                              label
){
  
  # prep cond effects mod 1
  
  mod_summary1 <- write.summ.table(mod1)
  model_geodist1 <- mod_summary1[rownames(mod_summary1) == "GeoDist",]
  
  if(network_type1 == "Obligate") {
    point_color1 <- "darkgreen"
    if(model_geodist1$Pgt0 >= 0.95){
      ribbon_color1 <- "darkgreen"
      linetype1 <- 'solid'
    } else if (model_geodist1$Pgt0 <= 0.05) {
      ribbon_color1 <- "darkgreen"
      linetype1 <- 'solid'
    } else {
      ribbon_color1 <- NA
      linetype1 <- 'dashed'
    }
  }
  
  if(network_type1 == "Transient") {
    point_color1 <- "darkorange"
    if(model_geodist1$Pgt0 >= 0.95){
      ribbon_color1 <- "darkorange"
      linetype1 <- 'solid'
    } else if (model_geodist1$Pgt0 <= 0.05) {
      ribbon_color1 <- "darkorange"
      linetype1 <- 'solid'
    } else {
      ribbon_color1 <- NA
      linetype1 <- 'dashed'
      }
  }
  
  # Extract the data from conditional_effects
  cond_effects_data1 <- conditional_effects(mod1, effects = this.effect, resp = this.resp, plot = FALSE)
  plot_data1 <- cond_effects_data1[[this.effect]]
  
  # now repeat for mod 2
  mod_summary2 <- write.summ.table(mod2)
  model_geodist2 <- mod_summary2[rownames(mod_summary2) == "GeoDist",]
  
  if(network_type2 == "Obligate") {
    point_color2 <- "darkgreen"
    if(model_geodist2$Pgt0 >= 0.95){
      ribbon_color2 <- "darkgreen"
      linetype2 <- 'solid'
    } else if (model_geodist2$Pgt0 <= 0.05) {
      ribbon_color2 <- "darkgreen"
      linetype2 <- 'solid'
    } else {
      ribbon_color2 <- NA
      linetype2 <- 'dashed'
      }
  }
  
  if(network_type2 == "Transient") {
    point_color2 <- "darkorange"
    if(model_geodist2$Pgt0 >= 0.95){
      ribbon_color2 <- "darkorange"
      linetype2 <- 'solid'
    } else if (model_geodist2$Pgt0 <= 0.05) {
      ribbon_color2 <- "darkorange"
      linetype2 <- 'solid'
    } else {
      ribbon_color2 <- NA
      linetype2 <- 'dashed'
      }
  }
  
 
  
  # Extract the data from conditional_effects
  cond_effects_data2 <- conditional_effects(mod2, effects = this.effect, resp = this.resp, plot = FALSE)
  plot_data2 <- cond_effects_data2[[this.effect]]
  
  # Plot using ggplot2 for credible intervals with geom_ribbon
  plot_obj <- ggplot(plot_data1, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.5, 
                fill = ribbon_color1,
                color = point_color1#, linetype='dotted'
                ) +
    # geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__),
    #                 ymax = upper__ - 0.1 * (upper__ - lower__)),
    #             alpha = 0.3, 
    #             fill=ribbon_color1, 
    #             color = point_color1
    #             #, linetype='dashed'
    #             ) +
    # geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__),
    #                 ymax = upper__ - 0.25 * (upper__ - lower__)),
    #             alpha = 0.4, 
    #             fill=ribbon_color1,
    #             color = point_color1
    #             #, linetype='solid'
    #             ) +
    #   # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(data=plot_data2, aes(ymin = lower__, ymax = upper__), alpha = 0.5, 
                fill = ribbon_color2,
                color = point_color2
                #, linetype='dotted'
                ) +
    # geom_ribbon(data=plot_data2, aes(ymin = lower__ + 0.1 * (upper__ - lower__),
    #                 ymax = upper__ - 0.1 * (upper__ - lower__)),
    #             alpha = 0.3, 
    #             fill=ribbon_color2, 
    #             color = point_color2
    #             #, linetype='dashed'
    #             ) +
    # geom_ribbon(data=plot_data2, aes(ymin = lower__ + 0.25 * (upper__ - lower__),
    #                 ymax = upper__ - 0.25 * (upper__ - lower__)),
    #             alpha = 0.4, 
    #             fill=ribbon_color2,
    #             color = point_color2#, linetype='solid'
    #             ) +
    #   # Add points for original data
    geom_point(data = this.network1, aes(x = .data[[this.effect]], y = .data[[this.resp]]),
               fill = point_color1, alpha = 0.5, pch=21, cex=3) +
    # Add points for original data
    geom_point(data = this.network2, aes(x = .data[[this.effect]], y = .data[[this.resp]]),
               fill = point_color2, alpha = 0.5, pch=21, cex=3) +
    #Add line for the estimates
    #geom_line(data = plot_data1, color = 'black', linetype=linetype1, linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data1, color = point_color1, linetype=linetype1, linewidth=3, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    #geom_line(data = plot_data2, color = 'black', linetype=linetype2, linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data2, color = point_color2, linetype=linetype2, linewidth=3, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    theme_classic()  +
    labs(x = "Geographic Distance (km)", y = label,
         fill = "Credible Interval") +
    theme_classic() +
    ylim(-0.25,1.25) +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none")
  
  combined_mods <- bind_rows(model_geodist1, model_geodist2)
  
  rownames(combined_mods) <- c(paste(network_type1,label), paste(network_type2,label))
  
  return(list(plot_obj, combined_mods))
}
