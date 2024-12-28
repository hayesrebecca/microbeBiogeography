# Function: plot_model_condeff_compare
# 
# Description:
# Creates a comparative plot of conditional effects from two Bayesian models, 
# including credible intervals and observed data points.
# 
# Parameters:
# - model.a, model.b: Bayesian models (from `brms`) to compare.
# - this.effect: Predictor variable to plot.
# - this.resp.a, this.resp.b: Response variables corresponding to `model.a` and `model.b`.
# - point.data.a, point.data.b: Dataframes with observed data for `model.a` and `model.b`.
# - axis.breaks: Numeric vector for x-axis tick positions.
# - axis.labs: Labels for x-axis ticks.
# - xlabel, ylabel: X and Y-axis labels for the plot.
# - mod1color, mod2color: Colors for `model.a` and `model.b` plot elements. 
#   Defaults are 'navy' and 'coral'.
# - fill.a, fill.b: Logical flags to fill the credible interval ribbons for 
#   `model.a` and `model.b`.
# 
# Output:
# Returns a ggplot object with:
# - Conditional effect lines for both models.
# - 95%, 80%, and 50% credible intervals as shaded ribbons.
# - Points for observed data.
# 
plot_model_condeff_compare <- function(model.a,
                                       model.b,
                                       this.effect,
                                       this.resp.a, 
                                       this.resp.b,
                                       point.data.a,
                                       point.data.b,
                                       axis.breaks,
                                       axis.labs,
                                       xlabel,
                                       ylabel,
                                       mod1color='navy',
                                       mod2color='coral',
                                       fill.a=FALSE,
                                       fill.b=FALSE
){
  
  ## PD obligate ~ Bee Diversity
  
  # Extract the data from conditional_effects
  cond_effects_data_a <- conditional_effects(model.a, effects = this.effect, resp = this.resp.a, plot = FALSE)
  plot_data_a <- cond_effects_data_a[[paste(this.resp.a,".",this.resp.a, "_",this.effect, sep='')]]
  
  ## fixing col names TODO: should fix the col names in prep
  if ('PD.obligate' %in% colnames(plot_data_a)){
    plot_data_a$PDobligate <- plot_data_a$PD.obligate
  }
  if ('PD.obligate.log' %in% colnames(plot_data_a)){
    plot_data_a$PDobligatelog <- plot_data_a$PD.obligate.log
  }
  if ('PD.transient' %in% colnames(plot_data_a)){
    plot_data_a$PDtransient <- plot_data_a$PD.transient
  }
  if ('PD.transient.log' %in% colnames(plot_data_a)){
    plot_data_a$PDtransientlog <- plot_data_a$PD.transient.log
  }
  
  # Extract the data from conditional_effects
  cond_effects_data_b <- conditional_effects(model.b, effects = this.effect, resp = this.resp.b, plot = FALSE)
  plot_data_b <- cond_effects_data_b[[paste(this.resp.b, ".", this.resp.b,"_",this.effect, sep='')]]
  
  ## fixing col names TODO: should fix the col names in prep
  if ('PD.obligate' %in% colnames(plot_data_b)){
    plot_data_b$PDobligate <- plot_data_b$PD.obligate
  }
  if ('PD.obligate.log' %in% colnames(plot_data_b)){
    plot_data_b$PDobligatelog <- plot_data_b$PD.obligate.log
  }
  if ('PD.transient' %in% colnames(plot_data_b)){
    plot_data_b$PDtransient <- plot_data_b$PD.transient
  }
  if ('PD.transient.log' %in% colnames(plot_data_b)){
    plot_data_b$PDtransientlog <- plot_data_b$PD.transient.log
  }
  
  ## model a fill
  if (fill.a==TRUE){
    mod1fill=mod1color
  } else {mod1fill=NA}
  
  ## model b fill
  if (fill.b==TRUE){
    mod2fill=mod2color
  } else {mod2fill=NA}
  
  
  # Plot using ggplot2 for credible intervals with geom_ribbon
  plot_obj <- ggplot(plot_data_a, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, 
                fill=mod1fill,
                color = mod1color, linetype='dotted') +
    geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__),
                    ymax = upper__ - 0.1 * (upper__ - lower__)),
                alpha = 0.3, 
                fill=mod1fill, 
                color = mod1color, linetype='dashed') +
    geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__),
                    ymax = upper__ - 0.25 * (upper__ - lower__)),
                alpha = 0.4, 
                fill=mod1fill,
                color = mod1color, linetype='solid') +
    # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(data=plot_data_b, aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill=mod2fill, color = mod2color, linetype='dotted') +
    geom_ribbon(data=plot_data_b, aes(ymin = lower__ + 0.1 * (upper__ - lower__),
                                      ymax = upper__ - 0.1 * (upper__ - lower__)),
                alpha = 0.3, fill=mod2fill, color = mod2color, linetype='dashed') +
    geom_ribbon(data=plot_data_b, aes(ymin = lower__ + 0.25 * (upper__ - lower__),
                                      ymax = upper__ - 0.25 * (upper__ - lower__)),
                alpha = 0.4, fill=mod2fill, color = mod2color, linetype='solid') +
    #Add line for the estimates
    geom_line(data = plot_data_a, color = 'black', linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data_a, color = mod1color, linewidth=2, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data_b, color = 'black', linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add line for the estimates
    geom_line(data=plot_data_b, linewidth=2, aes(x = .data[[this.effect]],  y = .data$estimate__), color = mod2color) +
    # Add points for original data
    geom_point(data = point.data.a, aes(x = .data[[this.effect]], y = .data[[this.resp.a]]),
               fill = mod1color, alpha = 0.9,color="black", pch=21, cex=3) +
    # Add points for original data
    geom_point(data = point.data.b, aes(x = .data[[this.effect]], y = .data[[this.resp.b]]),
               fill = mod2color, alpha = 0.9, color="black", pch=21, cex=3) +
    coord_cartesian(xlim = if_else(range(point.data.a[[this.effect]]) > range(point.data.b[[this.effect]]), range(point.data.a[[this.effect]]), range(point.data.b[[this.effect]]))) +
    # Labels and theme
    labs(x = xlabel, y = ylabel) +
    scale_x_continuous(breaks = axis.breaks, labels = axis.labs) +
    theme_classic() +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none")
  
  plot_obj
}

# Function: plot_model_condeff_single
#
# Description:
# Creates a plot of conditional effects from a single Bayesian model, 
# including credible intervals and observed data points.
#
# Parameters:
# - model: A Bayesian model (from `brms`) to visualize.
# - this.effect: Predictor variable to plot.
# - this.resp: Response variable corresponding to the model.
# - point.data: Dataframe with observed data points for the model.
# - axis.breaks: Numeric vector for x-axis tick positions.
# - axis.labs: Labels for x-axis ticks.
# - xlabel, ylabel: Labels for x and y axes.
# - mod1color: Color for model plot elements. Default is 'navy'.
#
# Output:
# Returns a ggplot object with:
# - Conditional effect line for the model.
# - 95%, 80%, and 50% credible intervals as shaded ribbons.
# - Points for observed data.
plot_model_condeff_single <- function(model,
                                       this.effect,
                                       this.resp,
                                       point.data,
                                       axis.breaks,
                                       axis.labs,
                                       xlabel,
                                       ylabel,
                                       mod1color='navy'
){

  # Extract the data from conditional_effects
  cond_effects_data <- conditional_effects(model, effects = this.effect, resp = this.resp, plot = FALSE)
  plot_data <- cond_effects_data[[paste(this.resp,".",this.resp, "_",this.effect, sep='')]]
  
  ## fixing col names TODO: should fix the col names in prep
  if ('PD.obligate' %in% colnames(plot_data)){
    plot_data$PDobligate <- plot_data$PD.obligate
  }
  if ('PD.obligate.log' %in% colnames(plot_data)){
    plot_data$PDobligatelog <- plot_data$PD.obligate.log
  }
  if ('PD.transient' %in% colnames(plot_data)){
    plot_data$PDtransient <- plot_data$PD.transient
  }
  if ('PD.transient.log' %in% colnames(plot_data)){
    plot_data$PDtransientlog <- plot_data$PD.transient.log
  }
  
  # Plot using ggplot2 for credible intervals with geom_ribbon
  plot_obj <- ggplot(plot_data, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add ribbons for the 95%, 80%, and 50% credible intervals
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, 
                fill=mod1color,
                color = mod1color, linetype='dotted') +
    geom_ribbon(aes(ymin = lower__ + 0.1 * (upper__ - lower__),
                    ymax = upper__ - 0.1 * (upper__ - lower__)),
                alpha = 0.3, 
                fill=mod1color, 
                color = mod1color, linetype='dashed') +
    geom_ribbon(aes(ymin = lower__ + 0.25 * (upper__ - lower__),
                    ymax = upper__ - 0.25 * (upper__ - lower__)),
                alpha = 0.4, 
                fill=mod1color,
                color = mod1color, linetype='solid') +
    #Add line for the estimates
    geom_line(data = plot_data, color = 'black', linewidth=2.5, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    #Add line for the estimates
    geom_line(data = plot_data, color = mod1color, linewidth=2, aes(x = .data[[this.effect]], y = .data$estimate__)) +
    # Add points for original data
    geom_point(data = point.data, aes(x = .data[[this.effect]], y = .data[[this.resp]]),
               fill = mod1color, alpha = 0.9,color="black", pch=21, cex=3) +
    coord_cartesian(xlim = range(point.data[[this.effect]])) +
    # Labels and theme
    labs(x = xlabel, y = ylabel) +
    scale_x_continuous(breaks = axis.breaks, labels = axis.labs) +
    theme_classic() +
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          text = element_text(size=16),
          legend.position = "none")
  
  plot_obj
}

