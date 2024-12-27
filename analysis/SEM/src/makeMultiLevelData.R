# Module: Data Preparation for Multi-Level Analysis
#
# This module provides functions for preparing individual-level data for multi-level 
# analyses. It includes functionality for splitting data by year, assigning unique 
# site identifiers, and creating weight columns used in statistical modeling.
#
# Functions:
# - makeDataMultiLevel: Processes individual-level data into a format suitable for 
#   multi-level modeling, including site-level weights and unique site identifiers.
# - addWeightCol: Adds weight and site identifier columns to a dataset for a given year.
#

# makeDataMultiLevel: Prepare individual-level data for multi-level analysis
#
# This function splits individual-level data by year, assigns unique site identifiers, 
# and creates weight columns for each year. The processed data is then combined into a 
# single dataset for multi-level analysis.
#
# Arguments:
#   indiv.data: A data frame containing individual-level data.
#   site.col: The name of the column specifying site identifiers (as a string).
#   year.col: The name of the column specifying years (default is "Year").
#   weight.col.name: The name of the column to store weight values (default is "Weights").
#   site.id.col.name: The name of the column to store unique site IDs (default is "SiteIDs").
#
# Returns:
#   A data frame with added weight and site identifier columns, prepared for multi-level analysis.
#
# Example:
#   indiv.data <- data.frame(Site = rep(c("A", "B"), each = 3), Year = c(2020, 2021, 2020, 2021, 2020, 2021))
#   makeDataMultiLevel(indiv.data, site.col = "Site")
#
makeDataMultiLevel <- function(indiv.data, site.col, year.col="Year", weight.col.name="Weights",
                               site.id.col.name="SiteIDs"){
    ## split data by year
    indiv.data.split <- split(indiv.data, indiv.data[, year.col])
    
    ## maybe in the future this will need to be an sapply
    out.indiv.data <- lapply(indiv.data.split, addWeightCol,
                             site.col=site.col,  weight.col.name= weight.col.name,
                             site.id.col.name=site.id.col.name)
    
    out.indiv.data <- do.call(rbind, out.indiv.data)
    
    return(out.indiv.data)
    
}

# addWeightCol: Add weight and site ID columns to data for a single year
#
# This function adds two new columns to the input dataset: a column for unique site 
# identifiers and a weight column for multi-level modeling. The weight column is set 
# to 0 for rows beyond the first occurrence of a site within the year.
#
# Arguments:
#   each.year.dat: A data frame containing data for a single year.
#   site.col: The name of the column specifying site identifiers (as a string).
#   weight.col.name: The name of the column to store weight values (default is "Weights").
#   site.id.col.name: The name of the column to store unique site IDs (default is "SiteIDs").
#
# Returns:
#   A data frame with added weight and site identifier columns.
#
# Example:
#   each.year.dat <- data.frame(Site = c("A", "A", "B", "B"))
#   addWeightCol(each.year.dat, site.col = "Site")
#
addWeightCol <- function(each.year.dat, site.col, 
                         weight.col.name="Weights",
                         site.id.col.name="SiteIDs"){
    site.ids <- unlist(tapply(each.year.dat[, site.col],
                              each.year.dat[, site.col],
                              function(x) 1:length(x)))
    
    
    names(site.ids) <- NULL
    each.year.dat[, site.id.col.name] <- site.ids
    each.year.dat[, weight.col.name] <- each.year.dat[, site.id.col.name]
    each.year.dat[, weight.col.name][each.year.dat[weight.col.name] > 1] <- 0
    return(each.year.dat)
}
