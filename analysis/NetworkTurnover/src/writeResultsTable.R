#' Write Model Summary and Posterior Samples to Tables
#'
#' This function extracts model coefficients and their posterior distributions from a fitted model object 
#' (e.g., from `brms` or `rstanarm`), computes the probability of being greater than and less than zero for each coefficient, 
#' and saves the results in both a formatted LaTeX table (for use in manuscripts) and a CSV file.
#'
#' Arguments:
#'   mod.output: A fitted model object (e.g., from `brms` or `rstanarm`) containing the model output, including 
#'               coefficient estimates and posterior samples.
#'   mod.name: A character string representing the name of the model, used to generate file names for the saved 
#'             tables. This name will be used to create both a LaTeX-formatted `.txt` table and a `.csv` file.
#'
#' Details:
#'   - The function first extracts the summary of fixed effect coefficients and formats them into a data frame.
#'   - Posterior samples for each coefficient are then retrieved, and the coefficients are filtered to match those 
#'     in the summary table.
#'   - The function computes the probabilities that each coefficient is greater than zero (`Pgt0`) and less than zero 
#'     (`Plt0`), and adds these as new columns in the summary table.
#'   - Two tables are saved:
#'     - A LaTeX-formatted table (`.txt`) for inclusion in LaTeX documents, with coefficients separated by an ampersand (`&`) for table formatting.
#'     - A CSV file (`.csv`) with the same summary table for easier manipulation or export.
#'
#' Returns:
#'   This function does not return any value. It writes two files:
#'     - A LaTeX-formatted table saved in the directory `"saved/tables/"` with a `.txt` extension.
#'     - A CSV file saved in the same directory with a `.csv` extension.
#'
#' Example:
#'   # Assuming 'mod' is a fitted model object (e.g., from `brms` or `rstanarm`):
#'   write.ms.table(mod.output = mod, mod.name = "my_model_results")
#'
#' See Also:
#'   `posterior_samples` (for extracting posterior samples)
write.ms.table <- function(mod.output, mod.name){
    sum.mod <- as.data.frame(round(summary(mod.output)$fixed,2))

    coeffs <- c(paste0("b_",
                       rownames(sum.mod)),
                paste0("bs_",
                       rownames(sum.mod)))

    samps.mod <- posterior_samples(mod.output)

    coeffs <- coeffs[coeffs %in% colnames(samps.mod)]

    samps.mod <- samps.mod[, coeffs]

    coeff.samps <- colnames(samps.mod)
    coeff.samps.sum <- sub("[a-z]*_", "", coeff.samps)

   samps.mod <- samps.mod[order(match( coeff.samps.sum, rownames(sum.mod)))]

    sum.mod$Pgt0  <- round(apply(samps.mod, 2, function(x)
        sum(x > 0)/length(x)), 2)

    sum.mod$Plt0  <- round(apply(samps.mod, 2, function(x)
        sum(x < 0)/length(x)),2)


    write.table(sum.mod,
                file=sprintf("saved/tables/%s.txt", mod.name),
                sep="&")

    write.csv(sum.mod,
              file=sprintf("saved/tables/%s.csv", mod.name))
}

# Function: write.summ.table
#
# Description:
# This function generates a summary table from the output of a Bayesian model. It extracts the fixed effect coefficients and their 
# corresponding posterior samples, calculates the proportion of posterior samples greater than zero (Pgt0) and less than zero (Plt0),
# and returns a data frame containing the fixed effects' estimates, standard errors, credible intervals, and the aforementioned probabilities.
# The table is sorted to match the order of the coefficients in the model's summary.
#
# Parameters:
# - mod.output: A fitted Bayesian model object (e.g., from `brms` or `rstanarm`) that contains the model's fixed effects and posterior samples.
#   The model must include posterior samples of the fixed effects.
#
# Output:
# - Returns a data frame with the following columns:
#   - Estimate: The mean of the posterior distribution for each fixed effect.
#   - Est.Error: The standard error of the posterior distribution for each fixed effect.
#   - l-95% CI: The lower bound of the 95% credible interval for each fixed effect.
#   - u-95% CI: The upper bound of the 95% credible interval for each fixed effect.
#   - Pgt0: The proportion of posterior samples for each coefficient that are greater than 0.
#   - Plt0: The proportion of posterior samples for each coefficient that are less than 0.
#
# Notes:
# - The function extracts the fixed effect summary from the model output and then calculates the posterior samples for these coefficients.
# - It computes the probabilities `Pgt0` and `Plt0` to show the likelihood of each coefficient being greater than or less than zero, respectively.
# - The function assumes that the model output contains posterior samples for the fixed effects.
# - All proportions (`Pgt0` and `Plt0`) are rounded to two decimal places.
#
# Example Usage:
# write.summ.table(mod_output)
#
# This will generate a summary table for the fitted model `mod_output`, including the fixed effect estimates, their credible intervals,
# and the proportions of posterior samples greater than and less than zero.
#
# The function can be used to evaluate the significance of the fixed effects based on their posterior distributions.
#
write.summ.table <- function(mod.output){
    sum.mod <- as.data.frame(round(summary(mod.output)$fixed,4))
    
    coeffs <- c(paste0("b_",
                       rownames(sum.mod)),
                paste0("bs_",
                       rownames(sum.mod)))
    
    samps.mod <- posterior_samples(mod.output)
    
    coeffs <- coeffs[coeffs %in% colnames(samps.mod)]
    
    samps.mod <- samps.mod[, coeffs]
    
    coeff.samps <- colnames(samps.mod)
    coeff.samps.sum <- sub("[a-z]*_", "", coeff.samps)
    
    samps.mod <- samps.mod[order(match( coeff.samps.sum, rownames(sum.mod)))]
    
    sum.mod$Pgt0  <- round(apply(samps.mod, 2, function(x)
        sum(x > 0)/length(x)), 2)
    
    sum.mod$Plt0  <- round(apply(samps.mod, 2, function(x)
        sum(x < 0)/length(x)),2)
    
    
    return(sum.mod)
}
