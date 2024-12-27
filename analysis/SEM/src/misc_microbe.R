#' Miscellaneous Helper Functions for Data Manipulation and Cleaning
#'
#' This module contains a collection of utility functions designed for common data 
#' manipulation, cleaning, and processing tasks in R. These functions support 
#' generating visual outputs, handling missing values, standardizing data, modifying 
#' colors, and cleaning strings. They aim to streamline workflows and provide 
#' reusable solutions for repetitive tasks.
#'
#' Functions Included:
#' - `pdf.f`: Save the output of a plotting function to a PDF file.
#' - `add.alpha`: Add transparency (alpha) to colors for use in visualizations.
#' - `standardize`: Standardize a numeric vector by subtracting the mean and dividing 
#'   by the standard deviation.
#' - `unstandardize`: Reverse standardization by restoring the original scale of a numeric vector.
#' - `standardize.axis`: Standardize a vector relative to a reference vector for 
#'   consistent scaling in visualizations.
#' - `fix.white.space`: Remove extra spaces and clean up character strings.
#' - `check_for_NA`: Check for missing values (NAs) in specific elements of a global 
#'   data frame and report counts.
#'
#' Usage:
#' These functions are intended to be modular and reusable across different projects.
#' Each function includes detailed inline documentation with examples to help users 
#' understand their behavior and applications.
#'
#' Note:
#' - The global object `spec.net` is required for the `check_for_NA` function to operate.
#' - Ensure appropriate arguments are provided for each function to prevent errors.
#'



# pdf.f: Save the output of a function to a PDF file
# 
# This function wraps around `pdf()` to create a PDF file from a user-provided 
# function `f`, ensuring the graphics device is properly closed.
#
# Arguments:
#   f: A function that generates the content to save in the PDF. 
#      It should take no arguments.
#   file: A string specifying the path to the PDF file.
#   ...: Additional arguments passed to `pdf()`.
#
# Example:
#   pdf.f(function() {
#     plot(1:10, 1:10, main = "Example Plot")
#   }, file = "example.pdf")
#
pdf.f <- function(f, file, ...) {
    cat(sprintf('Writing %s\n', file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
}

# add.alpha: Add transparency to colors
#
# This function modifies a vector of colors by adding an alpha (transparency) value.
#
# Arguments:
#   col: A vector of colors, specified as names (e.g., "red") or hex codes (e.g., "#FF0000").
#   alpha: A numeric value between 0 and 1 specifying the transparency level, 
#          where 0 is fully transparent and 1 is fully opaque. Default is 0.2.
#
# Example:
#   add.alpha(c("red", "blue", "green"), alpha = 0.5)
#   # Returns a vector of semi-transparent versions of the input colors.
#
add.alpha <- function(col, alpha=0.2){
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
              rgb(x[1], x[2], x[3],
                  alpha=alpha))
}

# standardize: Standardize a numeric vector
#
# This function standardizes a numeric vector by subtracting the mean and dividing 
# by the standard deviation, resulting in a vector with a mean of 0 and a standard 
# deviation of 1.
#
# Arguments:
#   x: A numeric vector to standardize.
#
# Details:
#   Missing values (NA) are ignored in the calculation of the mean and standard deviation.
#
# Example:
#   standardize(c(1, 2, 3, 4, 5))
#   # Returns: c(-1.264911, -0.632456, 0, 0.632456, 1.264911)
#
standardize <- function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# unstandardize: Reverse standardization of a numeric vector
#
# This function reverses the standardization of a numeric vector by scaling it 
# back to the original mean and standard deviation.
#
# Arguments:
#   x: A standardized numeric vector to unstandardize.
#   orig: The original numeric vector from which the mean and standard deviation 
#         will be computed to reverse the standardization.
#
# Details:
#   The function uses the mean and standard deviation of the `orig` vector to 
#   transform `x` back to the original scale. Missing values (NA) in `orig` are 
#   ignored during the calculation.
#
# Example:
#   standardized <- standardize(c(1, 2, 3, 4, 5))
#   unstandardize(standardized, c(1, 2, 3, 4, 5))
#   # Returns: c(1, 2, 3, 4, 5)
#
unstandardize <- function(x, orig) {
  (x * sd(orig, na.rm = TRUE)) + mean(orig, na.rm = TRUE)
}

# standardize.axis: Standardize a numeric vector relative to a reference vector for plotting
#
# This function standardizes a numeric vector `x` using the mean and standard deviation 
# of a reference vector `orig`. This is useful for ensuring that `x` is standardized 
# relative to an external dataset when plotting scaled results.
#
# Arguments:
#   x: A numeric vector to standardize.
#   orig: A reference numeric vector used to calculate the mean and standard deviation.
#
# Details:
#   The function subtracts the mean of `orig` from `x` and divides by the standard deviation 
#   of `orig`. Missing values (NA) in `orig` are ignored during the calculations.
#
# Example:
#   standardize.axis(c(10, 12, 14), c(1, 2, 3, 4, 5))
#   # Returns: c(15.81139, 19.76393, 23.71648)
#
standardize.axis <- function(x, orig) {
  (x - mean(orig, na.rm = TRUE)) / sd(orig, na.rm = TRUE)
}


# fix.white.space: Remove extra spaces from character strings
#
# This function cleans up character strings by removing extra spaces, reducing 
# all consecutive spaces to a single space, and ensuring no leading spaces remain.
#
# Arguments:
#   d: A vector of character strings to process.
#
# Details:
#   - Consecutive spaces (e.g., 2 or more spaces in a row) are replaced with a single space.
#   - Leading spaces at the beginning of the string are removed.
#   - Handles both character and non-character inputs by converting them to character.
#
# Example:
#   fix.white.space(c("  Hello   world  ", "This    is  R!"))
#   # Returns: c("Hello world", "This is R!")
#
fix.white.space <- function(d) {
  d <- as.character(d)
  remove.first <- function(s) substr(s, 2, nchar(s))
  d <- gsub("      ", " ", d, fixed=TRUE)
  d <- gsub("     ", " ", d, fixed=TRUE)
  d <- gsub("    ", " ", d, fixed=TRUE)
  d <- gsub("   ", " ", d, fixed=TRUE)
  d <- gsub("  ", " ", d, fixed=TRUE)

  tmp <- strsplit(as.character(d), " ")
  d <- sapply(tmp, function(x) paste(x, collapse=" "))

  first <- substr(d, 1, 1)
  d[first==" "] <- remove.first(d[first==" "])
  d
}


# check_for_NA: Check for missing values (NA) in specific elements of a list
#
# This function checks for missing values (NAs) in specific elements of a global 
# data frame `spec.net` and prints the count of NAs for each element in `which.col`.
#
# Arguments:
#   which.col: A character vector specifying the names of the elements in `spec.net` 
#               to check for missing values.
#
# Details:
#   - For each name in `which.col`, the function calculates the number of NAs in the 
#     corresponding element of `spec.net`.
#   - If NAs are present, the function prints a message specifying the name of the 
#     element and the count of NAs.
#   - The global object `spec.net` must be defined and contain the specified elements.
#
# Example:
#   spec.net <- list(
#     col1 = c(1, 2, NA, 4),
#     col2 = c(5, 6, 7),
#     col3 = c(NA, NA, 3)
#   )
#   check_for_NA(c("col1", "col3"))
#   # Output:
#   # [1] "microbe col1 NAs: 1"
#   # [1] "microbe col3 NAs: 2"
#
check_for_NA <- function(which.col){
  for (x in which.col){
    anyNA <- length(spec.net[[x]][is.na(spec.net[[x]])])
    if (anyNA > 0){
    print(paste('microbe', x, "NAs:", length(spec.net[[x]][is.na(spec.net[[x]])])))
    }
  }
}

