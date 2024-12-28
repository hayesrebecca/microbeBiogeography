
# Function: chao
#
# Description:
# Calculates species richness estimates and related statistics based on the 
# provided data. Supports abundance and incidence-based data types, including 
# transformations for frequency count formats. Provides details about rare species 
# groups and species tables.
#
# Parameters:
# - data: Input dataset for species richness estimation. Can be:
#   - Numeric vector (abundance data).
#   - Matrix or data frame (incidence raw data with rows as species and columns as sampling units).
# - datatype: Type of input data. Options:
#   - "abundance" (default): Abundance-based data.
#   - "incidence_raw": Raw incidence data (presence-absence matrix).
#   - "abundance_freq_count": Frequency count format for abundance data.
#   - "incidence_freq_count": Frequency count format for incidence data.
# - k: Cutoff value to define "rare" species. Non-negative integer (default: 10).
# - conf: Confidence level for credible intervals. Numeric between 0 and 1 (default: 0.95).
#
# Output:
# Returns a list of class "ChaoSpecies" with:
# - Basic_data_information: Summary statistics about the input data.
# - Rare_species_group / Infreq_species_group: Details about rare or infrequent species.
# - Species_table: Richness estimates, standard errors, and confidence intervals.
#
chao <- function (data, datatype = c("abundance"),
                  k = 10, conf = 0.95){
    if (is.matrix(data) == T || is.data.frame(data) == T) {
        if (datatype != "incidence_raw") {
            if (ncol(data) == 1) {
                data <- data[, 1]
            }
            else {
                data <- data[1, ]
            }
        }
        else {
            t <- ncol(data)
            dat <- rowSums(data)
            dat <- as.integer(dat)
            t_infreq <- sum(colSums(data[which(dat < k), ]) >=
                1)
            data <- dat
            data <- c(t_infreq, t, data)
        }
    }
    if (datatype == "abundance_freq_count") {
        data <- as.integer(data)
        length4b <- length(data)
        data <- rep(data[seq(1, length4b, 2)], data[seq(2, length4b,
            2)])
        names(data) <- paste("x", 1:length(data), sep = "")
        datatype <- "abundance"
    }
    if (datatype == "incidence_freq_count") {
        t <- as.integer(data[1])
        data <- data[-c(1)]
        data <- as.integer(data)
        lengthdat <- length(data)
        data <- rep(data[seq(1, lengthdat, 2)], data[seq(2, lengthdat,
            2)])
        data <- c(t, data)
        names(data) <- c("T", paste("y", 1:(length(data) - 1),
            sep = ""))
        datatype <- "incidence_freq"
    }
    method <- "all"
    if (k != round(k) || k < 0)
        stop("Error: The cutoff t to define less abundant species must be non-negative integer!")
    if (is.numeric(conf) == FALSE || conf > 1 || conf < 0)
        stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
    if (datatype == "abundance") {
        f <- function(i, data) {
            length(data[which(data == i)])
        }
        if (f(1, data) == sum(data)) {
            stop("Error: The information of data is not enough.")
        }
        z <- (list(Basic_data_information = basicAbun(data, k)[[1]],
            Rare_species_group = RareSpeciesGroup(data, k), Species_table = round(SpecAbunOut(data,
                method, k, conf), 3)))
    }
    else if (datatype == "incidence_raw") {
        dat <- data[-1]
        Q <- function(i, data) {
            length(data[which(data == i)])
        }
        if (Q(1, dat) == sum(dat)) {
            stop("Error: The information of data is not enough.")
        }
        z <- (list(Basic_data_information = basicInci(data[-1],
            k)[[1]], Infreq_species_group = InfreqSpeciesGroup(data[-1],
            k), Species_table = round(SpecInciOut_raw(data, method,
            k, conf), 3)))
    }
    else if (datatype == "incidence_freq") {
        dat <- data[-1]
        Q <- function(i, data) {
            length(data[which(data == i)])
        }
        if (Q(1, dat) == sum(dat)) {
            stop("Error: The information of data is not enough.")
        }
        z <- (list(Basic_data_information = basicInci(data, k)[[1]],
            Infreq_species_group = InfreqSpeciesGroup(data, k),
            Species_table = round(SpecInciOut(data, method, k,
                conf), 3)))
    }
    else {
        stop("Error: The data type is wrong.")
    }
    class(z) <- c("ChaoSpecies")
    z
}
