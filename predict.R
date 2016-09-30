#!/usr/bin/env Rscript

#
# Shyam Saladi (saladi@caltech.edu)
# Bil Clemons Laboratory (clemonslab.caltech.edu)
#
# input: svmlight model dat file, .RData file, name of transformation, input stats file
# output: predictions
#

library(readr)
library(dplyr)
library(magrittr)
library(caret)


main <- function(args) {
    if (is.na(args[4])) {
        cat("Need to supply a the following arguments:
            1 - SVMlight model .dat file
            2 - .RData file
            3 - name of transformation variable (see preProcessing, R caret), e.g. *xtrans
            4 - input .allstats.csv file\n", file = stderr())
        quit()
    }
    
    # load the workspace that has the transformation data (for preprocessing)
    load(args[2])
    
    # read in csvfile with the calculated features
    test_data <- read.csv(file = args[4], header = TRUE, as.is = TRUE)
    
    xtrans <- get(args[3])
    weights <- read_svmlight_model(args[1], names(xtrans$mean))
    scores <- prediction_fn(test_data, xtrans, weights)
    
    options(digits = 9)
    cat(scores, sep = "\n")
}

prediction_fn <- function(test_data, xtrans, weights) {
    require(data.table, quietly = TRUE)
    # In case there were features that were not calculated,
    # e.g. RNAss, add columns with NA for them
    # Contribution will not be taken into account for score calculation
    
    # features we use for the ML
    prediction_features <- names(weights)
    
    # features requested for the ML but not in the input dataset
    # keep in mind setdiff gives the asymmetric difference
    # (http://stat.ethz.ch/R-manual/R-patched/library/base/html/sets.html)
    missing_features <- setdiff(prediction_features, colnames(test_data))
    
    for (x in missing_features) {
        test_data[[x]] <- NA
    }
    
    # keep only the columns we want for the ML (get rid of the extras)
    test_data <- test_data[, prediction_features]
    
    # do the preprocessing (scaling and centering)
    if (!is.null(xtrans)) {
        library(caret)
        test_data <- predict(xtrans, test_data)
    }
    
    # remove features without a value
    test_data[is.na(test_data)] <- 0
    
    # reorder vector
    weights <- weights[colnames(test_data)]
    for (i in seq_along(test_data))
        set(test_data, j = i, value = test_data[[i]] * weights[[i]])
    
    # prediction
    test_data %>% rowSums(na.rm = TRUE)
}

read_svmlight_model <- function(filename, feat_names) {
    weights <- read_delim(filename, skip = 11, delim = " ", col_names = FALSE,
                          col_types = cols(.default = col_character(),
                                           X1 = col_integer())) %>%
        select(-X1) %>%
        gather(id, weight) %>%
        filter(weight != "#") %>%
        separate(weight, into = c("idx", "weight"), sep = ":") %>%
        mutate(weight = as.numeric(weight)) %>%
        select(-id, -idx) %>% unlist
    names(weights) <- feat_names
    weights
}


# Adapted from klaR
# http://www.inside-r.org/packages/cran/klaR/docs/svmlight
svmlight.file <- function(x, train = FALSE, ...)
{
    if (is.vector(x)) x <- t(x)
    erg <- x
    sn <- 1:nrow(x)
    if (!train) erg[sn, 1] <- paste("1:", x[sn, 1], sep = "")
    if (ncol(x) > 1) {
        j <- 2:ncol(x)
        erg[ , -1] <- matrix(paste(j - train, t(x[,j]), sep = ":"),
                             ncol = ncol(x) - 1, byrow = TRUE)
    }
    return(erg)
}

if (!interactive()) {
    args <- commandArgs(TRUE)
    main(args)
}
