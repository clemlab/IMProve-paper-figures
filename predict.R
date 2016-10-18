#!/usr/bin/env Rscript

#
# Shyam Saladi (saladi@caltech.edu)
# Bil Clemons Laboratory (clemonslab.caltech.edu)
#
# input: svmlight model dat file, .RData file, name of transformation, input stats file
# output: predictions
#

library(tidyverse)
library(magrittr)

# install("myutils)
library(myutils)


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

if (!interactive()) {
    args <- commandArgs(TRUE)
    main(args)
}
