#!/temp/saladi/R311/bin/Rscript

#
# Shyam Saladi (saladi@caltech.edu)
# Bil Clemons Laboratory (clemonslab.caltech.edu)
# Last Revised: April 17, 2016
#
# input: svmlight model dat file, .RData file, name of transformation, input stats file
# output: predictions
#

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
    test_data <- read.csv(file = args[4], header = TRUE, as.is = TRUE) #, row.names = 1)

    # if we want to capture a specific set of features defined by args[5]
    if (length(args) > 4) {
        feature_type <- read.csv("/ul/saladi/nycomps_fluman_raw/feature_type.csv", header = TRUE)
        featureselection <- feature_type[feature_type$type == args[5],]$feature
        columns <- data.frame(names = colnames(trans_data))
        columns <- columns[columns$names %in% as.character(featureselection),]
        trans_data <- subset(trans_data, select = as.character(columns))
    }

    scores <- prediction_fn(test_data, xtrans = get(args[3]), model_dat_fn = args[1])
    options(digits = 9)
    cat(scores, sep = "\n")
}

# either svm_classify or svm_rank_classify will work fine
# See: svmlight.joachims.org
# classify_cmd <- "/temp/saladi/svm_light/svm_classify"

prediction_fn <- function(test_data, xtrans, model_dat_fn,
    classify_cmd = "/temp/saladi/svm_rank/svm_rank_classify") {

    # Using caret_6.0-52
    # with R version 3.2.2 (2015-08-14), "Fire Safety" on x86_64-redhat-linux-gnu
    require(caret, quietly = TRUE)

    # set up files for svmlight
    data_fn <- tempfile("svmrank_", fileext = ".dat")
    out_fn <- paste0(data_fn, ".predictions")

    # In case there were features that were not calculated,
    # e.g. RNAss, add columns with NA for them
    # Contribution will not be taken into account for score calculation

    # features we use for the ML
    prediction_features <- names(xtrans$mean)

    # features in the test dataset
    test_data_features <- colnames(test_data)

    # features requested for the ML but not in the input dataset
    # keep in mind setdiff gives the asymmetric difference
    # (http://stat.ethz.ch/R-manual/R-patched/library/base/html/sets.html)
    missing_features <- setdiff(prediction_features, colnames(test_data))

    for (x in missing_features) {
        test_data[[x]] <- NA
    }
    rm(x)

    # keep only the columns we want for the ML (get rid of the extras)
    test_data <- test_data[, prediction_features]

    # do the preprocessing (scaling and centering)
    test_data <- predict(xtrans, test_data)

    # put the file in svmlight file format
    trans_data <- svmlight.file(as.matrix(test_data))

    # write the file in svmlight format
    write.table(data.frame(outcome = rep(1, nrow(test_data)),
                           vectorid = rep("qid:1", nrow(test_data)),
                           trans_data),
                file = data_fn,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    # remove features with no value (NA)
    system(paste("perl -p -i -e 's/ \\d+:NA//g'", data_fn))
    system(paste("perl -p -i -e 's/ \\d+:NaN//g'", data_fn))

    # run the classification
    temp <- system2(classify_cmd,
                    args = c(data_fn, model_dat_fn, out_fn),
                    stdout = NULL, stderr = "")

    preds <- read.delim(out_fn, header = FALSE)$V1

    # remove svmlight input and predictions file
    unlink(c(data_fn, out_fn))

    preds
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
        erg[ , -1] <- matrix(paste(j - train, t(x[,j]), sep = ":"), ncol = ncol(x) - 1, byrow = TRUE)
    }
    return(erg)
}

if (!interactive()) {
    args <- commandArgs(TRUE)
    main(args)
}
