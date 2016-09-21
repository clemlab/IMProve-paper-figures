# if on brmaster, need to copy over from homedir to data drive
# cp ~/homedir/ml-ecoli-paper/training_bootstraps_with_bca_ci_envinput.RData ~/compute_temp/
# cp ~/homedir/ml-ecoli-paper/training_bootstraps_with_bca_ci.R ~/compute_temp/

bs_count <- 10 ^ 6
cores <- 48

library(magrittr)
library(dplyr)
library(boot)
library(parallel)

setwd("~/compute_temp")

load("training_bootstraps_with_bca_ci_envinput.RData")
summary_fn <- function(...) {
    cor(... , method = "spearman", use = "pairwise.complete.obs")
}

# bootstrap for confidence intervals
cor_boot <- function(data, indices,  ...) {
    summary_fn(data[indices,]$ml21, data[indices,]$outcome, ...)
}

# multi
cor_boot_wrap <- function(data, file_name, FUN = cor_boot,
                          count = bs_count, mc.cores = cores, ...) {
    boot_data <- boot(data = data %>% select(ml21, outcome), statistic = FUN,
                      R = bs_count, parallel = "multicore", ncpus = cores)
    save(boot_data, file = paste0(file_name, ".RData"))
    boot.ci(boot_data, type = "bca")
}

ecoli_xy_boot_in_ci <-
    cor_boot_wrap(filter(ecoli_training, cterm == "in"), "ecoli_xy_boot_in")
gc()
ecoli_xy_boot_out_ci <-
    cor_boot_wrap(filter(ecoli_training, cterm == "out"), "ecoli_xy_boot_out")
gc()
ecoli_xy_boot_fluman_ci <-
    cor_boot_wrap(filter(ecoli_training, cterm == "fluman"), "ecoli_xy_boot_fluman")
gc()
ecoli_daley_fluman_boot_ci <-
    cor_boot_wrap(rename(ecoli_daley_fluman, ml21 = fluman_avg, outcome = daley_avg),
                  "ecoli_daley_fluman_boot")
gc()

# p-value by shuffling


# single
shuff_cor_pval <- function(index, data, x_col="ml21", y_col="outcome", ...) {
     y_shuff <- base::sample(data[, y_col], size = nrow(data), replace = FALSE)
     summary_fn(data[, x_col], y_shuff, ...)
}

# summary statistic
calc_pval <- function(orig_df, shuf, x_col="ml21", y_col="outcome", ...) {
    stat <- summary_fn(orig_df[, x_col], orig_df[, y_col], ...)
    
    # adjustment to account for the original data that does meet the criteria
    (sum(shuf >= stat) + 1) / (length(shuf) + 1)
}

# multi
shuff_cor_pval_wrap <- function(data, file_name, FUN = shuff_cor_pval,
                                count = bs_count, mc.cores = cores, ...) {
    shuf_data <- mclapply(seq(bs_count), FUN,
                          data = data %>% select(ml21, outcome),...,
                          mc.cores = cores) %>% unlist
    save(shuf_data, file = paste0(file_name, ".RData"))
    calc_pval(data, shuf_data)
}

ecoli_xy_shuf_in_pval <- filter(ecoli_training, cterm == "in") %>%
    shuff_cor_pval_wrap(file_name = "ecoli_xy_shuf_in")
gc()

ecoli_xy_shuf_out_pval <- filter(ecoli_training, cterm == "out") %>%
    shuff_cor_pval_wrap(file_name = "ecoli_xy_shuf_out")
gc()

ecoli_xy_shuf_fluman_pval <- filter(ecoli_training, cterm == "fluman") %>%
    shuff_cor_pval_wrap(file_name = "ecoli_xy_shuf_fluman")
gc()

ecoli_daley_fluman_shuf_pval <- rename(ecoli_daley_fluman,ml21 = fluman_avg,
                                       outcome = daley_avg) %>%
    shuff_cor_pval_wrap(file_name = "ecoli_daley_fluman_shuf")
gc()

save(ecoli_xy_boot_in_ci, ecoli_xy_boot_out_ci, ecoli_xy_boot_fluman_ci,
     ecoli_daley_fluman_boot_ci,
     ecoli_xy_shuf_in_pval, ecoli_xy_shuf_out_pval, ecoli_xy_shuf_fluman_pval,
     ecoli_daley_fluman_shuf_pval,
     file = "training_bootstraps_with_bca_ci_summary.RData")

