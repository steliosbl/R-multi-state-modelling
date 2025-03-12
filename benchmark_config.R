library(tidyverse)
library(glue)

# Test 1 - Fixed Testing Static Training #######################################
test1 <- function(out_name = "bench_test1.tsv") {
    param_grid <- expand.grid(
        test = "1",
        model = c("survival", "mstate", "flexsurv"),
        elements = 2:4,
        n_covs = 4,
        n_train = 2^(11:16),
        n_test = 16,
        times = -1,
        cpus = 1,
        rep = 1:3
    ) %>%
        arrange(n_train) %>%
        mutate(row = row_number()) %>%
        select(
            row, test, model, elements, n_covs,
            n_train, n_test, times, cpus, rep
        )
}

# Test 2 - Fit Scaling  ########################################################
test2 <- function(out_name = "bench_test2.tsv") {
    param_grid <- expand.grid(
        test = "2",
        model = c("survival", "mstate", "flexsurv"),
        elements = 2:4,
        n_covs = c(2, 4, 8, 16),
        n_train = 2^(11:21),
        n_test = 0,
        times = -1,
        cpus = 1,
        rep = 1:3
    ) %>%
        arrange(n_train) %>%
        mutate(row = row_number()) %>%
        select(
            row, test, model, elements, n_covs,
            n_train, n_test, times, cpus, rep
        )
}

# Test 3 - Prediction Parallelism ##############################################
test3 <- function() {
    param_grid <- expand.grid(
        test = "3",
        model = c("survival", "mstate", "flexsurv"),
        elements = 2:3,
        n_covs = 4,
        n_train = 2^13,
        n_test = c(2^4, 2^5, 2^6, 2^7, 2^8, 2^9, 2^10),
        times = -1,
        cpus = c(2^0, 2^1, 2^2, 2^3, 2^4, 2^5, 2^6),
        rep = 1:3
    ) %>%
        arrange(cpus) %>%
        mutate(row = row_number()) %>%
        select(
            row, test, model, elements, n_covs,
            n_train, n_test, times, cpus, rep
        )
}

# Test 4 - Coarsening ##########################################################
test4 <- function() {
    param_grid <- expand.grid(
        test = "4",
        model = c("survival", "mstate", "flexsurv"),
        elements = 3,
        n_covs = 4,
        n_train = 2^14,
        n_test = c(2^4, 2^6),
        times = 2^(4:19),
        cpus = 1,
        rep = 1:3
    ) %>%
        arrange(times) %>%
        mutate(row = row_number()) %>%
        select(
            row, test, model, elements, n_covs,
            n_train, n_test, times, cpus, rep
        )
}

export_tsv <- function(param_grid, out_name, out_dir = "data/config") {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    out_path <- file.path(out_dir, out_name)
    write.table(
        param_grid, out_path,
        sep = "\t",
        row.names = FALSE, quote = FALSE
    )
}

export_tsv(test1(), "test1.tsv")
export_tsv(test2(), "test2.tsv")
export_tsv(test3(), "test3.tsv")
export_tsv(test4(), "test4.tsv")
