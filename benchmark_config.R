library(tidyverse)
library(glue)

# Test 1 - Fixed Testing Static Training #######################################
test1 <- function(out_name = "bench_test1.tsv") {
    param_grid <- expand.grid(
        test = "1",
        model = c("survival", "mstate"),
        n_elements = 2:4,
        n_covariates = 4,
        n_train = c(2^8, 2^9, 2^10, 2^11),
        n_test = 32,
        n_cpu = 1,
        rep = 1:3
    ) %>%
        arrange(n_train) %>%
        mutate(row = row_number()) %>%
        select(row, test, model, n_elements, n_covariates, n_train, n_test, n_cpu, rep)
}

# Test 2 - Fit Scaling  ########################################################
test2 <- function(out_name = "bench_test2.tsv") {
    param_grid <- expand.grid(
        test = "2",
        model = c("survival", "mstate"),
        n_elements = 2:4,
        n_covariates = c(2, 4, 8, 16),
        n_train = c(2^8, 2^9, 2^10, 2^11, 2^12, 2^13, 2^14, 2^15, 2^16, 2^17, 2^18, 2^19),
        n_test = 0,
        n_cpu = 1,
        rep = 1:3
    ) %>%
        arrange(n_train) %>%
        mutate(row = row_number()) %>%
        select(row, test, model, n_elements, n_covariates, n_train, n_test, n_cpu, rep)
}

# Test 3 - Prediction Parallelism ##############################################
test3 <- function() {
    n_test <- c(2^4, 2^5, 2^6, 2^7, 2^8, 2^9, 2^10)
    n_cpu <- c(2^0, 2^1, 2^2, 2^3, 2^4, 2^5, 2^6)

    fixed_params <- list(
        test = "3",
        model = c("survival", "mstate"),
        n_elements = 2:3,
        n_covariates = 4,
        n_train = 2^13,
        rep = 1:3
    )

    param_grids <- list()
    for (i in seq_along(n_test)) {
        param_grids[[i]] <- do.call(
            expand.grid,
            c(
                fixed_params,
                n_test = n_test[i],
                n_cpu = n_cpu[i]
            )
        ) %>%
            arrange(n_train) %>%
            mutate(row = row_number()) %>%
            select(row, test, model, n_elements, n_covariates, n_train, n_test, n_cpu, rep)
    }
    param_grids
}

export_tsv <- function(param_grid, out_name, out_dir = "data/config") {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
    }
    out_path <- file.path(out_dir, out_name)
    write.table(param_grid, out_path, sep = "\t", row.names = FALSE, quote=FALSE)
}

#export_tsv(test1(), "bench_test1.tsv")
export_tsv(test2(), "bench_test2.tsv")
parallel_param_grids <- test3()
lapply(seq_along(parallel_param_grids), function(i) {
    export_tsv(parallel_param_grids[[i]], glue("bench_test3_{2^(i-1)}.tsv"))
})
