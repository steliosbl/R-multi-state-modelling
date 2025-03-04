# HEADER #######################################################################
#
# Author: Stelios Boulitsakis Logothetis
# Email: sb2690@medschl.cam.ac.uk
#
# Date: 2025-02-27
# Name:	benchmark.R
# Description: This script benchmarks the performance of different
# multi-state model implementations on simulated datasets.
#
# SETUP ########################################################################
# Load required libraries
library(tidyverse)
library(survival)
library(mstate)
library(flexsurv)

source("mythoslib/convert.R")
source("mythoslib/hasse.R")


# Set seed for reproducibility
set.seed(123)

# FUNCTIONS ####################################################################
prep_hasse_counting <- function(walks, covariates, tmat) {
    result <- walks %>%
        rename(event = vertex) %>%
        augment_timeline_censoring("id", tmat, 1.0, "censor") %>%
        timeline_to_counting(., "id") %>%
        arrange(id, tstart) %>%
        mutate(from = factor(from), to = factor(to)) %>%
        mutate(to = fct_relevel(as.factor(to), "censor")) %>%
        left_join(covariates, by = "id")

    attr(result, "covariates") <- colnames(covariates)[
        colnames(covariates) != "id"
    ]

    result
}

prep_hasse_mstate <- function(walks, covariates, tmat) {
    cov_names <- colnames(covariates)[colnames(covariates) != "id"]

    df <- walks %>%
        rename(event = vertex) %>%
        augment_timeline_censoring("id", tmat, 1.0, "censor") %>%
        timeline_to_mstate("id", tmat, "censor") %>%
        select(-time_censor, -status_censor) %>%
        left_join(covariates, by = "id")

    msdata <- msprep(
        data = df,
        time = paste0("time_", rownames(tmat)),
        status = paste0("status_", rownames(tmat)),
        trans = tmat,
        id = "id",
        start = list(
            time = df$itime,
            state = match(df$istate, dimnames(tmat)[[1]])
        ),
        keep = cov_names
    )

    covariates_expanded <- msdata %>%
        expand.covs(cov_names, append = FALSE, longnames = TRUE)

    result <- bind_cols(msdata, covariates_expanded)
    attr(result, "covariates") <- cov_names
    attr(result, "covariates_expanded") <- colnames(covariates_expanded)

    result
}

fit_survival <- function(data) {
    covariates <- attr(data, "covariates")

    formula <- as.formula(
        paste0("Surv(tstart, tstop, to) ~ ", paste(covariates, collapse = "+"))
    )

    time_start <- Sys.time()
    model <- coxph(formula, data, istate = from, id = id)
    time_stop <- Sys.time()
    list(
        model = model,
        time = as.numeric(time_stop - time_start)
    )
}

predict_survival <- function(model, newdata, tmat) {
    newdata <- newdata %>%
        filter(tstart == 0)

    result <- by(newdata, newdata$from, function(df) {
        p0 <- as.integer(model$states == df[1, ]$from)
        df <- df %>% select(id, attr(newdata, "covariates"))

        time_start <- Sys.time()
        pred <- survfit(model, newdata = df, p0 = p0)
        time_stop <- Sys.time()

        list(
            pred = pred,
            time = time_stop - time_start
        )
    })

    list(
        time = sum(unlist(map(result, ~ as.numeric(.x$time)))),
        pred = lapply(result, function(x) x$pred)
    )
}

predict_survival_parallel <- function(model, newdata, tmat, n_cpu) {
    # Get unique IDs and partition them into N groups
    ids <- unique(newdata$id)
    partitions <- split(ids, cut(seq_along(ids), n_cpu, labels = FALSE))

    # Apply function in parallel
    time_start <- Sys.time()
    results <- parallel::mclapply(partitions, function(id_subset) {
        subset <- newdata %>% filter(id %in% id_subset)
        predict_survival(model, subset, tmat)
    }, mc.cores = n_cpu)
    time_stop <- Sys.time()

    list(
        time = as.numeric(time_stop - time_start),
        pred = NULL
    )
}

fit_mstate <- function(data) {
    covariates_expanded <- attr(data, "covariates_expanded")

    formula <- as.formula(paste0(
        "Surv(Tstart, Tstop, status) ~ ",
        paste(covariates_expanded, collapse = "+"),
        " + cluster(id) + strata(trans)"
    ))

    time_start <- Sys.time()
    model <- coxph(formula, data, method = "breslow")
    time_stop <- Sys.time()
    list(
        model = model,
        time = as.numeric(time_stop - time_start)
    )
}

predict_mstate <- function(model, newdata, tmat) {
    newdata <- newdata %>%
        rename(strata = trans) %>%
        select(strata, attr(newdata, "covariates_expanded"))

    time_start <- Sys.time()
    pred <- msfit(model, newdata, trans = tmat) %>%
        probtrans(predt = 0) %>%
        .[[1]]
    time_stop <- Sys.time()

    list(
        time = time,
        pred = pred
    )
}

predict_mstate_parallel <- function(model, newdata, tmat, n_cpu) {
    # Get unique IDs and partition them into N groups
    ids <- unique(newdata$id)
    partitions <- split(ids, cut(seq_along(ids), n_cpu, labels = FALSE))

    # Apply function in parallel
    time_start <- Sys.time()
    results <- parallel::mclapply(partitions, function(id_subset) {
        subset <- newdata %>% filter(id %in% id_subset)
        predict_mstate(model, subset, tmat)
    }, mc.cores = n_cpu)
    time_stop <- Sys.time()

    list(
        time = as.numeric(time_stop - time_start),
        pred = NULL
    )
}

results_out <- function(data, file = "results.csv") {
    exists <- file.exists(file) && file.info(file)$size != 0
    write.table(as.data.frame(data), file,
        sep = ",", row.names = FALSE,
        col.names = !exists, append = exists
    )
}

pipeline_survival <- function(n_elements, n_train, n_covs, n_test, n_cpu) {
    g <- poset_hasse(n_elements)
    w <- hasse_walks(g, n_train + n_test)
    c <- hasse_covariates(g, w, n_covs)
    tmat <- hasse_tmat(g)

    df <- prep_hasse_counting(w, c, tmat)

    df_train <- df %>% filter(id %in% 1:n_train)
    fit <- fit_survival(df_train)

    if (n_test > 0) {
        df_test <- df %>% filter(id %in% (n_train + 1):(n_train + n_test))
        if (n_cpu > 1) {
            pred <- predict_survival_parallel(fit$model, df_test, tmat, n_cpu)
        } else {
            pred <- predict_survival(fit$model, df_test, tmat)
        }
    } else {
        pred <- list(pred = NULL, time = 0)
    }

    list(
        type = "survival",
        model = fit$model,
        pred = pred$pred,
        time_fit = fit$time,
        time_pred = pred$time
    )
}

pipeline_mstate <- function(n_elements, n_train, n_covs, n_test, n_cpu) {
    g <- poset_hasse(n_elements)
    w <- hasse_walks(g, n_train + n_test)
    c <- hasse_covariates(g, w, n_covs)
    tmat <- hasse_tmat(g)

    df <- prep_hasse_mstate(w, c, tmat)

    df_train <- df %>% filter(id %in% 1:n_train)
    fit <- fit_mstate(df_train)

    if (n_test > 0) {
        df_test <- df %>% filter(id %in% (n_train + 1):(n_train + n_test))
        if (n_cpu > 1) {
            pred <- predict_mstate_parallel(fit$model, df_test, tmat, n_cpu)
        } else {
            pred <- predict_mstate(fit$model, df_test, tmat)
        }
    } else {
        pred <- list(pred = NULL, time = 0)
    }

    list(
        type = "mstate",
        model = fit$model,
        pred = pred$pred,
        time_fit = fit$time,
        time_pred = pred$time
    )
}

pipeline <- list(
    "survival" = pipeline_survival,
    "mstate" = pipeline_mstate
)

# RUN ##########################################################################

# Test 1 - Fixed Testing Static Training #######################################
test1 <- function(n_elements, n_covs, n_test, n_rep = 3) {
    param_models <- c("survival", "mstate")
    param_n_train <- c(2^8, 2^9, 2^10, 2^11)
    param_n_rep <- 1:n_rep

    param_grid <- expand.grid(
        test = "1",
        model = param_models,
        n_elements = n_elements,
        n_train = param_n_train,
        n_covs = n_covs,
        n_test = n_test,
        n_cpu = 1,
        n_rep = param_n_rep
    )

    params <- list(
        model = "survival",
        n_elements = 4,
        n_train = 1024,
        n_covs = 8,
        n_test = 16,
        n_cpu = 1,
        n_rep = 1
    )

    lapply(seq_len(nrow(param_grid)), function(i) {
        params <- param_grid[i, ]
        print(params)
        tryCatch(
            {
                result <- pipeline[[params$model]](
                    n_elements = params$n_elements,
                    n_train = params$n_train,
                    n_covs = params$n_covs,
                    n_test = params$n_test,
                    n_cpu = params$n_cpu
                )
            },
            error = function(e) {
                result <- list(
                    time_fit = -1,
                    time_pred = -1
                )
            }
        )

        output <- c(
            params,
            list(
                time_fit = result$time_fit,
                time_pred = result$time_pred
            )
        )

        results_out(output)
    })
}

# Test 2 - Fit Scaling  ########################################################
test2 <- function(n_rep = 3) {
    param_models <- c("survival", "mstate")
    param_n_elements <- c(3, 4, 5)
    param_n_train <- c(2^9, 2^10, 2^11, 2^12, 2^13, 2^14, 2^15)
    param_n_covs <- c(2, 4, 8, 16)
    param_n_rep <- 1:n_rep

    param_grid <- expand.grid(
        test = "2",
        model = param_models,
        n_elements = param_n_elements,
        n_train = param_n_train,
        n_covs = param_n_covs,
        n_test = 0,
        n_cpu = 1,
        n_rep = param_n_rep
    )

    lapply(seq_len(nrow(param_grid)), function(i) {
        params <- param_grid[i, ]
        print(params)
        tryCatch(
            {
                result <- pipeline[[params$model]](
                    n_elements = params$n_elements,
                    n_train = params$n_train,
                    n_covs = params$n_covs,
                    n_test = params$n_test,
                    n_cpu = params$n_cpu
                )
            },
            error = function(e) {
                result <- list(
                    time_fit = NULL,
                    time_pred = NULL
                )
            }
        )

        output <- c(
            params,
            list(
                time_fit = result$time_fit,
                time_pred = result$time_pred
            )
        )

        results_out(output)
    })
}

# Test 3 - Predict Parallelisation  ############################################
test3 <- function(n_elements = 3, n_train = 2^13, n_covs = 4, n_rep = 3) {
    param_models <- c("survival", "mstate")
    param_n_test <- c(2^4, 2^5, 2^6, 2^7, 2^8)
    param_n_cpu <- c(2^0, 2^1, 2^2, 2^3, 2^4)
    param_n_rep <- 1:n_rep

    param_grid <- expand.grid(
        test = "3",
        model = param_models,
        n_elements = n_elements,
        n_train = n_train,
        n_covs = n_covs,
        n_test = param_n_test,
        n_cpu = param_n_cpu,
        n_rep = param_n_rep
    ) %>% arrange(n_test)

    lapply(seq_len(nrow(param_grid)), function(i) {
        params <- param_grid[i, ]
        print(params)
        result <- list(
            time_fit = -1,
            time_pred = -1
        )
        tryCatch(
            {
                result <- pipeline[[params$model]](
                    n_elements = params$n_elements,
                    n_train = params$n_train,
                    n_covs = params$n_covs,
                    n_test = params$n_test,
                    n_cpu = params$n_cpu
                )
            },
            error = function(e) NULL
        )

        output <- c(
            params,
            list(
                time_fit = result$time_fit,
                time_pred = result$time_pred
            )
        )

        results_out(output)
    })
}
