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
library(qs)

source("mythoslib/convert.R")
source("mythoslib/hasse.R")


# Set seed for reproducibility
set.seed(123)

# FITTING ######################################################################
fit_survival <- function(data) {
    # Dynamically create formula using the covariate names
    covs <- attr(data, "covariates")
    formula <- as.formula(
        paste0("Surv(tstart, tstop, to) ~ ", paste(covs, collapse = "+"))
    )

    # Fit the Cox PH model, tracking how long it takes
    time_start <- Sys.time()
    model <- coxph(formula, data, istate = from, id = id)
    time_stop <- Sys.time()

    # Return the model and timings
    list(
        model = model,
        time = as.numeric(time_stop - time_start, units = "secs")
    )
}

fit_mstate <- function(data) {
    # Dynamically create formula using the expanded covariate names
    covs_expanded <- attr(data, "covariates_expanded")
    formula <- as.formula(paste0(
        "Surv(Tstart, Tstop, status) ~ ",
        paste(covs_expanded, collapse = "+"),
        " + cluster(id) + strata(trans)"
    ))

    # Fit the Cox PH model, tracking how long it takes
    time_start <- Sys.time()
    model <- coxph(formula, data, method = "breslow")
    time_stop <- Sys.time()

    # Return the model and timings
    list(
        model = model,
        time = as.numeric(time_stop - time_start, units = "secs")
    )
}

fit_flexsurv <- function(data, tmat, knots = 2) {
    # Dynamically create formula using the covariate names
    covs <- attr(data, "covariates")
    formula <- as.formula(
        paste0("Surv(Tstart, Tstop, status) ~ ", paste(covs, collapse = "+"))
    )

    # Identify the unique transitions
    transitions <- tmat[!is.na(tmat)]

    time_start <- Sys.time()
    # Create cause-specific models
    models <- by(data, data$trans, function(subset) {
        flexsurvspline(
            formula,
            data = subset,
            k = knots
        )
    })
    model <- do.call(fmsm, c(models, list(trans = tmat)))
    time_stop <- Sys.time()

    list(
        model = model,
        time = as.numeric(time_stop - time_start, units = "secs"),
        knots = knots
    )
}

# PREDICTION ###################################################################
predict_survival <- function(model, data, tmat) {
    time_start <- Sys.time()
    # Subset according to initial state
    pred <- by(data, data$from, function(subset) {
        # For each initial state, compose the corresponding p0 vector
        p0 <- as.integer(model$states == subset[1, ]$from)

        # Select only the covariates
        subset <- subset %>% select(id, attr(data, "covariates"))

        # Call survfit
        survfit(model, newdata = subset, p0 = p0)
    })
    time_stop <- Sys.time()

    list(
        predictions = pred,
        time = as.numeric(time_stop - time_start, units = "secs")
    )
}

predict_mstate <- function(model, data, tmat) {
    time_start <- Sys.time()
    pred <- msfit(object = model, newdata = data, trans = tmat) %>%
        probtrans(predt = 0)
    time_stop <- Sys.time()

    list(
        predictions = pred,
        time = as.numeric(time_stop - time_start, units = "secs")
    )
}

predict_flexsurv <- function(model, data, tmat, tgrid) {
    time_start <- Sys.time()
    pred <- lapply(seq_len(nrow(data)), function(row) {
        pmatrix.fs(
            model,
            newdata = data[row, ],
            trans = tmat,
            t = tgrid,
            tidy = TRUE
        )
    })
    time_stop <- Sys.time()

    list(
        predictions = pred,
        time = as.numeric(time_stop - time_start, units = "secs")
    )
}

predict_parallel <- function(predict_func, data, cpus = 1) {
    # Partition data into $cpus total groups
    if (cpus > 1) {
        partitions <- data %>%
            group_by(cut(seq_along(id), cpus, labels = FALSE)) %>%
            group_split()
        partitions <- lapply(partitions, function(partition) {
            attr(partition, "covariates") <- attr(data, "covariates")
            partition
        })
    } else {
        partitions <- list(data)
    }

    # Parallelise prediction over $cpus processes
    time_start <- Sys.time()
    pred <- parallel::mclapply(partitions, function(partition) {
        predict_func(partition)
    }, mc.cores = cpus)
    time_stop <- Sys.time()

    # Return predictions and timings
    list(
        predictions = pred,
        time = as.numeric(time_stop - time_start, units = "secs")
    )
}

# PIPELINES ####################################################################
pipeline_survival <- function(walks, covariates, tmat,
                              n_train, n_test, cpus = 1) {
    # 1. Prepare the dataset
    data <- convert_hasse_counting(walks, covariates, tmat)
    df_train <- data %>% filter(id %in% 1:n_train)
    df_test <- data %>%
        filter(id %in% (n_train + 1):(n_train + n_test)) %>%
        filter(tstart == 0)

    # 2. Fit the model
    fit <- fit_survival(df_train)

    # 3. Predict
    if (n_test > 0) {
        pred <- predict_parallel(
            function(df) predict_survival(fit$model, df, tmat),
            df_test, cpus
        )
    } else {
        pred <- list(predictions = NULL, time = 0)
    }

    list(
        model = fit$model,
        predictions = pred$predictions,
        time_fit = fit$time,
        time_pred = pred$time
    )
}

pipeline_mstate <- function(walks, covariates, tmat,
                            n_train, n_test, cpus = 1) {
    # 1. Prepare the dataset
    data <- convert_hasse_msdata(walks, covariates, tmat)

    # Split the data into training and testing sets
    df_train <- data %>% filter(id %in% 1:n_train)
    df_test <- data %>%
        filter(id %in% (n_train + 1):(n_train + n_test)) %>%
        rename(strata = trans)

    # 2. Fit the model

    fit <- fit_mstate(df_train)

    # 3. Predict
    if (n_test > 0) {
        pred <- predict_parallel(
            function(df) predict_mstate(fit$model, df, tmat),
            df_test, cpus
        )
    } else {
        pred <- list(predictions = NULL, time = 0)
    }

    list(
        model = fit$model,
        pred = pred$predictions,
        time_fit = fit$time,
        time_pred = pred$time
    )
}

pipeline_flexsurv <- function(walks, covariates, tmat,
                              n_train, n_test, cpus = 1) {
    # 1. Prepare the dataset
    data <- convert_hasse_msdata(walks, covariates, tmat, expand = FALSE)

    # Split the data into training and testing sets
    df_train <- data %>% filter(id %in% 1:n_train)
    df_test <- data %>%
        filter(id %in% (n_train + 1):(n_train + n_test)) %>%
        filter(Tstart == 0)

    # Obtain the time grid for prediction
    tgrid <- df_train %>%
        filter(status == 1) %>%
        pull(Tstop) %>%
        unique() %>%
        sort()

    # 2. Fit the model - iterating K downwards until convergence
    k <- 5
    while (k > 0) {
        tryCatch(
            {
                print(paste("Fitting flexsurvspline with", k, "knots"))
                fit <- fit_flexsurv(df_train, tmat, knots = k)
                break
            },
            error = function(e) {
                msg <- conditionMessage(e)
                err_vmin <- grepl("initial value in 'vmmin' is not finite", msg)
                if (!err_vmin || k == 1) {
                    stop(e)
                }
                print(msg)
            }
        )
        k <- k - 1
    }

    # 3. Predict
    if (n_test > 0) {
        pred <- predict_parallel(
            function(df) predict_flexsurv(fit$model, df, tmat, tgrid),
            df_test, cpus
        )
    } else {
        pred <- list(predictions = NULL, time = 0)
    }

    list(
        model = fit$model,
        pred = pred$predictions,
        time_fit = fit$time,
        time_pred = pred$time,
        knits = fit$knots
    )
}

# RUN ##########################################################################
pipelines <- list(
    survival = pipeline_survival,
    mstate = pipeline_mstate,
    flexsurv = pipeline_flexsurv
)

generate_hasse <- function(elements, n_covs, n_train, n_test, times) {
    graph <- poset_hasse(elements)
    walks <- hasse_walks(graph, n_train + n_test, times)
    covariates <- hasse_covariates(graph, walks, n_covs)
    tmat <- hasse_tmat(graph)

    list(
        graph = graph,
        walks = walks,
        covariates = covariates,
        tmat = tmat
    )
}

bench <- function(model, elements, n_covs, n_train, n_test, times, cpus, rep) {
    # Generate the hasse graph and associated data
    print("Generating hasse graph.")
    time_start <- Sys.time()
    hasse_data <- generate_hasse(elements, n_covs, n_train, n_test, times)
    time_stop <- Sys.time()
    print(time_stop - time_start)

    # Run the pipeline
    print("Running pipeline")
    time_start <- Sys.time()
    result <- tryCatch(
        {
            pipelines[[model]](
                hasse_data$walks,
                hasse_data$covariates,
                hasse_data$tmat,
                n_train,
                n_test,
                cpus
            )
        },
        error = function(e) {
            print(conditionMessage(e))
            list(
                error = e,
                error_message = conditionMessage(e),
                hasse = hasse_data
            )
        }
    )
    time_stop <- Sys.time()
    print(time_stop - time_start)

    result[["parameters"]] <- list(
        model = model,
        elements = elements,
        n_covs = n_covs,
        n_train = n_train,
        n_test = n_test,
        times = times,
        cpus = cpus,
        rep = rep
    )

    result
}

main <- function(result_dir = "data/results", model_dir = "data/models") {
    # Access command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) == 0) {
        stop("No arguments provided")
    }

    # Obtain test type
    test_type <- noquote(args[1])

    # Obtain model type
    model <- noquote(args[2])
    if (!(model %in% c("survival", "mstate", "flexsurv"))) {
        stop(paste("Invalid model type", model))
    }

    # Obtain number of elements
    elements <- as.integer(args[3])
    if (is.na(elements) || elements < 2) {
        stop(paste("Invalid number of elements", elements))
    }

    # Obtain number of covariates
    n_covs <- as.integer(args[4])
    if (is.na(n_covs) || n_covs < 1) {
        stop(paste("Invalid number of covariates", n_covs))
    }

    # Obtain number of training and testing samples
    n_train <- as.integer(args[5])
    n_test <- as.integer(args[6])
    if (is.na(n_train) || is.na(n_test) || n_train < 1 || n_test > n_train) {
        stop(paste("Invalid number of samples", n_train, "/", n_test))
    }

    # Obtain number of distinct event times
    # Allow -1 to indicate arbitrary event times
    times <- as.integer(args[7])
    if (is.na(times) || !(times > 0 || times == -1)) {
        stop(paste("Invalid number of event times", times))
    }

    # Obtain number of CPU cores
    cpus <- as.integer(args[8])
    if (is.na(cpus) || cpus < 1) {
        stop(paste("Invalid number of CPU cores", cpus))
    }

    # Obtain number of repetitions
    rep <- as.integer(args[9])
    if (is.na(rep) || rep < 1) {
        stop(paste("Invalid repetition number", rep))
    }

    # Print parameters
    print("Parameters:")
    print(paste(
        "Test type:", test_type, "Model:", model, "Elements:", elements,
        "Covariates:", n_covs, "Train:", n_train, "Test:", n_test,
        "Times:", times, "CPUs:", cpus, "Rep:", rep,
        sep = " "
    ))

    # Run the benchmark
    result <- bench(
        model, elements, n_covs, n_train,
        n_test, times, cpus, rep
    )

    print(paste("Time fit:", result$time_fit, "Time pred:", result$time_pred))
    # Save the results
    print("Saving results")
    if (!dir.exists(result_dir)) {
        dir.create(result_dir, recursive = TRUE)
    }
    if (!dir.exists(model_dir)) {
        dir.create(model_dir, recursive = TRUE)
    }
    filename <- paste(
        "results_", test_type, "_", model, "_", elements, "_",
        n_covs, "_", n_train, "_", n_test, "_", cpus, "_", rep,
        sep = ""
    )
    result$parameters[["test"]] <- test_type
    write_csv(data.frame(list(
        parameters = result$parameters,
        time_fit = result$time_fit,
        time_pred = result$time_pred
    )), file.path(result_dir, paste0(filename, ".csv")))
    qsave(result, file.path(model_dir, paste0(filename, ".qs")))
    print("Done")
}

main()
