# HEADER #######################################################################
#
# Author: Stelios Boulitsakis Logothetis
# Email: sb2690@medschl.cam.ac.uk
#
# Date: 2025-03-13
# Name:	mythoslib/survfit.R
# Description: Helper functions for working with survfit.coxphms objects
#
# SETUP ########################################################################
library(dplyr)
library(abind)
library(data.table)

# FUNCTIONS ####################################################################

#' Merge a list of survfit.coxphms objects into a single survfit object
#'
#' This function takes a list of survfit.coxphms objects and merges them into a
#' single survfit object as neatly as possible. Most of the components of the
#' input objects are invariant as they are derived from the fitting set, not
#' the prediction set. The only components that are not invariant are
#' 'pstate', 'cumhaz', 'newdata', and 'call'
#' @param fit_list A list of survfit.coxphms objects. Names are assumed to be
#'                the names of the initial state (p0) for each survfit object
#' @return A survfit object containing all the data from the input objects
#'
#' @export
merge_survfit_list <- function(fit_list) {
    # 'pstate' - Estimated probability of each state at each time
    # Stack 'pstate' elements along the 2nd dimension
    pstate <- abind(lapply(fit_list, function(fit) fit$pstate), along = 2)


    # 'cumhaz' - Cumulative hazard for each possible transition
    # Stack identically to 'pstate'
    cumhaz <- abind(lapply(fit_list, function(fit) fit$cumhaz), along = 2)

    # 'newdata' - Covariate values for the curves
    # Combine newdata data frames by row-binding
    newdata <- do.call(rbind, lapply(fit_list, function(fit) fit$newdata))

    # 'call' - Image of the call that produced the object
    # Collate model fitting calls into a list
    call <- lapply(fit_list, function(fit) fit$call)

    # Create an additional "istate" attribute to distinguish the p0
    # of the input objects
    n_ids <- unlist(lapply(fit_list, function(fit) max(dim(fit$pstate)[2], 0)))
    istate <- rep(names(fit_list), n_ids)

    # For components assumed to be invariant, take from the first non-null fit.
    first_fit <- (fit_list[!sapply(fit_list, is.null)])[[1]]

    # Build the merged result list.
    merged_fit <- list(
        # Number of observations in each curve
        n            = first_fit$n, # (invariant)
        # Time points at which the curve has a step
        time         = first_fit$time, # (invariant)
        # Number of subjects at risk at t
        n.risk       = first_fit$n.risk, # (invariant)
        # Number of events that occur at time t
        n.event      = first_fit$n.event, # (invariant)
        # Number of subjects who exit the risk set, without an event, at time t
        n.censor     = first_fit$n.censor, # (invariant)
        # Estimated probability of each state at each time
        pstate       = pstate,
        # Number of occurrences of each transition at time t
        n.transition = first_fit$n.transition, # (invariant)
        # Number of IDs that were used to fit the model
        n.id         = first_fit$n.id, # (invariant)
        # Cumulative hazard for each possible transition
        cumhaz       = cumhaz,
        # Initial state proportions - Impossible to merge
        p0           = first_fit$p0,
        # Number of elements of the time vector corresponding to each curve
        strata       = first_fit$strat, # (invariant)
        # Time of first event
        start.time   = first_fit$start.time,
        # Transition matrix of fitting set
        transitions  = first_fit$transitions,
        # Names of states
        states       = first_fit$states,
        # Type of data used to fit the model (expected to be "mcounting")
        type         = first_fit$type,
        # Starting time for the curves
        t0           = first_fit$t0,
        # Covariate values for the curves
        newdata      = newdata,
        # Images of the calls that produced the objects
        call         = call
    )

    # Set the class attribute to retain the original class
    class(merged_fit) <- class(first_fit)

    # Include the "istate" attribute
    attr(merged_fit, "istate") <- istate

    merged_fit
}

#' Convert a survfit.coxphms object to a long-format data frame of state
#' occupancy probabilities.
#'
#' This function takes a survfit.coxphms object and converts it to a tidy
#' long-format data frame. The resulting data frame will have one row per
#' combination of ID, state, and time point, and will contain the corresponding
#' state occupancy probability for that sample. This accounts for which stratum
#' and initial state each ID belongs to. The input object is assumed to either
#' have a single p0 initial state, or to be the result of a merge_survfit_list
#' operation (and therefore has an "istate" attribute).
#'
#' @param obj A survfit.coxphms object. Must have an "istate" attribute or
#'           a p0 matrix with a single column of all 1 elements.
#' @param strata_column The name of the column in the newdata data frame that
#'                     contains the strata values for each ID.
#' @param id_column The name of the column in the newdata data frame that
#'                contains the ID values.
#' @return A tidy long-format data frame containing the state occupancy
#'        probabilities for each ID, state, and time point. Includes attributes:
#' - "strata": A data frame containing the strata values for each ID
#' - "states": A character vector containing the names of the states
#' - "n.strata": A table containing the number of IDs in each stratum
#' - "istate": A character vector containing the initial state for each ID
#'
#' @export
tidy_pstate <- function(obj, strata_column, id_column = "id") {
    # Identify the strata values
    strata_values <- names(obj$strata) %>%
        sapply(function(x) sub(".*=", "", x)) %>%
        as.numeric()

    # Expand the strata labels along the "time" vector
    # The first strata[[1]] time points belong to the first stratum, the next
    # obj$strata[[2]] time points belong to the second stratum, and so on
    time_strata <- rep(strata_values, obj$strata)

    # Melt the pstate matrix into a long format
    dt <- as.data.table(obj$pstate)
    setnames(dt, names(dt), c("time_idx", "id_idx", "state_idx", "value"))

    # Map the indices to proper values
    dt$time <- obj$time[dt$time_idx]
    dt$id <- obj$newdata[dt$id_idx, id_column]
    dt$state <- obj$states[dt$state_idx]

    # Filter - keep only rows where the time stratum equals the ID's stratum
    mask <- time_strata[dt$time_idx] == obj$newdata[dt$id_idx, strata_column]

    # Keep only the relevant rows and columns
    dt <- dt[mask, c("id", "state", "time", "value")]

    attr(dt, "strata") <- obj$newdata[, c(id_column, strata_column)]
    attr(dt, "states") <- obj$states
    attr(dt, "n.strata") <- obj$newdata[, strata_column] %>% table()

    if ("istate" %in% names(attributes(obj))) {
        attr(dt, "istate") <- data.frame(
            id = obj$newdata[, id_column],
            istate = attr(obj, "istate")
        )
    } else {
        # Identify the first column to have all 1 elements and return its name
        attr(dt, "istate") <- obj$p0 %>%
            apply(2, function(col) all(col == 1)) %>%
            obj$states[.]
    }

    dt
}

#' Given a survfit.coxphms object, recovers the ordering of transitions.
#'
#' Multi-state Cox models in the survival package provide the ordering of
#' transitions and coefficients in their 'cmap' attribute. E.g.:
#' r$> model$cmap
#'     1:2 2:3 4:3 1:4
#' C_1   1   5   9  13
#' C_2   2   6  10  14
#' C_3   3   7  11  15
#' C_4   4   8  12  16
#' The ordering of transitions is necessary for interpreting the estimated
#' cumulative hazards in a survfit object. This helper function recovers this
#' ordering and returns a list containing "from" and "to" character vectors
#'
#' @param obj A survfit.coxphms object
#' @return A list containing "from" and "to" vectors
#'
#' @export
transition_ordering <- function(obj) {
    obj_transitions <- obj$transitions %>%
        as.matrix() %>%
        .[, colnames(.) != "(censored)"]

    # Find nonzero indices
    nonzero_indices <- which(obj_transitions != 0, arr.ind = TRUE)

    # Construct from and to vectors
    list(
        from = rownames(obj_transitions)[nonzero_indices[, 1]],
        to = colnames(obj_transitions)[nonzero_indices[, 2]]
    )
}

#' Convert a survfit.coxphms object to a long-format data frame of cumulative
#' transition probabilities.
#'
#' This function takes a survfit.coxphms object and converts it to a tidy
#' long-format data frame. The resulting data frame will have one row per
#' combination of ID, transition, and time point, and will contain the
#' corresponding transition probability for that sample. This accounts for which
#' stratum and initial state each ID belongs to. The input object is assumed to
#' either have a single p0 initial state, or to be the result of a
#' merge_survfit_list operation (and therefore has an "istate" attribute).
#'
#' @param obj A survfit.coxphms object. Must have an "istate" attribute or
#'           a p0 matrix with a single column of all 1 elements.
#' @param strata_column The name of the column in the newdata data frame that
#'                     contains the strata values for each ID.
#' @param id_column The name of the column in the newdata data frame that
#'                contains the ID values.
#' @return A tidy long-format data frame containing the transition probabilities
#'  for each ID, transition, and time point. Includes attributes:
#' - "strata": A data frame containing the strata values for each ID
#' - "states": A character vector containing the names of the states
#' - "n.strata": A table containing the number of IDs in each stratum
#' - "istate": A character vector containing the initial state for each ID
#'
#' @export
tidy_probtrans <- function(obj, strata_column, id_column = "id") {
    # Identify the strata values
    strata_values <- names(obj$strata) %>%
        sapply(function(x) sub(".*=", "", x)) %>%
        as.numeric()

    # Expand the strata labels along the "time" vector
    # The first strata[[1]] time points belong to the first stratum, the next
    # obj$strata[[2]] time points belong to the second stratum, and so on
    time_strata <- rep(strata_values, obj$strata)

    # Melt the cumhaz matrix into a long format
    dt <- as.data.table(obj$cumhaz)
    setnames(dt, names(dt), c("time_idx", "id_idx", "trans_idx", "cumhaz"))

    # Map the indices to proper values
    dt$time <- obj$time[dt$time_idx]
    dt$id <- obj$newdata[dt$id_idx, id_column]

    transitions <- transition_ordering(obj)
    trans_names <- paste(transitions$from, transitions$to, sep = ":")
    dt$transition <- trans_names[dt$trans_idx]
    dt$from <- transitions$from[dt$trans_idx]
    dt$to <- transitions$to[dt$trans_idx]

    # Filter - keep only rows where the time stratum equals the ID's stratum
    mask <- time_strata[dt$time_idx] == obj$newdata[dt$id_idx, strata_column]

    # Keep only the relevant rows and columns
    dt <- dt[mask, c("id", "time", "from", "to", "transition", "cumhaz")]

    attr(dt, "strata") <- obj$newdata[, c(id_column, strata_column)]
    attr(dt, "states") <- obj$states
    attr(dt, "transitions") <- transitions
    attr(dt, "n.strata") <- obj$newdata[, strata_column] %>% table()

    if ("istate" %in% names(attributes(obj))) {
        attr(dt, "istate") <- data.frame(
            id = obj$newdata[, id_column],
            istate = attr(obj, "istate")
        )
    } else {
        # Identify the first column to have all 1 elements and return its name
        attr(dt, "istate") <- obj$p0 %>%
            apply(2, function(col) all(col == 1)) %>%
            obj$states[.]
    }

    dt
}

#' Convert a survfit summary object to a tidy long-format data frame.
#'
#' Polyfill for the `data.frame` argument in the `summary.survfitms` function
#' Converts the given summary object of a survfit object to a dataframe.
#'
#' @param fit A survfit object's summary (output of `summary.survfitms`)
#' @return A tidy long-format data frame containing the summary data.
#'
#' @export
survfit_summary_to_df <- function(fit) {
    # We will copy over the following columns
    cols <- match(
        c("time", "pstate", "std.err", "lower", "upper"),
        names(fit),
        nomatch = 0
    )

    # Identify the number of points in the survival curves (length of pstate)
    # And the number of states (width of pstate)
    n_points <- dim(fit$pstate)[1]
    n_states <- dim(fit$pstate)[2]

    # Convert matrices to vectors
    ndata <- lapply(fit[cols], as.vector)

    # Define columns
    ndata$time <- rep(ndata$time, n_states)
    ndata$strata <- rep(fit$strata, n_states)
    ndata$state <- rep(fit$states, each = n_points)

    data.frame(ndata)
}
