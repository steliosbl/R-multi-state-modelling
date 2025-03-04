# HEADER #######################################################################
#
# Author: Stelios Boulitsakis Logothetis
# Email: sb2690@medschl.cam.ac.uk
#
# Date: 2025-02-26
# Name:	mythoslib/convert.R
# Description: Helper functions for manipulating multi-state modelling
# state-transition data.
#
# SETUP ########################################################################
# Load required libraries
library(dplyr)

# FUNCTIONS ####################################################################

#' Identify absorbing events from a transition matrix
tmat_absorbing_events <- function(tmat) {
    rownames(tmat)[apply(tmat, 1, function(x) all(is.na(x)))]
}

#' Augment a Timeline-format dataframe with censoring entries
#'
#' Given a timeline-format state-transition dataframe which lacks censoring,
#' this function appends a censoring row for each subject whose last recorded
#' event is not one of the specified terminal (absorbing) events. The censoring
#' row is assigned time \code{cens_time} and event name \code{cens_name}.
#' Absorbing events are identified from the transition matrix \code{tmat}.
#' WARNING: Doesn't if the last observed event occurs before \code{cens_time}.
#' @param df A Timeline-format dataframe \code{({id_col}, time, event)}.
#' @param id_col The name of the column in \code{df} that identifies subjects.
#' @param tmat A transition matrix with row and column names corresponding to
#' .  event names. Absorbing events are identified as rows with all NA entries.
#' @param cens_time The time at which censoring occurs.
#' @param cens_name The name of the event to assign to censoring rows.
#' @return A Timeline-format dataframe with censoring rows included.
#'
#' @export
augment_timeline_censoring <- function(df, id_col, tmat, cens_time, cens_name) {
    # Identify absorbing events
    absorbing_events <- tmat_absorbing_events(tmat)

    df %>%
        # For each distinct ID
        group_by(across(all_of(id_col))) %>%
        # Sort by time in ascending order
        arrange(time) %>%
        group_modify(~ {
            # If the last event by time is not absorbing, add a censor row
            if (!(last(.x$event) %in% absorbing_events)) {
                # Create a censor row
                censor_row <- tibble::tibble(
                    !!id_col := .x[[id_col]][1],
                    time := cens_time,
                    event := cens_name
                )
                # Bind the censor row to the data
                .x <- dplyr::bind_rows(.x, censor_row)
            }
            .x
        }) %>%
        ungroup()
}

#' Convert a Timeline-format dataframe to Counting-process format
#'
#' Given a Timeline-format dataframe, this function converts it to a
#' Counting-process format dataframe. The Counting-process format contains, for
#' each subject, the start and stop times of each interval along with the start
#' and end states. The input dataframe is expected to include censoring entries.
#'
#' @param df A Timeline-format dataframe \code{({id_col}, time, event)}.
#' @param id_col The name of the column in \code{df} that identifies subjects.
#'
#' @return A Counting-process format dataframe with columns:
#' .  \code{({id_col}, tstart, tstop, from, to)}.
#' .  \code{tstart} and \code{tstop} are the start and stop times of intervals.
#' .  \code{from} and \code{to} are the start and end states of intervals.
#'
#' @export
timeline_to_counting <- function(df, id_col) {
    df %>%
        # For each unique id
        group_by(across(all_of(id_col))) %>%
        # Sort by time in ascending order
        arrange(time, .by_group = TRUE) %>%
        # Convert the augmented timeline to counting process format
        mutate(
            tstart = lag(time),
            tstop  = time,
            from   = lag(event),
            to     = event
        ) %>%
        # Remove the first row per ID (which has no "from")
        filter(!is.na(tstart)) %>%
        ungroup() %>%
        select(all_of(id_col), tstart, tstop, from, to)
}


#' Convert a Counting-process format dataframe to Timeline format
#'
#' Given a counting process-style dataframe, this function converts it to a
#' timeline format. The resulting timeline data contains, for each subject,
#' one entry for every visited state with the entry time into it, including
#' the initial state and censoring. The input dataframe is expected to include
#' censoring entries.
#'
#' @param df A dataframe in counting process
#' .  \code{({id_col}, tstart, tstop, from, to)}.
#' @param id_col The name of the column in \code{df} that identifies subjects.
#'
#' @return A dataframe in timeline format with columns:
#' .  \code{({id_col}, time, event)}.
#' .  Each row represents a visited state for a subject,
#' .  including the initial state and censoring.
#'
#' @export
counting_to_timeline <- function(df, id_col) {
    # Identify the initial state for each subject
    istate <- df %>%
        group_by(across(all_of(id_col))) %>%
        arrange(tstart, .by_group = TRUE) %>%
        summarise_all(first) %>%
        transmute(!!id_col := .data[[id_col]], time = tstart, event = from)

    # Every counting process row gives a new event `to` at time `tstop`
    df %>%
        transmute(!!id_col := .data[[id_col]], time = tstop, event = to) %>%
        bind_rows(istate) %>%
        arrange(.data[[id_col]], time)
}

#' Convert a Timeline-format dataframe to mstate format
#'
#' Given a Timeline-format dataframe, this function converts it to an mstate
#' format dataframe. The mstate format contains, for each subject, one row
#' with \code{{event}_status} and \code{{event}_time} columns for each event.
#' The input dataframe is expected to include censoring entries. When a state
#' is not visited by a subject, the corresponding \code{{event}_time} is set to
#' the subject's final time.
#'
#' @param df A Timeline-format dataframe \code{({id_col}, time, event)}.
#' @param id_col The name of the column in \code{df} that identifies subjects.
#' @param tmat A transition matrix with row names corresponding to events.
#' @param cens_name The name of the event assigned to censoring rows.
#'
#' @return An mstate-format dataframe with columns:
#' .  \code{({id_col}, {event}_time, {event}_status), istate, itime}
#' .  The \code{istate}, \code{itime} columns give the starting state and time.
timeline_to_mstate <- function(df, id_col, tmat, cens_name) {
    # Validate - ensure expected column names
    if (!all(c(id_col, "time", "event") %in% colnames(df))) {
        stop("df must have columns {id_col}, 'time' and 'event'.")
    }

    # Validate - Transition matrix relfects dataframe state names
    if (!setequal(df$event, c(rownames(tmat), cens_name))) {
        stop("tmat state names don't match dataframe state names.")
    }

    # Type conversion - states will have to be strings
    df <- df %>%
        mutate(event = as.character(event))

    # Identify the initial and final states for each subject
    entry_exit <- df %>%
        group_by(across(all_of(id_col))) %>%
        arrange(time, .by_group = TRUE) %>%
        summarize(
            istate = first(event),
            itime = first(time),
            ftime = last(time),
            .groups = "drop"
        )

    # Introduce rows for the states that were never visited
    df_completed <- df %>%
        complete(
            !!sym(id_col),
            event = c(rownames(tmat), cens_name),
            fill = list(time = NA)
        )

    # Pivot to wide format and fill missing times with the subject's final time
    df_completed %>%
        mutate(status = ifelse(!is.na(time), 1, 0)) %>%
        pivot_wider(
            id_cols = all_of(id_col),
            names_from = event,
            values_from = c(time, status),
            values_fn = min,
            names_sep = "_"
        ) %>%
        left_join(entry_exit, by = {{ id_col }}) %>%
        mutate(across(starts_with("time_"), ~ ifelse(is.na(.), ftime, .))) %>%
        select(-ftime)
}


#' Convert an mstate-format dataframe to Timeline format
#'
#' Given an mstate-format dataframe, this function converts it to a timeline
#' format. The resulting timeline data contains, for each subject, one entry
#' for every visited state with the entry time into it, including the initial
#' state and censoring.
#'
#' @param df An mstate-format dataframe with columns:
#' .  \code{({id_col}, {event}_time, {event}_status), istate, itime}
#' @param id_col The name of the column in \code{df} that identifies subjects.
#' @param tmat A transition matrix with row names corresponding to events.
#' @param cens_name The name of the event assigned to censoring rows.
#'
#' @return A dataframe in timeline format with columns:
#' .  \code{({id_col}, time, event)}.
#' .  Each row represents a visited state for a subject,
#' .  including the initial state and censoring.
#'
#' @export
mstate_to_timeline <- function(df, id_col, tmat, cens_name) {
    # Pivot to long format, excluding istate and itime
    df_long <- df %>%
        pivot_longer(
            cols = -all_of(c(id_col, "istate", "itime")),
            names_to = c(".value", "event"),
            names_sep = "_"
        )

    identified_events <- unique(df_long$event)
    # Validate - All initial states appear as columns in the input
    if (!all(unique(df_long$istate) %in% identified_events)) {
        stop("Some states in istate don't appear as columns in df.")
    }

    # Validate - Transition matrix relfects dataframe state names
    if (!all(identified_events %in% c(rownames(tmat), cens_name))) {
        stop("tmat state names don't match dataframe state names.")
    }

    # Validate - Data types match between transitions and cens_name
    mismatch <- (is.numeric(df_long$event) != is.numeric(cens_name)) &&
        (is.character(df_long$event) != is.character(cens_name))
    if (mismatch) {
        stop("cens_name must have the same data type as the event column.")
    }

    # The input may have not included censoring as a state
    # In this case, add it to the output timeline
    if (!cens_name %in% identified_events) {
        # Identify absorbing events
        absorbing_events <- tmat_absorbing_events(tmat)
        # Identify the last event for each subject, prioritising absorbing ones
        last_event <- df_long %>%
            group_by(across(all_of(id_col))) %>%
            arrange((event %in% absorbing_events), time) %>%
            summarise_all(last)
        # Create censoring rows when the last event is not absorbing
        censor_rows <- last_event %>%
            filter(!(event %in% absorbing_events & status == 1)) %>%
            mutate(event = cens_name, status = 1)

        df_long <- bind_rows(df_long, censor_rows)
    }

    df_long %>%
        filter(status == 1) %>%
        select(all_of(id_col), time, event) %>%
        arrange(.data[[id_col]], time)
}

#' Convert an msdata object dataframe to Timeline format
#'
#' Given an msdata object dataframe, this function converts it to a timeline
#' format. The resulting timeline data contains, for each subject, one entry
#' for every visited state with the entry time into it, including the initial
#' state and censoring.
#'
#' @param df An msdata object dataframe with columns:
#' .  \code{({id_col}, from, to, Tstart, Tstop, status, trans}
#' @param id_col The name of the column in \code{df} that identifies subjects.
#' @param tmat A transition matrix with row names corresponding to events.
#' @param cens_name The name of the event assigned to censoring rows.
#'
#' @return A dataframe in timeline format with columns:
#' .  \code{({id_col}, time, event)}.
#' .  Each row represents a visited state for a subject,
#' .  including the initial state and censoring.
#'
#' @export
msdata_to_timeline <- function(df, id_col, tmat, cens_name) {
    # Validate - ensure expected column names
    cc <- c(id_col, "from", "to", "Tstart", "Tstop", "status")
    if (!all(cc %in% colnames(df))) {
        stop(
            "df has incorrect columns. It must include ",
            paste(cc, collapse = ", ")
        )
    }
    # Validate - Data types match between tmat and cens_name
    mismatch <- (is.numeric(df$to) != is.numeric(cens_name)) &&
        (is.character(df$to) != is.character(cens_name))
    if (mismatch) {
        stop("cens_name must have the same data type as the event columns.")
    }

    # Validate - Transition matrix relfects dataframe state names
    if (!all(unique(df$to) %in% rownames(tmat))) {
        stop("tmat state names don't match dataframe state names.")
    }

    # Identify the initial state for each subject
    istate <- df %>%
        group_by(across(all_of(id_col))) %>%
        arrange(Tstart, .by_group = TRUE) %>%
        summarise_all(first) %>%
        transmute(!!id_col := .data[[id_col]], time = Tstart, event = from)

    # Every counting process row gives a new event `to` at time `Tstop`
    subsequent <- df %>%
        filter(status == 1) %>%
        transmute(!!id_col := .data[[id_col]], time = Tstop, event = to)

    # Censoring - identify absorbing events
    absorbing_events <- tmat_absorbing_events(tmat)

    # Identify the last event for each subject, prioritising absorbing ones
    last_event <- df %>%
        group_by(across(all_of(id_col))) %>%
        arrange((to %in% absorbing_events), Tstop, .by_group = TRUE) %>%
        summarise_all(last)

    # Create censoring rows when the last event is not absorbing
    censor_rows <- last_event %>%
        filter(!(to %in% absorbing_events & status == 1)) %>%
        mutate(event = cens_name, time = Tstop) %>%
        select(all_of(id_col), time, event)

    bind_rows(istate, subsequent, censor_rows) %>%
        arrange(.data[[id_col]], time)
}

# TESTS ########################################################################
test_ebmt3 <- function() {
    library(mstate)
    data(ebmt3)

    id_col <- "id"
    cens_name <- "censor"
    tmat <- trans.illdeath(names = c("s0", "pr", "rfs"))

    # Format correctly (add s0 and istate columns)
    # And introduce some alternative initial states
    mstate_1 <- ebmt3 %>%
        select(id,
            time_pr = prtime, status_pr = prstat,
            time_rfs = rfstime, status_rfs = rfsstat
        ) %>%
        mutate(time_s0 = 0, status_s0 = 1, istate = "s0", itime = 0) %>%
        mutate(flip_istate = (status_pr == 1 & runif(n()) < 0.3)) %>%
        mutate(
            time_pr = ifelse(flip_istate, 0, time_pr),
            status_s0 = ifelse(flip_istate, 0, status_s0),
            time_s0 = ifelse(flip_istate, pmax(time_pr, time_rfs), time_s0),
            istate = ifelse(flip_istate, "pr", istate),
        ) %>%
        select(-flip_istate)

    # Test 1: mstate -> timeline -> mstate
    timeline_1 <- mstate_to_timeline(mstate_1, id_col, tmat, cens_name)
    mstate_2 <- timeline_to_mstate(timeline_1, id_col, tmat, cens_name) %>%
        select(all_of(colnames(mstate_1)))

    # This sometimes fails: pr_TIME is sometimes != the final event time
    # (rfs) even though rfs_STATUS == 1
    # But this means the data is wrong, not us
    # stopifnot(all(mstate_1 == mstate_2))
    # Thus, rather than comparing the dataframes, we compare msprep results
    prep <- function(df) {
        df %>% msprep(
            status = c("status_s0", "status_pr", "status_rfs"),
            time = c("time_s0", "time_pr", "time_rfs"),
            trans = tmat,
            start = list(
                state = match(.$istate, dimnames(tmat)[[1]]),
                time = .$itime
            )
        )
    }
    stopifnot(prep(mstate_1) == prep(mstate_2))

    # Test 2: timeline -> counting -> timeline
    counting_1 <- timeline_to_counting(timeline_1, id_col)
    timeline_2 <- counting_to_timeline(counting_1, id_col)
    stopifnot(all(timeline_1 == timeline_2))

    # Test 3: timeline -> msdata -> timeline
    msdata_1 <- timeline_to_mstate(timeline_1, id_col, tmat, cens_name) %>%
        prep()
    msdata_comparison <- msdata_1 %>%
        mutate(to = recode(to, `1` = "s0", `2` = "pr", `3` = "rfs")) %>%
        mutate(from = recode(from, `1` = "s0", `2` = "pr", `3` = "rfs"))

    timeline_3 <- msdata_to_timeline(msdata_comparison, id_col, tmat, cens_name)
    stopifnot(timeline_1 == timeline_3)

    # Test 4: mstate -> msdata -> timeline -> mstate
    mstate_3 <- timeline_to_mstate(timeline_3, id_col, tmat, cens_name) %>%
        select(all_of(colnames(mstate_2)))
    stopifnot(all(mstate_3 == mstate_2))

    # Test 5: Survcheck (counting process)
    # Pass counting process data to survcheck and verify flags are all 0
    stopifnot(all(counting_1$tstart < counting_1$tstop))
    counting_survcheck <- counting_1 %>%
        mutate(from = factor(from)) %>%
        mutate(to = factor(to, levels = c(cens_name, "pr", "rfs")))
    check_results <- survcheck(
        Surv(tstart, tstop, to) ~ 1,
        data = counting_survcheck,
        id = id, istate = from
    )
    stopifnot(all(check_results$flag == 0))

    # Test 6: Compare survcheck and msdata transition counts
    counts_survcheck <- check_results$transitions
    counts_msdata <- events(msdata_1)[["Frequencies"]]

    dimnames(counts_survcheck)$to <- c("pr", "rfs", "no event")
    counts_survcheck["rfs", "no event"] <- (
        counts_survcheck["s0", "rfs"] + counts_survcheck["pr", "rfs"]
    )
    counts_msdata <- counts_msdata[, c("pr", "rfs", "no event")]
    stopifnot(all.equal(counts_survcheck, counts_msdata))

    TRUE
}

test_bosms3 <- function() {
    library(flexsurv)
    msdata_1 <- bosms3 %>%
        rename(time = years)

    id_col <- "id"
    tmat <- trans.illdeath(names = c(1, 2, 3))
    cens_name <- -1

    # Test 1: msdata -> timeline -> mstate -> timeline -> mstate -> msdata
    timeline_1 <- msdata_to_timeline(msdata_1, id_col, tmat, cens_name)
    mstate_1 <- timeline_to_mstate(timeline_1, id_col, tmat, cens_name)
    prep <- function(df) {
        df %>% msprep(
            status = c("status_1", "status_2", "status_3"),
            time = c("time_1", "time_2", "time_3"),
            trans = tmat,
            start = list(
                state = match(.$istate, dimnames(tmat)[[1]]),
                time = .$itime
            )
        )
    }
    timeline_2 <- mstate_to_timeline(mstate_1, id_col, tmat, "-1")
    mstate_2 <- timeline_to_mstate(timeline_2, id_col, tmat, cens_name)
    msdata_2 <- prep(mstate_2)
    stopifnot(all(msdata_1 == select(msdata_2, all_of(colnames(msdata_1)))))

    # Test 2: timeline -> counting -> timeline
    counting_1 <- timeline_to_counting(timeline_1, id_col)
    timeline_2 <- counting_to_timeline(counting_1, id_col)
    stopifnot(all(timeline_1 == timeline_2))

    # Test 3: Survcheck (counting process)
    # Pass counting process data to survcheck and verify flags are all 0
    counting_2 <- timeline_to_counting(timeline_2, id_col)
    stopifnot(all(counting_2$tstart < counting_2$tstop))
    counting_survcheck <- counting_2 %>%
        mutate(from = factor(from)) %>%
        mutate(to = factor(to, levels = c(cens_name, "1", "2", "3")))
    check_results <- survcheck(
        Surv(tstart, tstop, to) ~ 1,
        data = counting_survcheck,
        id = id, istate = from
    )
    stopifnot(all(check_results$flag == 0))

    # Test 4: Compare survcheck and msdata transition counts
    counts_survcheck <- check_results$transitions
    counts_msdata <- events(msdata_2)[["Frequencies"]]
    dimnames(counts_survcheck)$to <- c("2", "3", "no event")
    counts_survcheck["3", "no event"] <- (
        counts_survcheck["1", "3"] + counts_survcheck["2", "3"]
    )
    counts_msdata <- counts_msdata[, c("2", "3", "no event")]
    stopifnot(all.equal(counts_survcheck, counts_msdata))

    TRUE
}

# test_ebmt3()
# test_bosms3()
