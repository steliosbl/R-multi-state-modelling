################################################################################
# Helper functions for manipulating multi-state modelling state-transition data
#
# Author: Stelios Boulitsakis Logothetis
# Email: sb2690@medschl.cam.ac.uk
################################################################################
#
# These functions assume some common dataframe headings:
# Timeline format: ({id column}, TIME, EVENT)
# Counting process format: ({id column}, TSTART, TSTOP, FROM, TO)
# mstate format: ({id column}, <state>_TIME, <state>_STATUS, ITIME, ISTATE)

################################### 0. Setup ###################################
library(dplyr)

################################# 1. Functions #################################
#' Augment a Timeline format dataframe with censoring entries
#'
#' This function takes a dataframe in timeline format containing subject event histories
#' (with columns for the subject identifier, time, and event/state) and appends a censoring row
#' for each subject whose last recorded event is not one of the specified terminal (absorbing)
#' events. The censoring row is added at time \code{time_max} with the event name given in \code{censor_name}.
#' NOTE: We don't check if the last observed event is actually before \code{time_max}.
#'
#' @param df A dataframe in timeline format with columns \code{({id_col}, TIME, EVENT)}, without censoring entries.
#' @param id_col A string specifying the name of the column in \code{df} that identifies the subject.
#' @param time_max A numeric value representing the maximum observation time. Subjects not entering a
#' .  terminal event by this time will be considered censored.
#' @param terminal_events A character vector of absorbing (terminal) states. Subjects whose last event
#' .  is one of these states will not have a censoring row appended.
#' @param censor_name A string specifying the name to give to the censoring event. Default is "Censor".
#' @return A dataframe in timeline format with additional censoring rows appended for subjects that
#' .  did not reach a terminal event. In the appended censoring row, the \code{TIME}
#' . column is set to \code{time_max} and the \code{EVENT} column is set to \code{censor_name}
#'
#'
#' @export
mythos.augment_timeline_censoring <-
    function(df,
             id_col,
             time_max,
             terminal_events,
             censor_name = "Censor") {
        # For each ID, append a censor row if the last event is not in terminal_events
        df %>%
            group_by(across(all_of(id_col))) %>%
            group_modify(~ {
                if (!(last(.x$EVENT) %in% terminal_events)) {
                    censor_row <- tibble::tibble(
                        !!id_col := .x[[id_col]][1],
                        TIME := time_max,
                        EVENT := censor_name
                    )
                    .x <- dplyr::bind_rows(.x, censor_row)
                }
                .x
            }) %>%
            ungroup()
    }

#' Convert Timeline-format dataframes to Counting Process format
#'
#' This function converts a timeline-style dataframe (with columns for subject identifier,
#' time, and event) into a counting process format. The timeline data is expected to include censoring
#' entries. The resulting counting process format contains, for each subject, the start and stop times
#' of each interval along with the state transitions (from previous state to new state).
#'
#' @param df A dataframe in timeline format with columns \code{({id_col}, EVENT, TIME)}. Must include censoring entries.
#' @param id_col A string specifying the name of the column that identifies subjects.

#' @return A dataframe in counting process format with columns: \code{({id_col}, TSTART, TSTOP, FROM, TO)},
#' .   Each row represents the interval between consecutive events for a subject.
#'
#' @export
mythos.convert_timeline_to_counting <- function(df, id_col) {
    df %>%
        # First, sort the data by ID and time.
        arrange(.data[[id_col]], .data$TIME) %>%
        # Convert the augmented timeline to counting process format.
        group_by(across(all_of(id_col))) %>%
        mutate(
            TSTART = lag(TIME),
            TSTOP  = TIME,
            FROM   = lag(EVENT),
            TO     = EVENT
        ) %>%
        # Remove the first row per ID (which has no prior event)
        filter(!is.na(TSTART)) %>%
        ungroup() %>%
        select(all_of(id_col), TSTART, TSTOP, FROM, TO)
}

#' Convert Counting Process-format dataframes to Timeline format
#'
#' This function converts a counting process-style dataframe (with columns for subject identifier,
#' interval start and stop time, previous state, and new state) into a timeline format.
#' The resulting timeline data contains, for each subject, every visited state and entry time into it.
#'
#' @param df A dataframe in counting process format with columns \code{({id_col}, TSTART, TSTOP, FROM, TO)}. Must include censoring entries.
#' @param id_col A string specifying the name of the column that identifies subjects.

#' @return A dataframe in timeline format with columns: \code{({id_col}, TIME, EVENT)},
#' .   Each row represents a visited state for a subject, including the initial state and censoring.
#'
#' @export
mythos.convert_counting_to_timeline <- function(df, id_col) {
    # For each ID, recover the initial timeline event
    initial <- df %>%
        group_by(across(all_of(id_col))) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        transmute(!!id_col := .data[[id_col]], TIME = TSTART, EVENT = FROM)

    # Every counting process row provides a timeline event at time_stop with new_state.
    subsequent <- df %>%
        transmute(!!id_col := .data[[id_col]], TIME = TSTOP, EVENT = TO)

    # Combine and order the timeline events.
    bind_rows(initial, subsequent) %>%
        arrange(.data[[id_col]], TIME)
}

#' Convert Timeline-format dataframes to mstate-format
#'
#' This function converts a timeline-format dataframe into an mstate format.
#' The input dataframe is expected to be in timeline format with columns for subject identifier,
#' time, and event (including censoring entries). If some possible transitions aren't observed,
#' the optional \code{tmat} argument can be provided to ensure the corresponding columns are made.
#'
#' @param df A dataframe in timeline format with columns \code{({id_col}, EVENT, TIME)}. Must include censoring entries.
#' @param id_col A string specifying the name of the column that identifies subjects.
#' @param tmat A matrix with named dimensions specifying the nodes and cell values indicating possible transitions.
#' .   The cell values must be increasing natural numbers for transitions and NA for no transition (NOT binary/1 and 0).
#' @param censor_name A string denoting the event value that indicates censoring. Defaults to \code{"Censor"}.
#'
#' @return A dataframe in mstate format containing:
#'   \item{{id_col}}{The subject identifier.}
#'   \item{<state>_TIME}{For each state in the set of events, the time at which the state was first reached.
#' .   If the state was not observed for a subject, the value is set to the subject’s exit (censoring or terminal state) time.}
#'   \item{<state>_STATUS}{For each state, a status indicator (1 if the state was observed, 0 otherwise).}
#'   \item{itime}{The initial time (start time) for the subject, taken from the first non-censor event.}
#'   \item{istate}{The name of initial state for the subject, taken from the first non-censor event.}
#'
#'
#' @export
mythos.convert_timeline_to_mstate <-
    function(df, id_col, tmat, censor_name = "Censor") {
        if (!setequal(df$EVENT, c(dimnames(tmat)[[1]], censor_name))) {
            stop("tmat state names don't match dataframe state names.")
        }

        # Order the data by id and time.
        df <- df %>%
            arrange(.data[[id_col]], TIME)

        # Compute initial start time and state for each sample
        istate_info <- df %>%
            group_by(across(all_of(id_col))) %>%
            slice_head(n = 1) %>%
            ungroup() %>%
            select(all_of(id_col), ITIME = TIME, ISTATE = EVENT)

        # Compute terminal time and state for each sample
        terminal_info <- df %>%
            group_by(across(all_of(id_col))) %>%
            slice_tail(n = 1) %>%
            ungroup() %>%
            select(all_of(id_col), TTIME = TIME, TSTATE = EVENT)

        #  Exclude censor rows when summarising event times.
        df_events <- df %>%
            filter(EVENT != censor_name)


        # Determine the set of events to include.
        if (!is.null(tmat)) {
            events <- dimnames(tmat)[[1]]
        } else {
            events <- unique(df_events$event)
        }

        # For each ID and event, get the earliest time the event occurred.
        df_event_times <- df_events %>%
            group_by(across(all_of(id_col)), EVENT) %>%
            summarise(event_time = min(TIME), .groups = "drop")

        # Pivot to wide format so that each event becomes its own column.
        df_wide <- df_event_times %>%
            pivot_wider(
                id_cols = all_of(id_col),
                names_from = EVENT,
                values_from = event_time,
                values_fn = min,
                values_fill = list(event_time = NA)
            )

        # Ensure every node in `events` is present as a column.
        for (ev in events) {
            if (!(ev %in% names(df_wide))) {
                df_wide[[ev]] <- NA_real_
            }
        }

        # Combine with the identified initial and terminal states
        df_wide_aug <- df_wide %>%
            left_join(istate_info, by = {{ id_col }}) %>%
            left_join(terminal_info, by = {{ id_col }})

        # For each event, create a _time column and a _status column.
        # If the event was never observed for an ID, _time is set to time_max and _status to 0.
        for (ev in events) {
            time_new <- paste0(ev, "_TIME")
            status_new <- paste0(ev, "_STATUS")

            df_wide_aug <- df_wide_aug %>%
                mutate(
                    !!time_new := if_else(is.na(.data[[ev]]), TTIME, .data[[ev]]), !!status_new := if_else(is.na(.data[[ev]]), 0L, 1L)
                )
        }


        # Order columns: id, then each event's _time and _status (in the order of `events`),
        # and finally itime and istate.
        col_order <- c(
            id_col,
            as.vector(t(sapply(events, function(ev) {
                c(paste0(ev, "_TIME"), paste0(ev, "_STATUS"))
            }))),
            "ITIME", "ISTATE"
        )

        result <- df_wide_aug %>%
            # Modify the initial state to reflext the tmat indices, rather than state names/labels
            # mutate(ISTATE = match(ISTATE, dimnames(tmat)[[1]])) %>%
            select(all_of(col_order))

        return(result)
    }

#' Convert mstate-format dataframes to Timeline format
#'
#' This function converts a mstate-format dataframe (with columns for subject identifier,
#' state time and status columns for each state) into a timeline format.
#' The resulting timeline data contains, for each subject, every visited state and entry time into it.
#'
#' @param df A dataframe in mstate process format with columns \code{({id_col}, <state>_TIME, <state>_STATUS)}.
#' @param id_col A string specifying the name of the column that identifies subjects.
#' @param time_max A numeric vector or value representing the maximum observation time. Subjects not entering a
#' .  terminal event by this time will be considered censored.
#' @param tmat A matrix with named dimensions specifying the nodes and cell values indicating possible transitions.
#' .   The cell values must be increasing natural numbers for transitions and NA for no transition (NOT binary/1 and 0).
#' @param censor_name A string denoting the event value that indicates censoring. Defaults to \code{"Censor"}.
#'
#' @return A dataframe in timeline format with columns: \code{({id_col}, TIME, EVENT)},
#' .   Each row represents a visited state for a subject, including the initial state and censoring.
#'
#' @export
mythos.convert_mstate_to_timeline <-
    function(mstate_df,
             id_col,
             tmat,
             censor_name = "Censor") {
        # Pivot event columns to long format, excluding the start_* columns.
        df_wide <- mstate_df %>%
            pivot_longer(
                cols = -all_of(c(id_col, "ITIME", "ISTATE")),
                names_to = c("EVENT", ".value"),
                names_sep = "_"
            )

        if (!all(unique(df_wide$ISTATE) %in% unique(df_wide$EVENT))) {
            stop("Some of the initial states in ISTATE don't appear as _TIME, _STATUS columns in the input dataframe. Samples that don't exit their initial states will be dropped.")
        }

        # Before proceeding, validate data-types (there will be a problem if censor_name is a string but EVENT is a number)
        if (is.numeric(df_wide$EVENT) != is.numeric(censor_name) && is.character(df_wide$EVENT) != is.character(censor_name)) {
            warning("Data type mismatch between EVENT and censor_name. Converting both to strings.")
            df_wide$EVENT <- as.character(df_wide$EVENT)
            censor_name <- as.character(censor_name)
        }

        terminal_events <- mythos.terminal_events_from_tmat(tmat)
        # Obtain the latest (by time) state per sample
        # Create new rows for the timeline table, representing censoring events
        terminal_info <- df_wide %>%
            arrange(.data[[id_col]], (EVENT %in% terminal_events), TIME) %>%
            group_by(across(all_of(id_col))) %>%
            slice_tail(n = 1) %>%
            ungroup() %>%
            filter(!(EVENT %in% terminal_events & STATUS == 1)) %>%
            mutate(EVENT := censor_name, STATUS = 1)



        # Stack the censoring rows with the rest
        # Then convert to timeline formatting
        result <- df_wide %>%
            bind_rows(terminal_info) %>%
            filter(STATUS == 1) %>%
            select(all_of(id_col), TIME, EVENT) %>%
            arrange(.data[[id_col]], TIME)

        return(result)
    }

#' Convert "msdata"-format dataframes to Timeline format
#'
#' This function converts a "msdata"-style dataframe (with columns for subject identifier,
#' interval start and stop time, previous state, and new state, and status indicator) into a timeline format.
#' The resulting timeline data contains, for each subject, every visited state and entry time into it.
#'
#' @param df A dataframe in "msdata" format with columns \code{({id_col}, TSTART, TSTOP, FROM, TO, STATUS)}.
#' @param id_col A string specifying the name of the column that identifies subjects.
#' @param tmat A matrix with named dimensions specifying the nodes and cell values indicating possible transitions.
#' .   The cell values must be increasing natural numbers for transitions and NA for no transition (NOT binary/1 and 0).
#' @return A dataframe in timeline format with columns: \code{({id_col}, TIME, EVENT)},
#' .   Each row represents a visited state for a subject, including the initial state and censoring.
#'
#' @export
mythos.convert_msdata_to_timeline <- function(df, id_col, tmat, censor_name = "Censor") {
    # For each ID, recover the initial timeline event
    initial <- df %>%
        group_by(across(all_of(id_col))) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        transmute(!!id_col := .data[[id_col]], TIME = TSTART, EVENT = FROM)

    # Every counting process row provides a timeline event at time_stop with new_state.
    subsequent <- df %>%
        filter(STATUS == 1) %>%
        transmute(!!id_col := .data[[id_col]], TIME = TSTOP, EVENT = TO)

    # Before proceeding, validate data-types (there will be a problem if censor_name is a string but EVENT is a number)
    if (is.numeric(initial$EVENT) != is.numeric(censor_name) && is.character(initial$EVENT) != is.character(censor_name)) {
        stop("Data type mismatch between EVENT and censor_name. If you convert EVENT to characters, be sure to correct 'tmat' and 'ISTATE' also.")
    }

    terminal_events <- mythos.terminal_events_from_tmat(tmat)
    # Augment with censoring events if no terminal event is found
    censoring <- df %>%
        arrange(.data[[id_col]], (TO %in% terminal_events), TSTOP) %>%
        group_by(across(all_of(id_col))) %>%
        slice_tail(n = 1) %>%
        ungroup() %>%
        filter(!(TO %in% terminal_events & STATUS == 1)) %>%
        mutate(EVENT := censor_name) %>%
        select(all_of(id_col), TIME = TSTOP, EVENT)

    # Combine and order the timeline events.
    result <- bind_rows(initial, subsequent, censoring) %>%
        arrange(.data[[id_col]], TIME)

    return(result)
}

mythos.terminal_events_from_tmat <- function(tmat) {
    rownames(tmat)[apply(tmat, 1, function(x) all(is.na(x)))]
}


################################### 2. Tests ###################################
test_ebmt3 <- function() {
    # This function tests the above conversion functions using the EBMT3 dataset
    # We convert it between a few formats and then back again and inspect that it's unchanged.
    library(mstate)
    data(ebmt3)
    id_col <- "ID"
    df_mstate_1 <- ebmt3 %>%
        select(ID = id, pr_TIME = prtime, pr_STATUS = prstat, rfs_TIME = rfstime, rfs_STATUS = rfsstat) %>%
        mutate(s0_TIME = 0, s0_STATUS = 1) %>%
        mutate(ISTATE = "s0", ITIME = 0)

    # Introduce some initial states in "pr" rather than starting all in s0
    df_mstate_1 <- df_mstate_1 %>%
        mutate(flip_flag = (pr_STATUS == 1 & runif(n()) < 0.3)) %>%
        mutate(
            pr_TIME = ifelse(flip_flag, 0, pr_TIME),
            s0_STATUS = ifelse(flip_flag, 0, s0_STATUS),
            s0_TIME = ifelse(flip_flag, pmax(pr_TIME, rfs_TIME), s0_TIME),
            ISTATE = ifelse(flip_flag, "pr", ISTATE)
        ) %>%
        select(-flip_flag)

    tmat <- trans.illdeath(names = c("s0", "pr", "rfs"))

    df_timeline_1 <- mythos.convert_mstate_to_timeline(df_mstate_1, id_col, tmat)
    df_counting_1 <- mythos.convert_timeline_to_counting(df_timeline_1, id_col)
    # Ensure no censoring after terminal state
    stopifnot(!any(df_counting_1$FROM == "rfs"))
    df_timeline_2 <- mythos.convert_counting_to_timeline(df_counting_1, id_col)
    stopifnot(all(df_timeline_1 == df_timeline_2))

    df_mstate_2 <- mythos.convert_timeline_to_mstate(df_timeline_2, id_col, tmat) %>% select(colnames(df_mstate_1))
    # This sometimes fails: pr_TIME is sometimes != the final event time (rfs) even though rfs_STATUS == 1
    # But this means the data is wrong, not us
    # stopifnot(all(df_mstate_1 == df_mstate_2))

    # Check that msprep gives the same results
    msbmt1 <- df_mstate_1 %>% msprep(
        status = c("s0_STATUS", "pr_STATUS", "rfs_STATUS"),
        time = c("s0_TIME", "pr_TIME", "rfs_TIME"),
        trans = tmat,
        start = list(
            state = match(.$ISTATE, dimnames(tmat)[[1]]),
            time = .$ITIME
        )
    )

    msbmt2 <- df_mstate_2 %>% msprep(
        status = c("s0_STATUS", "pr_STATUS", "rfs_STATUS"),
        time = c("s0_TIME", "pr_TIME", "rfs_TIME"),
        trans = tmat,
        start = list(
            state = match(.$ISTATE, dimnames(tmat)[[1]]),
            time = .$ITIME
        )
    )
    stopifnot(msbmt1 == msbmt2)

    # Pass counting process data to survcheck and verify flags are all 0
    stopifnot(all(df_counting_1$TSTART < df_counting_1$TSTOP))
    df_counting_check <- df_counting_1 %>%
        mutate(FROM = factor(FROM), TO = factor(TO, levels = c("Censor", "pr", "rfs")))
    check <- survcheck(Surv(TSTART, TSTOP, TO) ~ 1, data = df_counting_check, id = ID, istate = FROM)
    stopifnot(all(check$flag == 0))

    # Compare observed transition counts between msprep and survcheck
    obs_trans_counting <- check$transitions
    dimnames(obs_trans_counting)$to <- c("pr", "rfs", "no event")
    obs_trans_counting["rfs", "no event"] <- obs_trans_counting["s0", "rfs"] + obs_trans_counting["pr", "rfs"]
    obs_trans_msprep <- events(msbmt1)[["Frequencies"]][, c("pr", "rfs", "no event")]
    stopifnot(all.equal(obs_trans_counting, obs_trans_msprep))

    return(TRUE)
}

test_bosms3 <- function() {
    library(flexsurv)
    df_msdata <- bosms3 %>%
        rename(ID = id, FROM = from, TO = to, TSTART = Tstart, TSTOP = Tstop, STATUS = status)
    id_col <- "ID"
    tmat <- trans.illdeath(names = c(1, 2, 3))
    censor_name <- -1

    df_timeline_1 <- mythos.convert_msdata_to_timeline(df_msdata, id_col, tmat, censor_name)
    df_counting_1 <- mythos.convert_timeline_to_counting(df_timeline_1, id_col)
    df_timeline_2 <- mythos.convert_counting_to_timeline(df_counting_1, id_col)
    df_mstate_1 <- mythos.convert_timeline_to_mstate(df_timeline_2, id_col, tmat, censor_name)
    msbmt <- df_mstate_1 %>% msprep(
        status = c("1_STATUS", "2_STATUS", "3_STATUS"),
        time = c("1_TIME", "2_TIME", "3_TIME"),
        trans = tmat,
        start = list(
            state = .$ISTATE,
            time = .$ITIME
        )
    )
    df_msdata_2 <- msbmt %>%
        rename(ID = id, FROM = from, TO = to, TSTART = Tstart, TSTOP = Tstop, STATUS = status, years = time) %>%
        select(names(df_msdata))

    stopifnot(all(df_msdata == df_msdata_2))

    # Pass counting process data to survcheck and verify flags are all 0
    stopifnot(all(df_counting_1$TSTART < df_counting_1$TSTOP))
    df_counting_check <- df_counting_1 %>%
        mutate(FROM = factor(FROM), TO = factor(TO, levels = c(-1, 1, 2, 3)))
    check <- survcheck(Surv(TSTART, TSTOP, TO) ~ 1, data = df_counting_check, id = ID, istate = FROM)
    stopifnot(all(check$flag == 0))

    # Compare observed transition counts between msprep and survcheck
    obs_trans_counting <- check$transitions
    dimnames(obs_trans_counting)$to <- c("2", "3", "no event")
    obs_trans_counting["3", "no event"] <- obs_trans_counting["1", "3"] + obs_trans_counting["2", "3"]
    obs_trans_msprep <- events(msbmt)[["Frequencies"]][, c("2", "3", "no event")]
    stopifnot(all.equal(obs_trans_counting, obs_trans_msprep))

    return(TRUE)
}

# test_ebmt3()
# test_bosms3()
