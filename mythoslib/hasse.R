# HEADER #######################################################################
#
# Author: Stelios Boulitsakis Logothetis
# Email: sb2690@medschl.cam.ac.uk
#
# Date: 2025-02-26
# Name:	mythoslib/hasse.R
# Description: Helper functions for generating health state-transition graphs
#
# SETUP ########################################################################
# Load required libraries
library(igraph)
library(MASS)
library(Matrix)
library(dplyr)

source("mythoslib/convert.R")

# Fix MASS conflict
select <- dplyr::select

# FUNCTIONS ####################################################################
#' Generate the Hasse diagram of the power set P({1,...,n})
#'
#' The Hasse diagram of a partially ordered set P is the directed acyclic graph
#' whose vertices are the elements of P and whose edges are the pairs (x, y)
#' for which x is a subset of y. This function generates the Hasse diagram of
#' the power set of the set of integers {1, ..., n} given size parameter n.
#' Thus, there is a vertex for each combination of {1, ..., n}, and an edge
#' between two vertices if one is a subset of the other.
#'
#' The edges are weighted to simulate state-transition probabilities.
#' Each vertex is assigned a probability of staying in the same state, Pr(stay),
#' and probabilities of transitioning to each of its children in the graph.
#' The sum of these probabilities is 1.
#'
#' @param n The size of the set {1, ..., n} to generate the lattice for
#' @param pr_stay_max The maximum probability of staying in the same state
#' @return An igraph object representing the Hasse diagram
#' of the power set P({1,...,n})
poset_hasse <- function(n, pr_stay_max = 0.25) {
    # Create vertices: all non-empty subsets of Poset({1,...,n})
    names <- c("0") # e.g., "0", "1", "2", "1,2"
    levels <- c(0) # e.g., 0 1 1 1 1 2 2 2 2 2 2 3 3 3 3 4

    # Generate nodes and levels they belong to
    for (i in 1:n) {
        combinations <- combn(1:n, i, simplify = FALSE)
        for (c in combinations) {
            names <- c(names, paste(sort(c), collapse = ","))
            levels <- c(levels, i)
        }
    }
    vertices <- data.frame(
        id = seq_along(names),
        name = names,
        level = levels,
        stringsAsFactors = FALSE
    )

    # Create edges: if a vertex X is a subset of vertex Y
    # i.e., if Y contains X plus one additional element
    # For each vertex, assign outgoing edge weights that sum to Pr(out) ~ [0,1]
    edge_dfs <- list()
    for (v in seq_along(names)) {
        vertex <- names[v]
        level <- levels[v]
        vertex_contents <- as.numeric(unlist(strsplit(vertex, ",")))

        # If this is the final level, no outgoing edges
        if (level >= n) {
            next
        }

        # Find all nodes that are exactly one level deeper
        candidate_children <- which(levels == level + 1)
        children <- c()

        # All nodes in level 0 are children of level 0 (the root)
        # For levels > 0, we pick the vertices that are supersets of the current
        if (level == 0) {
            children <- candidate_children
        } else {
            # Within this level, find the vertices that contain the current one
            for (c in candidate_children) {
                child_contents <- as.numeric(unlist(strsplit(names[c], ",")))
                # If all the current vertex contents are in the candidate
                if (all(vertex_contents %in% child_contents)) {
                    # Then it is a neighbour
                    children <- c(children, c)
                }
            }
        }

        # Assign random but valid probability weights to outgoing edges
        # First sample the probability of leaving this vertex at all, Pr(out)
        pr_out <- runif(1, min = 1 - pr_stay_max, max = 1)

        # Then, for the valid children, all the weights must sum to Pr(out)
        if (length(children) == 1) {
            weights <- pr_out
        } else {
            weights <- runif(length(children))
            weights <- pr_out * weights / sum(weights)
        }

        edge_df <- data.frame(
            from = rep(v, length(children)),
            to = children,
            weight = weights,
            stringsAsFactors = FALSE
        )
        edge_dfs[[v]] <- edge_df
    }
    edges <- do.call(rbind, edge_dfs)

    graph_from_data_frame(edges, vertices = vertices, directed = TRUE)
}

#' Plot the Hasse diagram of a power set
#' @param graph An igraph object representing the Hasse diagram
#' @return A plot of the Hasse diagram
plot.hasse <- function(graph) {
    sug_layout <- layout_with_sugiyama(graph)$layout

    plot(graph,
        layout = sug_layout, vertex.label = V(graph)$name,
        vertex.size = 10, vertex.label.cex = 0.8,
        edge.label = round(E(graph)$weight, 2),
        edge.label.cex = 0.8, edge.arrow.size = 0.5,
        edge.curved = 0, edge.width = 2
    )
}

#' Extract the transition matrix from the Hasse diagram
#' @param graph An igraph object representing the Hasse diagram
#' @return The transition matrix of the Hasse diagram
hasse_tmat <- function(graph) {
    n <- vcount(graph)

    # Initialise an n*n null matrix with vertex names as row/column labels
    result <- matrix(NA, n, n, dimnames = list(V(graph)$name, V(graph)$name))

    # Get the edge list as numeric vertex indices
    edges <- as_edgelist(graph, names = FALSE)

    # Use vectorised assignment to fill in the transition matrix
    # with increasing edge indices (as required by mstate)
    result[cbind(edges[, 1], edges[, 2])] <- seq_len(nrow(edges))

    result
}

#' Generate n random walks on the Hasse diagram of a power set P({1,...,s})
#'
#' This function generates n random walks on the Hasse diagram of the power set
#' P({1,...,s}). The walks are generated by starting at the root of the Hasse
#' diagram, i.e. the empty set. We also permit walks to start at level 1 (i.e.
#' the singletons) with probability proba_layer1 The walks then proceed by
#' randomly selecting a child node to transition to, with probabilities
#' determined by the edge weights. At each step, there is also a probability
#' Pr(stay) of staying in the same state. The walks continue until either
#' drawing a "stay" or reaching the final level of the Hasse diagram.
#' Each step is also assigned a monotonically increasing time.
#'
#' @param graph An igraph object representing the Hasse diagram
#' @param n The number of walks to generate
#' @param proba_layer1 The probability of starting at level 1 (the singletons)
#'
#' @return A list of n random walks on the Hasse diagram. Given as a dataframe
#' .   \code{id, vertex, time}.
hasse_walks <- function(graph, n, proba_layer1 = 0.25) {
    names <- V(graph)$name
    levels <- V(graph)$level

    walks <- list()
    for (i in 1:n) {
        # Start either at the root or a random vertex in layer 1 with
        # with probability given by proba_layer1.
        layer <- sample(c(0, 1), 1,
            prob = c(1 - proba_layer1, proba_layer1)
        )
        layer_vertices <- which(levels == layer)
        current <- layer_vertices[sample(length(layer_vertices), 1)]
        chain <- c(current)

        repeat {
            # Get valid children of the current vertex
            children <- as.numeric(neighbors(graph, current, mode = "out"))

            if (length(children) == 0) break # No children, end of walk

            # Get edge weights for outgoing edges
            weights <- E(graph)[.from(current)]$weight

            # The probability to stop here is 1- the sum of the weights
            p_stop <- 1 - sum(weights)

            # Sample which step to take next - either a child vertex or staying
            choices <- c(children, NA)
            next_step <- choices[
                sample(length(choices), 1, prob = c(weights, p_stop))
            ]

            if (is.na(next_step)) {
                break # Stop here
            } else {
                chain <- c(chain, next_step)
                current <- next_step
            }
        }

        # For each step in the chain, assign a monotonically increasing time
        if (length(chain) == 1) {
            times <- c(0)
        } else {
            times <- sort(c(0, runif(length(chain) - 1)))
        }

        chain_df <- data.frame(
            id = rep(i, length(chain)),
            vertex = names[chain],
            time = times,
            stringsAsFactors = FALSE
        )
        walks[[i]] <- chain_df
    }

    do.call(rbind, walks)
}

#' Generate a n*n correlation matrix that's skewed towards 0
#' @param n The size of the correlation matrix
#' @param alpha The alpha shape parameter of the Beta distribution
#' @param beta The beta shape parameter of the Beta distribution
#' @return A n*n random correlation matrix with skew towards 0
skewed_corr <- function(n, alpha = 1, beta = 5) {
    if (n == 1) {
        return(matrix(1, nrow = 1, ncol = 1))
    }

    r <- diag(1, n)
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            # Sample from a Beta to get a value near 0 most of the time.
            r_val <- rbeta(1, alpha, beta)
            # Randomly assign a sign.
            corr_val <- sample(c(-1, 1), 1) * r_val
            r[i, j] <- corr_val
            r[j, i] <- corr_val
        }
    }

    # Project onto the space of positive-definite matrices.
    as.matrix(nearPD(r)$mat)
}

#' Generate synthetic covariates for the simulated walks on the Hasse diagram
#'
#' Given a Hasse diagram and a set of random walks on it, this function
#' generates n synthetic covariates for the walks. The covariates are sampled
#' from Gaussian distributions and linked via a random correlation matrix that's
#' skewed towards 0. Further, each covariate is dependent on a singleton
#' element in the Poset. If that element is present in the walk, the covariate
#' is sampled from a different (higher) distribution than if it is not present.
#'
#' @param graph An igraph object representing the Hasse diagram
#' @param walks A dataframe of random walks on the Hasse diagram with columns
#' .   \code{id, vertex, time}
#' @param n The number of covariates to generate
#' @return A dataframe of synthetic covariates for the walks with columns:
#' .   \code{id, C_1, C_2, ..., C_n}
hasse_covariates <- function(graph, walks, n) {
    # Identify the singleton nodes (layer 1) in the Hasse graph
    # and randomly select one for each covariate to depend on
    singletons <- V(graph)$name[
        which(V(graph)$level == 1) %>% sample(n, replace = TRUE)
    ]

    # Generate base distribution parameters for each covariate
    # When the singleton is present in the walk, the covariate tends
    # to take higher values than if it was absent
    absent_means <- runif(n, 0.2, 0.8)
    present_means <- absent_means + runif(n, 0.5, 1.0)
    absent_sds <- runif(n, 0.1, 0.3)
    present_sds <- absent_sds + runif(n, 0.1, 0.3)

    # Create a random correlation matrix skewed towards 0
    corr <- skewed_corr(n)

    # Pick the latest step in each walk as we know that vertex will contain
    # all the singletons that are present in the walk
    walks_max <- walks %>%
        group_by(id) %>%
        slice_max(time, n = 1, with_ties = FALSE) %>%
        ungroup()


    result <- apply(walks_max, 1, function(walk) {
        singleton_mask <- sapply(singletons, function(s) {
            any(grepl(paste0("\\b", s, "\\b"), walk[["vertex"]]))
        })

        mean_vec <- ifelse(singleton_mask, present_means, absent_means)
        sd_vec <- ifelse(singleton_mask, present_sds, absent_sds)

        # Build the covariance matrix from the standard deviations
        # and the generated correlation matrix
        sd_mat <- diag(as.numeric(sd_vec), nrow = n, ncol = n)
        sigma <- sd_mat %*% corr %*% sd_mat

        # Sample a correlated vector of covariate values for this walk
        as.numeric(
            MASS::mvrnorm(1, mean_vec, sigma)
        )
    })

    result_df <- data.frame(
        id = walks_max$id,
        t(result),
        stringsAsFactors = FALSE
    )
    colnames(result_df)[-1] <- paste0("C_", 1:n)

    # Attach the chosen singletons and distribution parameters
    # as attributes for reference
    attr(result_df, "singletons") <- singletons
    attr(result_df, "absent_means") <- absent_means
    attr(result_df, "present_means") <- present_means
    attr(result_df, "absent_sds") <- absent_sds
    attr(result_df, "present_sds") <- present_sds
    attr(result_df, "correlation_matrix") <- corr

    result_df
}

#' Generate longitudinal covariates for the simulated walks on the Hasse diagram
#'
#' Given a Hasse diagram and a set of random walks on it, this function
#' generates n synthetic covariates for the walks. The covariates are sampled
#' from Gaussian distributions and linked via a random correlation matrix that's
#' skewed towards 0. Further, each covariate is dependent on a singleton
#' element in the Poset. If that element is present in the walk, the covariate
#' is sampled from a different (higher) distribution than if it is not present.
#' We generate longitudinal covariates, meaning one value per step in the walk.
#' The values for a given ID are correlated across steps. A covariate may be
#' constant, increasing, or decreasing, and the trends upward or downward follow
#' and exponential curve.
#'
#' @param graph An igraph object representing the Hasse diagram
#' @param walks A dataframe of random walks on the Hasse diagram with columns
#' .   \code{id, vertex, time}
#' @param n The number of covariates to generate
#' @return A dataframe of synthetic covariates for the walks with columns:
#' .   \code{id, vertex, time, C_1, C_2, ..., C_n}
hasse_covariates_longitudinal <- function(graph, walks, n) {
    # Identify the singleton nodes (layer 1) in the Hasse graph
    # and randomly select one for each covariate to depend on
    singletons <- V(graph)$name[
        which(V(graph)$level == 1) %>% sample(n, replace = TRUE)
    ]

    # Generate base distribution parameters for each covariate
    # When the singleton is present in the walk, the covariate tends
    # to take higher values than if it was absent
    absent_means <- runif(n, 0.2, 0.8)
    present_means <- absent_means + runif(n, 0.5, 1.0)
    absent_sds <- runif(n, 0.1, 0.3)
    present_sds <- absent_sds + runif(n, 0.1, 0.3)

    # Generate trend types for covariates. Each one may be constant, increasing,
    # or decreasing, depending on if its dependant singleton is present
    trend_present <- sample(c(0, +1, -1), n, replace = TRUE)
    trend_absent <- sample(c(0, +1, -1), n, replace = TRUE)

    # Create a random correlation matrix skewed towards 0
    corr <- skewed_corr(n)

    # Pick the latest step in each walk as we know that vertex will contain
    # all the singletons that are present in the walk
    walks_max <- walks %>%
        group_by(id) %>%
        slice_max(time, n = 1, with_ties = FALSE) %>%
        ungroup()

    # Generate longitudinal series per step in each walk
    long <- group_by(walks, id) %>%
        group_map(function(walk, id) {
            id <- id$id
            n_steps <- nrow(walk)

            exit <- slice_max(walk, time, n = 1, with_ties = FALSE)$vertex
            singleton_mask <- sapply(singletons, function(s) {
                any(grepl(paste0("\\b", s, "\\b"), exit))
            })
            # For each covariate, identify the appropriate mean, std, and trend
            mean_vec <- ifelse(singleton_mask, present_means, absent_means)
            sd_vec <- ifelse(singleton_mask, present_sds, absent_sds)
            trends <- ifelse(singleton_mask, trend_present, trend_absent)

            # Build the covariance matrix from the standard deviations
            # and the generated correlation matrix
            sd_mat <- diag(as.numeric(sd_vec), nrow = n, ncol = n)
            sigma <- sd_mat %*% corr %*% sd_mat

            # Sample a correlated vector of covariate values for this walk
            baseline <- as.numeric(
                MASS::mvrnorm(1, mean_vec, sigma)
            )

            # Generate a time series for each covariate
            covs <- sapply(1:n, function(c) {
                # Generate a small offset delta for the trend
                delta <- runif(1, 0, sd_vec[c])

                # Generate a linear trend from the baseline to the final
                # value over the steps. The size of each step depends on
                # the time gap between the successive steps
                if (n_steps == 1) {
                    w <- w_exp <- c(0)
                } else {
                    w <- (walk$time - walk$time[1]) /
                        (walk$time[n_steps] - walk$time[1])
                    k <- 3 # exponential curvature parameter
                    w_exp <- (exp(k * w) - 1) / (exp(k) - 1)
                }

                # The value at each time step is the baseline,
                # plus the delta in the correct direction (0 if constant trend)
                # weighted by the time difference between steps
                baseline[c] + (trends[c]) * delta * (1 - w_exp)
            }) %>%
                matrix(nrow = n_steps, ncol = n) %>%
                as.data.frame()
        }) %>%
        bind_rows()

    # Name the covariate columns as "C_1", "C_2", ..., "C_n".
    colnames(long) <- paste0("C_", 1:n)

    result_df <- walks %>% bind_cols(long)
    # Attach attributes for reference.
    attr(result_df, "singletons") <- singletons
    attr(result_df, "absent_means") <- absent_means
    attr(result_df, "present_means") <- present_means
    attr(result_df, "absent_sds") <- absent_sds
    attr(result_df, "present_sds") <- present_sds
    attr(result_df, "correlation_matrix") <- corr
    attr(result_df, "trend_present") <- trend_present
    attr(result_df, "trend_absent") <- trend_absent

    result_df
}

#' Convert walks and covariates along the Hasse graph to counting process data
#'
#' Given a set of random walks on the Hasse diagram and a set of covariates
#' generated for these walks, this function converts the data to a counting
#' process format suitable for use with the survival package.
#'
#' @param walks A dataframe of random walks on the Hasse diagram with columns
#' .   \code{id, vertex, time}
#' @param covariates A dataframe of covariates for the walks with columns:
#' .   \code{id, C_1, C_2, ..., C_n}
#' @param tmat The transition matrix of the Hasse diagram
#' @return A dataframe of counting process data with columns:
#' .   \code{id, start, stop, event, C_1, C_2, ..., C_n}
#' .   And attribute "covariates" giving the covariate column names
convert_hasse_counting <- function(walks, covariates, tmat, time_max = 1.0) {
    result <- walks %>%
        rename(event = vertex) %>%
        # Convert to a timeline
        augment_timeline_censoring("id", tmat, time_max, "censor") %>%
        # Convert to counting process
        timeline_to_counting(., "id") %>%
        # Event columns to factors with censoring as the 1st level
        mutate(from = factor(from), to = factor(to)) %>%
        mutate(to = fct_relevel(as.factor(to), "censor")) %>%
        # Augment with (time-fixed/baseline) covariate columns
        left_join(covariates, by = "id")

    # Include covariate column names as an attribute
    cov_names <- colnames(covariates)[colnames(covariates) != "id"]
    attr(result, "covariates") <- cov_names

    result
}

#' Convert walks and covariates along the Hasse graph to mstate data
#'
#' Given a set of random walks on the Hasse diagram and a set of covariates
#' generated for these walks, this function converts the data to a multi-state
#' format suitable for use with the mstate package. Includes the call to msprep.
#'
#' @param walks A dataframe of random walks on the Hasse diagram with columns
#' .   \code{id, vertex, time}
#' @param covariates A dataframe of covariates for the walks with columns:
#' .   \code{id, C_1, C_2, ..., C_n}
#' @param tmat The transition matrix of the Hasse diagram
#' @param time_max The maximum time to consider for each walk
#' @param expand Whether to expand the covariates to transition-specific columns
#' @return A msdata object with columns:
#' .   \code{id, from, to, Tstart, Tstop, time, status, C_1, C_2, ..., C_n,
#'     C_1.1, C_2.1, C_3.1, ... C_S.1, C_2.2, ..., C_S.n}
#' .   Includes attributes "covariates" giving the covariate column names and
#' .   "covariates_expanded" giving for transition-specific covariate names
convert_hasse_msdata <- function(walks, covariates, tmat, time_max = 1.0, expand = TRUE) {
    cov_names <- colnames(covariates)[colnames(covariates) != "id"]

    df_mstate <- walks %>%
        # Convert to a timeline
        rename(event = vertex) %>%
        augment_timeline_censoring("id", tmat, time_max, "censor") %>%
        # Convert to mstate format
        timeline_to_mstate("id", tmat, "censor") %>%
        select(-time_censor, -status_censor) %>%
        # Augment with (time-fixed/baseline) covariate columns
        left_join(covariates, by = "id")

    # Call msprep to create msdata object
    df_msdata <- df_mstate %>% msprep(
        time = paste0("time_", rownames(tmat)),
        status = paste0("status_", rownames(tmat)),
        trans = tmat,
        id = "id",
        start = list(
            time = .$itime,
            state = match(.$istate, dimnames(tmat)[[1]])
        ),
        keep = cov_names
    )

    if (!expand) {
        attr(df_msdata, "covariates") <- cov_names
        return(df_msdata)
    }
    
    # Create the table of expanded covariates (do not join with states yet)
    df_covs_expanded <- df_msdata %>%
        expand.covs(cov_names, append = FALSE, longnames = TRUE)

    # Join with states
    result <- bind_cols(df_msdata, df_covs_expanded)

    # Include covariate column names as attributes
    attr(result, "covariates") <- cov_names
    attr(result, "covariates_expanded") <- colnames(df_covs_expanded)

    result
}
