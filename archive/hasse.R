# This diagram is typically known as the Boolean lattice (or power set lattice).
# It represents the partially ordered set (poset) of all nonempty subsets of an S-element set
# (ordered by inclusion). In academic texts, you'll also find it referred to as the Hasse
# diagram of the Boolean lattice, often denoted as B_S or 2^[S].
# https://math.stackexchange.com/q/1622109/1153762
# The Hasse diagram of a partially ordered set P is the (directed) graph whose
# vertices are the elements of P and whose edges are the pairs (x, y) for which y covers x. It
# is usually drawn so that elements are placed higher than the elements they cover.
# FROM: https://pi.math.cornell.edu/~levine/18.312/alg-comb-lecture-7.pdf

library(igraph)
library(MASS)
library(Matrix)

## Part 1. Generate the Boolean lattice (Hasse diagram) with custom edge weights.
generate_boolean_lattice_old <- function(S) {
    # Create nodes: all non-empty subsets of {1,...,S}
    node_names <- c("0") # e.g., "1", "1,2", etc.
    levels <- c(0) # Level = size of the subset

    for (k in 1:S) {
        combs <- combn(1:S, k, simplify = FALSE)
        for (comb in combs) {
            node_names <- c(node_names, paste(sort(comb), collapse = ","))
            levels <- c(levels, k)
        }
    }

    vertices_df <- data.frame(
        id = 1:length(node_names),
        name = node_names,
        level = levels,
        stringsAsFactors = FALSE
    )

    # Create edges: from node A to node B if B is A plus one extra element.
    # For each node, assign outgoing edge weights that sum to p_total drawn uniformly from [0,1].
    edges <- data.frame(from = integer(), to = integer(), weight = numeric(), stringsAsFactors = FALSE)

    for (i in 1:nrow(vertices_df)) {
        # Get node i's subset and level
        current_level <- vertices_df$level[i]
        # Only consider edges from a node that is not in the last layer.
        if (current_level < S) {
            # Find all nodes that are exactly one level deeper.
            candidate_indices <- which(vertices_df$level == current_level + 1)

            if (current_level == 0) {
                valid_children <- candidate_indices
            } else {
                # A valid child must contain all elements in node i.
                A <- as.numeric(unlist(strsplit(vertices_df$name[i], split = ",")))
                valid_children <- c()
                for (j in candidate_indices) {
                    B <- as.numeric(unlist(strsplit(vertices_df$name[j], split = ",")))
                    if (all(A %in% B)) {
                        valid_children <- c(valid_children, j)
                    }
                }
            }

            # If there are valid children, assign weights.
            if (length(valid_children) > 0) {
                # Sample a total continuation probability in [0,1]
                p_total <- runif(1)
                # For the valid children, sample random numbers then normalize them so they sum to p_total.
                if (length(valid_children) == 1) {
                    weights <- p_total
                } else {
                    raw_weights <- runif(length(valid_children))
                    weights <- p_total * raw_weights / sum(raw_weights)
                }

                # Append an edge for each valid child with its corresponding weight.
                for (k in seq_along(valid_children)) {
                    edges <- rbind(edges, data.frame(
                        from = i,
                        to = valid_children[k],
                        weight = weights[k],
                        stringsAsFactors = FALSE
                    ))
                }
            }
        }
    }

    # Build and return the igraph object (directed from lower to higher layer).
    g <- graph_from_data_frame(d = edges, directed = TRUE, vertices = vertices_df)
    return(g)
}

## Part 2. Generate random walks down the lattice for each worker.
generate_lattice_walks <- function(g, num_workers, istate = c(0.75, 0.25)) {
    vertex_levels <- V(g)$level
    vertex_names <- V(g)$name

    walks <- list()

    for (worker in 1:num_workers) {
        # Start either at node 0 or at a random node in layer 1 (i.e., singletons).
        start_layer <- sample(c(0, 1), 1, prob = c(0.75, 0.25))
        layer1_nodes <- which(vertex_levels == start_layer)
        current <- layer1_nodes[sample(length(layer1_nodes), 1)]
        chain <- c(current)

        repeat {
            # Get valid children: out-neighbors at the next level.
            children <- as.numeric(neighbors(g, current, mode = "out"))
            valid_children <- children[vertex_levels[children] == vertex_levels[current] + 1]

            # If no children exist, break (terminal node).
            if (length(valid_children) == 0) break

            # Retrieve edge weights for transitions from current node.
            edge_ids <- sapply(valid_children, function(child) {
                get_edge_ids(g, vp = c(current, child))
            })
            edge_weights <- E(g)$weight[edge_ids]

            # The probability to continue is the sum of edge weights,
            # and the probability to stop is the remainder: 1 - sum(edge_weights).
            p_continue <- sum(edge_weights)
            p_stop <- 1 - p_continue

            # Form the probability vector for children and the "stop" option.
            choices <- c(valid_children, NA) # NA will indicate termination.
            probs <- c(edge_weights, p_stop)

            # Normalize probabilities in case of rounding issues.
            probs <- probs / sum(probs)

            # Sample the next action.
            next_step <- choices[sample(length(choices), 1, prob = probs)]

            if (is.na(next_step)) {
                # Terminate the walk.
                break
            } else {
                # Continue the walk.
                current <- next_step
                chain <- c(chain, current)
            }
        }

        # For each step in the chain, assign a random "time" ensuring monotonic increase.
        if (length(chain) == 1) {
            chain_times <- c(0)
        } else {
            chain_times <- sort(c(0, runif(length(chain) - 1)))
        }

        chain_names <- vertex_names[chain]

        chain_df <- data.frame(
            Worker_ID = rep(worker, length(chain)),
            time = chain_times,
            node = chain_names,
            stringsAsFactors = FALSE
        )

        walks[[worker]] <- chain_df
    }

    result_df <- do.call(rbind, walks)
    return(result_df)
}

plot_boolean_lattice <- function(g) {
    # Compute a layered layout using the Sugiyama algorithm.
    sug_layout <- layout_with_sugiyama(g)
    lay <- sug_layout$layout

    # Round edge weights for labeling.
    E(g)$label <- round(E(g)$weight, 2)

    # Plot the graph.
    plot(g,
        layout = lay,
        vertex.label = V(g)$name,
        vertex.size = 30,
        vertex.color = "skyblue",
        vertex.label.cex = 0.8,
        edge.arrow.size = 0.5,
        edge.label = E(g)$label,
        edge.label.cex = 0.8,
        edge.label.color = "black",
        edge.curved = 0.2, # Slight curvature helps separate labels
        edge.label.dist = 1.5, # Adjust distance of label from edge
        main = paste("Boolean Lattice (S =", S, ")"),
        margin = 0.2
    )
}

plot_boolean_lattice_2 <- function(g) {
    sug_layout <- layout_with_sugiyama(g)
    lay <- sug_layout$layout

    # Adjust plotting parameters based on S to avoid overcrowding.
    if (S <= 4) {
        vsize <- 30
        vcex <- 0.8
        ecex <- 0.8
        edge_label_dist <- 1.5
    } else {
        vsize <- max(15, 30 - (S - 4) * 3)
        vcex <- max(0.5, 0.8 - (S - 4) * 0.1)
        ecex <- max(0.5, 0.8 - (S - 4) * 0.1)
        edge_label_dist <- 2.0
    }

    E(g)$label <- round(E(g)$weight, 2)
    plot(g,
        layout = lay,
        vertex.label = V(g)$name, vertex.size = vsize, vertex.color = "skyblue",
        vertex.label.cex = vcex, edge.arrow.size = 0.5,
        edge.label = E(g)$label, edge.label.cex = ecex, edge.label.color = "black",
        edge.curved = 0.2, edge.label.dist = edge_label_dist,
        main = paste("Boolean Lattice (S =", S, ")"), margin = 0.2
    )
}


# Helper: Generate a skewed correlation matrix.
generate_skewed_corr <- function(C, shape1 = 1, shape2 = 5) {
    if (C == 1) {
        return(matrix(1, nrow = 1, ncol = 1))
    }
    R <- diag(1, C)
    for (i in 1:(C - 1)) {
        for (j in (i + 1):C) {
            # Sample from a Beta to get a value near 0 most of the time.
            r_val <- rbeta(1, shape1, shape2)
            # Randomly assign a sign.
            corr_val <- sample(c(-1, 1), 1) * r_val
            R[i, j] <- corr_val
            R[j, i] <- corr_val
        }
    }
    # Project onto the space of positive-definite matrices.
    as.matrix(nearPD(R)$mat)
}



generate_covariates <- function(g, df, C) {
    # Extract singleton node names (layer 1) from the graph.
    vertices <- as.data.frame(vertex_attr(g))
    layer1_nodes <- vertices$name[vertices$level == 1]

    # Randomly select one singleton per covariate.
    chosen_singletons <- sample(layer1_nodes, C, replace = TRUE)

    # Generate base distribution parameters for each covariate.
    absent_means <- runif(C, 0.2, 0.8)
    present_means <- absent_means + runif(C, 0.5, 1.0)
    absent_sds <- runif(C, 0.1, 0.3)
    present_sds <- runif(C, 0.1, 0.3)

    # Create a random correlation matrix.
    # R <- cor(matrix(rnorm(C * C), nrow = C))
    R <- generate_skewed_corr(C, shape1 = 1, shape2 = 5)


    workers <- unique(df$Worker_ID)
    n_workers <- length(workers)
    covariate_mat <- matrix(NA, nrow = n_workers, ncol = C)

    for (i in seq_along(workers)) {
        # Get the worker's walk as a vector of node labels.
        worker_walk <- df[df$Worker_ID == workers[i], "node"]

        # For each covariate, check if its chosen singleton is present anywhere in the walk.
        mean_vec <- sapply(1:C, function(j) {
            singleton <- chosen_singletons[j]
            # Use regex to detect the singleton as a distinct number within comma‐separated nodes.
            pattern <- paste0("(^|,)", singleton, "($|,)")
            if (any(grepl(pattern, worker_walk))) present_means[j] else absent_means[j]
        })
        sd_vec <- sapply(1:C, function(j) {
            singleton <- chosen_singletons[j]
            pattern <- paste0("(^|,)", singleton, "($|,)")
            if (any(grepl(pattern, worker_walk))) present_sds[j] else absent_sds[j]
        })

        # Build the covariance matrix from the standard deviations and the random correlation.
        # Ensure proper dimensions even if C==1
        Sigma <- diag(as.numeric(sd_vec), nrow = C, ncol = C) %*% R %*% diag(as.numeric(sd_vec), nrow = C, ncol = C)
        # Sample a correlated vector of covariate values for this worker.
        covariate_mat[i, ] <- as.numeric(mvrnorm(1, mu = mean_vec, Sigma = Sigma))
    }

    # Assemble the results into a dataframe.
    covariate_df <- data.frame(Worker_ID = workers, covariate_mat, stringsAsFactors = FALSE)
    colnames(covariate_df)[-1] <- paste0("C_", 1:C)

    # Attach the chosen singletons and distribution parameters as attributes for reference.
    attr(covariate_df, "chosen_singletons") <- chosen_singletons
    attr(covariate_df, "absent_means") <- absent_means
    attr(covariate_df, "present_means") <- present_means
    attr(covariate_df, "absent_sds") <- absent_sds
    attr(covariate_df, "present_sds") <- present_sds
    attr(covariate_df, "correlation_matrix") <- R

    return(covariate_df)
}

get_graph_tmat <- function(g) {
    n <- vcount(g)
    # Create an n x n matrix filled with NA, optionally with vertex names as row/column names
    mat <- matrix(NA, n, n, dimnames = list(V(g)$name, V(g)$name))

    # Get the edge list as numeric vertex indices
    edges <- as_edgelist(g, names = FALSE)

    # Use vectorized assignment to fill in increasing edge indices
    mat[cbind(edges[, 1], edges[, 2])] <- seq_len(nrow(edges))

    return(mat)
}

S <- 2
g <- generate_boolean_lattice(S)
df <- generate_lattice_walks(10, S)
covs <- generate_covariates(g, df, 5)
