# Install and load igraph if not already installed
library(igraph)

## Part 1. Generate the Boolean lattice (Hasse diagram)
generate_boolean_lattice <- function(S) {
    node_names <- c()  # Labels like "1", "1,2", etc.
    levels <- c()      # Level indicates the size of the subset
    
    # For each layer k (from singletons up to the full set)
    for (k in 1:S) {
        combs <- combn(1:S, k, simplify = FALSE)
        for (comb in combs) {
            node_names <- c(node_names, paste(sort(comb), collapse = ","))
            levels <- c(levels, k)
        }
    }
    
    vertices_df <- data.frame(id = 1:length(node_names),
                              name = node_names,
                              level = levels,
                              stringsAsFactors = FALSE)
    
    # Create edges: from node A to node B if B is A with one extra element.
    edges <- data.frame(from = integer(), to = integer())
    for (i in 1:nrow(vertices_df)) {
        A <- as.numeric(unlist(strsplit(vertices_df$name[i], split = ",")))
        for (j in 1:nrow(vertices_df)) {
            if (vertices_df$level[j] == vertices_df$level[i] + 1) {
                B <- as.numeric(unlist(strsplit(vertices_df$name[j], split = ",")))
                if (all(A %in% B)) {
                    edges <- rbind(edges, data.frame(from = vertices_df$id[i], to = vertices_df$id[j]))
                }
            }
        }
    }
    
    # Build and return the igraph object (directed from lower to higher layer)
    g <- graph_from_data_frame(d = edges, directed = TRUE, vertices = vertices_df)
    return(g)
}

## Part 2. Generate random walks down the lattice for each worker.
generate_lattice_walks <- function(num_workers, S) {
    # First, create the Boolean lattice for base set of size S.
    g <- generate_boolean_lattice(S)
    vertex_levels <- V(g)$level
    vertex_names <- V(g)$name
    
    walks <- list()
    
    for (worker in 1:num_workers) {
        # Decide a random chain length (from 1 up to S)
        L <- sample(1:S, 1)
        
        # Start at a random node in layer 1 (i.e. singletons)
        layer1_nodes <- which(vertex_levels == 1)
        current <- sample(layer1_nodes, 1)
        chain <- c(current)
        
        # Traverse downward: at each step choose a child (a superset with one extra element)
        for (i in 2:L) {
            # Get out-neighbors (children) of the current node
            children <- as.numeric(neighbors(g, current, mode = "out"))
            # Filter for those exactly one level deeper
            valid_children <- children[ vertex_levels[children] == vertex_levels[current] + 1 ]
            if (length(valid_children) == 0) break  # Stop if no children exist
            current <- valid_children[sample(length(valid_children), 1)] # Modified from sample(valid_children) because it bugs out at the last layer when the length is 1
            chain <- c(chain, current)
        }
        
        # For each step, assign a random "time" ensuring monotonic increase.
        # Generate random times and then sort them so that they increase along the walk.
        chain_times <- sort(runif(length(chain)))
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

# Example usage:
set.seed(123)    # For reproducibility
num_workers <- 500 # Number of workers
S <- 4           # Parameter: size of the base set (number of base elements)
lattice_walks_df <- generate_lattice_walks(num_workers, S)
print(lattice_walks_df)

"1,2,3,4" %in% lattice_walks_df$node
