library(cppRoutingCCH)
library(data.table)

parse_csv_arg <- function(value) {
  items <- strsplit(value, ",", fixed = TRUE)[[1]]
  items <- trimws(items)
  items[nzchar(items)]
}

read_benchmark_edges <- function(dataset) {
  if (dataset == "chicago") {
    file <- "data_readme/chicagoregional_net.csv"
    cost_col <- "length"
  } else if (dataset == "roads") {
    file <- "data_readme/roads.csv"
    cost_col <- "weight"
  } else {
    stop("unknown dataset: ", dataset)
  }

  raw <- data.table::fread(file, nThread = 1)
  edges <- data.frame(from = raw[[1]], to = raw[[2]], cost = raw[[cost_col]])
  edges[is.finite(edges$cost) & edges$cost >= 0, ]
}

make_benchmark_graph <- function(dataset, capacity = NULL, perturb = FALSE) {
  edges <- read_benchmark_edges(dataset)
  if (perturb) edges$cost <- edges$cost + seq_len(nrow(edges)) * 1e-9

  if (is.null(capacity)) {
    graph <- makegraph(edges, directed = TRUE)
  } else {
    graph <- makegraph(
      edges,
      directed = TRUE,
      alpha = 0.15,
      beta = 4,
      capacity = rep(capacity, nrow(edges))
    )
  }

  graph
}

make_random_od <- function(nodes, pairs_n, seed = 42L, demand = FALSE) {
  set.seed(seed)
  od <- data.frame(
    from = sample(nodes, pairs_n, replace = TRUE),
    to = sample(nodes, pairs_n, replace = TRUE)
  )
  if (demand) od$demand <- sample(1:10, pairs_n, replace = TRUE)
  od
}

avoid_same_node <- function(from, to, nodes) {
  same <- which(from == to)
  while (length(same) > 0) {
    to[same] <- sample(nodes, length(same), replace = TRUE)
    same <- same[from[same] == to[same]]
  }
  list(from = from, to = to)
}

make_shaped_od <- function(shape, nodes, pairs_n, seed) {
  set.seed(seed)
  if (shape == "sparse") {
    pairs <- avoid_same_node(
      sample(nodes, pairs_n, replace = TRUE),
      sample(nodes, pairs_n, replace = TRUE),
      nodes
    )
    demand <- sample(1:10, pairs_n, replace = TRUE)
  } else if (shape == "repeat_origins") {
    origin_n <- max(1L, min(25L, ceiling(sqrt(pairs_n))))
    origins <- sample(nodes, origin_n, replace = FALSE)
    pairs <- avoid_same_node(
      sample(origins, pairs_n, replace = TRUE),
      sample(nodes, pairs_n, replace = TRUE),
      nodes
    )
    demand <- sample(1:10, pairs_n, replace = TRUE)
  } else if (shape == "repeat_destinations") {
    dest_n <- max(1L, min(25L, ceiling(sqrt(pairs_n))))
    destinations <- sample(nodes, dest_n, replace = FALSE)
    pairs <- avoid_same_node(
      sample(nodes, pairs_n, replace = TRUE),
      sample(destinations, pairs_n, replace = TRUE),
      nodes
    )
    demand <- sample(1:10, pairs_n, replace = TRUE)
  } else if (shape == "volume_ranked") {
    candidate_n <- pairs_n * 3L
    pairs <- avoid_same_node(
      sample(nodes, candidate_n, replace = TRUE),
      sample(nodes, candidate_n, replace = TRUE),
      nodes
    )
    od <- data.table(from = pairs$from, to = pairs$to, demand = rlnorm(candidate_n, 0, 2))
    od <- od[, list(demand = sum(demand)), by = list(from, to)]
    setorder(od, -demand)
    od[, cum_share := cumsum(demand) / sum(demand)]
    od[, flow_rank := .I]
    od <- od[cum_share <= 0.95 | flow_rank == 1]
    if (nrow(od) > pairs_n) od <- od[seq_len(pairs_n)]
    return(od[, list(from, to, demand = demand / 0.95)])
  } else {
    stop("unknown OD shape: ", shape)
  }

  data.table(from = pairs$from, to = pairs$to, demand = as.numeric(demand))
}

compare_dist <- function(reference, candidate, tolerance = 1e-7) {
  both_inf <- is.infinite(reference) & is.infinite(candidate)
  diff <- abs(reference - candidate)
  diff[both_inf] <- 0
  finite_diff <- diff[is.finite(diff)]
  c(
    mismatches = sum(!(both_inf | diff <= tolerance), na.rm = TRUE),
    max_abs_diff = if (length(finite_diff)) max(finite_diff) else 0
  )
}

make_assignment_like_weights <- function(base_cost, iteration) {
  n <- length(base_cost)
  edge_pressure <- ((seq_len(n) * 1103515245 + iteration * 12345) %% 1000) / 1000
  base_cost * (1 + 0.15 * (0.25 + 1.75 * edge_pressure)^4)
}

sort_flow_table <- function(x) {
  x[order(x$from, x$to), ]
}
