library(cppRoutingCCH)
library(data.table)

parse_csv_arg <- function(value) {
  items <- strsplit(value, ",", fixed = TRUE)[[1]]
  items <- trimws(items)
  items[nzchar(items)]
}

decompress_gzip <- function(source, target) {
  input <- gzfile(source, "rb")
  output <- file(target, "wb")
  on.exit(close(input))
  on.exit(close(output), add = TRUE)

  repeat {
    chunk <- readBin(input, what = "raw", n = 1024L * 1024L)
    if (!length(chunk)) break
    writeBin(chunk, output)
  }

  invisible(target)
}

read_benchmark_edges <- function(dataset) {
  dimacs_file <- dataset
  if (file.exists(dataset) && grepl("[.]gr[.]gz$", dataset, ignore.case = TRUE)) {
    dimacs_file <- tempfile(fileext = ".gr")
    on.exit(unlink(dimacs_file), add = TRUE)
    decompress_gzip(dataset, dimacs_file)
  }

  if (file.exists(dimacs_file) && grepl("[.]gr$", dimacs_file, ignore.case = TRUE)) {
    raw <- data.table::fread(dimacs_file, skip = "a ", header = FALSE, nThread = 1)
    if (ncol(raw) != 4L || !all(raw[[1]] == "a")) {
      stop("invalid DIMACS .gr file: ", dataset)
    }
    edges <- data.frame(from = raw[[2]], to = raw[[3]], cost = raw[[4]])
    return(edges[is.finite(edges$cost) & edges$cost >= 0, ])
  } else if (dataset == "chicago") {
    file <- "data_readme/chicagoregional_net.csv"
    cost_col <- "length"
  } else if (dataset == "roads") {
    file <- "data_readme/roads.csv"
    cost_col <- "weight"
  } else {
    stop("unknown dataset or DIMACS .gr path: ", dataset)
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
  pairs <- avoid_same_node(
    sample(nodes, pairs_n, replace = TRUE),
    sample(nodes, pairs_n, replace = TRUE),
    nodes
  )
  od <- data.frame(from = pairs$from, to = pairs$to)
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

is_prepared_cch_method <- function(method) {
  identical(method, "cch_prepared")
}

cch_aon_method <- function(method) {
  if (is_prepared_cch_method(method)) "cch" else method
}

flow_check <- function(reference, candidate, tolerance = 1e-7) {
  reference <- sort_flow_table(reference)
  candidate <- sort_flow_table(candidate)
  flow_diff <- abs(reference$flow - candidate$flow)
  c(
    flow_mismatches = sum(flow_diff > tolerance),
    max_flow_diff = max(flow_diff),
    max_cost_diff = max(abs(reference$cost - candidate$cost))
  )
}
