source("benchmarks/bench_utils.R")

RcppParallel::setThreadOptions(numThreads = 1)

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[[1]] else "chicago"
iterations <- if (length(args) >= 2) as.integer(args[[2]]) else 3L
pairs_n <- if (length(args) >= 3) as.integer(args[[3]]) else if (dataset == "roads") 100L else 1000L
matrix_n <- if (length(args) >= 4) as.integer(args[[4]]) else 25L
seed <- if (length(args) >= 5) as.integer(args[[5]]) else 42L

graph <- make_benchmark_graph(dataset, perturb = TRUE)
base_cost <- graph$data$dist
od <- make_random_od(graph$dict$ref, pairs_n, seed = seed)
matrix_from <- unique(od$from)[seq_len(min(matrix_n, length(unique(od$from))))]
matrix_to <- unique(od$to)[seq_len(min(matrix_n, length(unique(od$to))))]

t_ch_prepare <- system.time(ch <- cpp_contract(graph, silent = TRUE))
t_cch_prepare <- system.time(cch <- cpp_cch_prepare(graph))

cat("DATASET", dataset, "\n")
cat("nodes=", graph$nbnode, "input_edges=", nrow(graph$data),
    "ch_edges=", nrow(ch$data), "cch_arcs=", length(cch$tail), "\n")
print(data.frame(
  phase = c("ch_prepare", "cch_prepare"),
  elapsed_sec = c(t_ch_prepare[["elapsed"]], t_cch_prepare[["elapsed"]])
))

timings <- data.table()
checks <- data.table()

distance_check_row <- function(iteration, comparison, reference, candidate) {
  check <- compare_dist(reference, candidate)
  data.table(
    iteration,
    comparison,
    mismatches = unname(check[["mismatches"]]),
    max_abs_diff = unname(check[["max_abs_diff"]])
  )
}

for (iteration in seq_len(iterations)) {
  weights <- make_assignment_like_weights(base_cost, iteration)
  graph_iter <- graph
  graph_iter$data$dist <- weights

  t_customize <- system.time(metric <- cpp_cch_customize(cch, weights))
  t_cch_pair <- system.time(cch_pair <- get_distance_pair(metric, od$from, od$to))
  t_dijkstra <- system.time(reference_pair <- get_distance_pair(
    graph_iter, od$from, od$to, algorithm = "Dijkstra"
  ))
  t_ch_contract <- system.time(ch_iter <- cpp_contract(graph_iter, silent = TRUE))
  t_ch_pair <- system.time(ch_pair <- get_distance_pair(ch_iter, od$from, od$to))
  t_cch_matrix <- system.time(cch_matrix <- get_distance_matrix(metric, matrix_from, matrix_to))
  t_dijkstra_matrix <- system.time(reference_matrix <- get_distance_matrix(
    graph_iter, matrix_from, matrix_to
  ))

  timings <- rbind(
    timings,
    data.table(iteration, phase = c(
      "cch_customize", "cch_pair", "dijkstra_pair", "ch_contract",
      "ch_pair", "cch_matrix", "dijkstra_matrix"
    ), elapsed_sec = c(
      t_customize[["elapsed"]], t_cch_pair[["elapsed"]], t_dijkstra[["elapsed"]],
      t_ch_contract[["elapsed"]], t_ch_pair[["elapsed"]],
      t_cch_matrix[["elapsed"]], t_dijkstra_matrix[["elapsed"]]
    ))
  )
  checks <- rbind(
    checks,
    distance_check_row(iteration, "cch_pair", reference_pair, cch_pair),
    distance_check_row(iteration, "ch_pair", reference_pair, ch_pair),
    distance_check_row(iteration, "cch_matrix", reference_matrix, cch_matrix)
  )
}

cat("timing means\n")
print(timings[, list(elapsed_sec = mean(elapsed_sec)), by = phase])
cat("correctness maxima\n")
print(checks[, list(mismatches = max(mismatches), max_abs_diff = max(max_abs_diff)), by = comparison])
