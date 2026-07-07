source("benchmarks/bench_utils.R")

RcppParallel::setThreadOptions(numThreads = 1)

compare_repeated <- function(dataset, iterations, pairs_n) {
  cat("DATASET", dataset, "\n")
  cat("iterations=", iterations, "pairs=", pairs_n, "\n")

  graph <- make_benchmark_graph(dataset)
  base_cost <- graph$data$dist
  od <- make_random_od(graph$dict$ref, pairs_n)

  t_cch_prepare <- system.time(cch <- cpp_cch_prepare(graph))

  timings <- data.frame(
    iteration = integer(),
    phase = character(),
    elapsed_sec = numeric()
  )

  checks <- data.frame(
    iteration = integer(),
    comparison = character(),
    mismatches = numeric(),
    max_abs_diff = numeric()
  )

  for (iteration in seq_len(iterations)) {
    weights <- make_assignment_like_weights(base_cost, iteration)

    cch_customize_time <- system.time(cch_metric <- cpp_cch_customize(cch, weights))
    cch_query_time <- system.time(cch_dist <- get_distance_pair(cch_metric, od$from, od$to))

    graph_iter <- graph
    graph_iter$data$dist <- weights
    ch_contract_time <- system.time(ch <- cpp_contract(graph_iter, silent = TRUE))
    ch_query_time <- system.time(ch_dist <- get_distance_pair(ch, od$from, od$to))

    dijkstra_time <- system.time(dijkstra <- get_distance_pair(graph_iter, od$from, od$to, algorithm = "Dijkstra"))

    timings <- rbind(
      timings,
      data.frame(iteration = iteration, phase = "cch_customize", elapsed_sec = cch_customize_time[["elapsed"]]),
      data.frame(iteration = iteration, phase = "cch_query", elapsed_sec = cch_query_time[["elapsed"]]),
      data.frame(iteration = iteration, phase = "ch_contract", elapsed_sec = ch_contract_time[["elapsed"]]),
      data.frame(iteration = iteration, phase = "ch_query", elapsed_sec = ch_query_time[["elapsed"]]),
      data.frame(iteration = iteration, phase = "dijkstra_query", elapsed_sec = dijkstra_time[["elapsed"]])
    )

    cch_check <- compare_dist(dijkstra, cch_dist)
    ch_check <- compare_dist(dijkstra, ch_dist)
    checks <- rbind(
      checks,
      data.frame(iteration = iteration, comparison = "cch_vs_dijkstra",
                 mismatches = cch_check[["mismatches"]], max_abs_diff = cch_check[["max_abs_diff"]]),
      data.frame(iteration = iteration, comparison = "ch_vs_dijkstra",
                 mismatches = ch_check[["mismatches"]], max_abs_diff = ch_check[["max_abs_diff"]])
    )
  }

  cat("one-time prepare\n")
  print(data.frame(phase = "cch_prepare", elapsed_sec = t_cch_prepare[["elapsed"]]))

  cat("per-iteration totals\n")
  print(aggregate(elapsed_sec ~ phase, timings, sum))

  cat("per-iteration means\n")
  print(aggregate(elapsed_sec ~ phase, timings, mean))

  cat("correctness\n")
  print(aggregate(cbind(mismatches, max_abs_diff) ~ comparison, checks, max))
}

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args)) args[[1]] else "chicago"
iterations <- if (length(args) >= 2) as.integer(args[[2]]) else 5L
pairs_n <- if (length(args) >= 3) as.integer(args[[3]]) else if (dataset == "roads") 100L else 1000L

if (dataset == "chicago") {
  compare_repeated("chicago", iterations, pairs_n)
} else if (dataset == "roads") {
  compare_repeated("roads", iterations, pairs_n)
} else {
  stop("unknown dataset: ", dataset)
}
