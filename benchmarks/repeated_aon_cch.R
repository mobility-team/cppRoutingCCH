source("benchmarks/bench_utils.R")

RcppParallel::setThreadOptions(numThreads = 1)

compare_repeated_aon <- function(dataset, iterations, pairs_n, perturb = TRUE) {
  cat("DATASET", dataset, "\n")
  cat("iterations=", iterations, "pairs=", pairs_n, "perturb=", perturb, "\n")

  graph <- make_benchmark_graph(dataset, perturb = perturb)
  base_cost <- graph$data$dist
  od <- make_random_od(graph$dict$ref, pairs_n, demand = TRUE)

  t_cch_prepare <- system.time(cch <- cpp_cch_prepare(graph))

  timings <- data.frame(
    iteration = integer(),
    phase = character(),
    elapsed_sec = numeric()
  )

  checks <- data.frame(
    iteration = integer(),
    flow_mismatches = numeric(),
    max_flow_diff = numeric(),
    cch_loaded_diff = numeric(),
    bi_loaded_diff = numeric()
  )

  for (iteration in seq_len(iterations)) {
    weights <- make_assignment_like_weights(base_cost, iteration)
    graph_iter <- graph
    graph_iter$data$dist <- weights

    t_cch_customize <- system.time(metric <- cpp_cch_customize(cch, weights))
    t_cch_aon <- system.time(cch_aon <- get_aon(metric, od$from, od$to, od$demand))
    t_bi_aon <- system.time(bi_aon <- get_aon(graph_iter, od$from, od$to, od$demand, algorithm = "bi"))

    distances <- get_distance_pair(graph_iter, od$from, od$to, algorithm = "Dijkstra")
    sptt <- sum(od$demand * distances, na.rm = TRUE)
    cch_loaded <- sum(cch_aon$flow * cch_aon$cost)
    bi_loaded <- sum(bi_aon$flow * bi_aon$cost)

    cch_aon <- sort_flow_table(cch_aon)
    bi_aon <- sort_flow_table(bi_aon)
    flow_diff <- abs(cch_aon$flow - bi_aon$flow)

    timings <- rbind(
      timings,
      data.frame(iteration = iteration, phase = "cch_customize", elapsed_sec = t_cch_customize[["elapsed"]]),
      data.frame(iteration = iteration, phase = "cch_aon", elapsed_sec = t_cch_aon[["elapsed"]]),
      data.frame(iteration = iteration, phase = "bi_aon", elapsed_sec = t_bi_aon[["elapsed"]])
    )

    checks <- rbind(
      checks,
      data.frame(
        iteration = iteration,
        flow_mismatches = sum(flow_diff > 1e-7),
        max_flow_diff = max(flow_diff),
        cch_loaded_diff = cch_loaded - sptt,
        bi_loaded_diff = bi_loaded - sptt
      )
    )
  }

  cat("one-time prepare\n")
  print(data.frame(phase = "cch_prepare", elapsed_sec = t_cch_prepare[["elapsed"]]))

  cat("per-iteration totals\n")
  print(aggregate(elapsed_sec ~ phase, timings, sum))

  cat("per-iteration means\n")
  print(aggregate(elapsed_sec ~ phase, timings, mean))

  cat("correctness maxima\n")
  print(data.frame(
    flow_mismatches = max(checks$flow_mismatches),
    max_flow_diff = max(checks$max_flow_diff),
    cch_loaded_diff = max(abs(checks$cch_loaded_diff)),
    bi_loaded_diff = max(abs(checks$bi_loaded_diff))
  ))
}

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args)) args[[1]] else "chicago"
iterations <- if (length(args) >= 2) as.integer(args[[2]]) else 5L
pairs_n <- if (length(args) >= 3) as.integer(args[[3]]) else if (dataset == "roads") 50L else 1000L
perturb <- if (length(args) >= 4) as.logical(args[[4]]) else TRUE

if (dataset == "chicago") {
  compare_repeated_aon("chicago", iterations, pairs_n, perturb)
} else if (dataset == "roads") {
  compare_repeated_aon("roads", iterations, pairs_n, perturb)
} else {
  stop("unknown dataset: ", dataset)
}
