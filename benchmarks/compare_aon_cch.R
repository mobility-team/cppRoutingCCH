source("benchmarks/bench_utils.R")

RcppParallel::setThreadOptions(numThreads = 1)

compare_aon <- function(dataset, pairs_n, perturb = FALSE) {
  cat("DATASET", dataset, "\n")
  cat("pairs=", pairs_n, "perturb=", perturb, "\n")

  graph <- make_benchmark_graph(dataset, perturb = perturb)
  weights <- graph$data$dist
  od <- make_random_od(graph$dict$ref, pairs_n, demand = TRUE)

  t_prepare <- system.time(cch <- cpp_cch_prepare(graph))
  t_customize <- system.time(metric <- cpp_cch_customize(cch, weights))
  t_cch_aon <- system.time(cch_aon <- get_aon(metric, od$from, od$to, od$demand))
  t_bi_aon <- system.time(bi_aon <- get_aon(graph, od$from, od$to, od$demand, algorithm = "bi"))

  distances <- get_distance_pair(graph, od$from, od$to, algorithm = "Dijkstra")
  sptt <- sum(od$demand * distances, na.rm = TRUE)
  cch_loaded <- sum(cch_aon$flow * cch_aon$cost)
  bi_loaded <- sum(bi_aon$flow * bi_aon$cost)

  cch_aon <- sort_flow_table(cch_aon)
  bi_aon <- sort_flow_table(bi_aon)
  flow_diff <- abs(cch_aon$flow - bi_aon$flow)

  print(data.frame(
    phase = c("cch_prepare", "cch_customize", "cch_aon", "bi_aon"),
    elapsed_sec = c(t_prepare[["elapsed"]], t_customize[["elapsed"]], t_cch_aon[["elapsed"]], t_bi_aon[["elapsed"]])
  ))

  print(c(
    sptt = sptt,
    cch_loaded = cch_loaded,
    cch_loaded_diff = cch_loaded - sptt,
    bi_loaded = bi_loaded,
    bi_loaded_diff = bi_loaded - sptt,
    flow_mismatches = sum(flow_diff > 1e-7),
    max_flow_diff = max(flow_diff)
  ))
}

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args)) args[[1]] else "chicago"
pairs_n <- if (length(args) >= 2) as.integer(args[[2]]) else if (dataset == "roads") 100L else 1000L
perturb <- if (length(args) >= 3) as.logical(args[[3]]) else FALSE

if (dataset == "chicago") {
  compare_aon("chicago", pairs_n, perturb)
} else if (dataset == "roads") {
  compare_aon("roads", pairs_n, perturb)
} else {
  stop("unknown dataset: ", dataset)
}
