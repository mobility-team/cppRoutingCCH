source("benchmarks/bench_utils.R")

RcppParallel::setThreadOptions(numThreads = 1)

compare_assignment <- function(dataset, pairs_n, max_it, algorithm, perturb = TRUE) {
  cat("DATASET", dataset, "\n")
  cat("algorithm=", algorithm, "max_it=", max_it, "pairs=", pairs_n, "perturb=", perturb, "\n")

  graph <- make_benchmark_graph(dataset, capacity = 50, perturb = perturb)
  od <- make_random_od(graph$dict$ref, pairs_n, demand = TRUE)

  t_bi <- system.time(bi <- assign_traffic(
    graph, od$from, od$to, od$demand,
    algorithm = algorithm, max_gap = 0, max_it = max_it, aon_method = "bi", verbose = FALSE
  ))

  t_cch <- system.time(cch <- assign_traffic(
    graph, od$from, od$to, od$demand,
    algorithm = algorithm, max_gap = 0, max_it = max_it, aon_method = "cch", verbose = FALSE
  ))

  bi_data <- sort_flow_table(bi$data)
  cch_data <- sort_flow_table(cch$data)

  print(data.frame(
    method = c("bi", "cch"),
    elapsed_sec = c(t_bi[["elapsed"]], t_cch[["elapsed"]]),
    gap = c(bi$gap, cch$gap),
    iteration = c(bi$iteration, cch$iteration),
    total_flow = c(sum(bi$data$flow), sum(cch$data$flow))
  ))

  print(c(
    max_flow_diff = max(abs(bi_data$flow - cch_data$flow)),
    max_cost_diff = max(abs(bi_data$cost - cch_data$cost)),
    flow_mismatches = sum(abs(bi_data$flow - cch_data$flow) > 1e-7)
  ))
}

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args)) args[[1]] else "chicago"
pairs_n <- if (length(args) >= 2) as.integer(args[[2]]) else if (dataset == "roads") 50L else 100L
max_it <- if (length(args) >= 3) as.integer(args[[3]]) else 3L
algorithm <- if (length(args) >= 4) args[[4]] else "msa"
perturb <- if (length(args) >= 5) as.logical(args[[5]]) else TRUE

if (dataset == "chicago") {
  compare_assignment("chicago", pairs_n, max_it, algorithm, perturb)
} else if (dataset == "roads") {
  compare_assignment("roads", pairs_n, max_it, algorithm, perturb)
} else {
  stop("unknown dataset: ", dataset)
}
