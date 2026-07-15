source("benchmarks/bench_utils.R")

RcppParallel::setThreadOptions(numThreads = 1)

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[[1]] else "chicago"
iterations <- if (length(args) >= 2) as.integer(args[[2]]) else 3L
pairs_n <- if (length(args) >= 3) as.integer(args[[3]]) else if (dataset == "roads") 300L else 100L
shapes <- if (length(args) >= 4) parse_csv_arg(args[[4]]) else c("sparse", "repeat_origins", "repeat_destinations")
seed <- if (length(args) >= 5) as.integer(args[[5]]) else 42L

graph <- make_benchmark_graph(dataset, perturb = TRUE)
base_cost <- graph$data$dist
t_prepare <- system.time(cch <- cpp_cch_prepare(graph))
timings <- data.table()
checks <- data.table()

cat("DATASET", dataset, "\n")
cat("iterations=", iterations, "pairs_n=", pairs_n,
    "shapes=", paste(shapes, collapse = ","), "\n")
cat("cch_prepare_elapsed_sec=", t_prepare[["elapsed"]], "\n")

for (iteration in seq_len(iterations)) {
  weights <- make_assignment_like_weights(base_cost, iteration)
  graph_iter <- graph
  graph_iter$data$dist <- weights
  t_customize <- system.time(metric <- cpp_cch_customize(cch, weights))
  timings <- rbind(timings, data.table(iteration, shape = "all", method = "cch_customize",
                                       elapsed_sec = t_customize[["elapsed"]]))

  for (shape_index in seq_along(shapes)) {
    shape <- shapes[[shape_index]]
    od <- make_shaped_od(shape, graph$dict$ref, pairs_n, seed + shape_index)
    t_d <- system.time(reference <- get_aon(
      graph_iter, od$from, od$to, od$demand, algorithm = "d"
    ))
    t_cch <- system.time(candidate <- get_aon(metric, od$from, od$to, od$demand))

    timings <- rbind(
      timings,
      data.table(iteration, shape, method = c("d", "cch"),
                 elapsed_sec = c(t_d[["elapsed"]], t_cch[["elapsed"]]))
    )
    checks <- rbind(
      checks,
      {
        check <- flow_check(reference, candidate)
        data.table(
          iteration,
          shape,
          flow_mismatches = unname(check[["flow_mismatches"]]),
          max_flow_diff = unname(check[["max_flow_diff"]]),
          max_cost_diff = unname(check[["max_cost_diff"]])
        )
      }
    )
  }
}

cat("timing means\n")
print(timings[, list(elapsed_sec = mean(elapsed_sec)), by = list(shape, method)])
cat("correctness maxima\n")
print(checks[, list(
  flow_mismatches = max(flow_mismatches),
  max_flow_diff = max(max_flow_diff),
  max_cost_diff = max(max_cost_diff)
), by = shape])
