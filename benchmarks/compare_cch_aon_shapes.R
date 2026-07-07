source("benchmarks/bench_utils.R")

RcppParallel::setThreadOptions(numThreads = 1)

compare_shape <- function(graph, metric, shape, od) {
  cat("SHAPE", shape, "\n")
  cat("od_rows=", nrow(od),
      "unique_from=", length(unique(od$from)),
      "unique_to=", length(unique(od$to)),
      "total_demand=", sum(od$demand), "\n")

  t_d <- system.time(aon_d <- get_aon(graph, od$from, od$to, od$demand, algorithm = "d"))
  t_cch <- system.time(aon_cch <- get_aon(metric, od$from, od$to, od$demand))

  aon_d <- sort_flow_table(aon_d)
  aon_cch <- sort_flow_table(aon_cch)

  print(data.table(
    method = c("d", "cch"),
    elapsed_sec = c(t_d[["elapsed"]], t_cch[["elapsed"]]),
    loaded_cost = c(
      sum(aon_d$flow * aon_d$cost),
      sum(aon_cch$flow * aon_cch$cost)
    )
  ))

  print(data.table(
    method = "cch",
    max_flow_diff_vs_d = max(abs(aon_cch$flow - aon_d$flow)),
    flow_mismatches_vs_d = sum(abs(aon_cch$flow - aon_d$flow) > 1e-7)
  ))
}

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[[1]] else "chicago"
pairs_n <- if (length(args) >= 2) as.integer(args[[2]]) else if (dataset == "roads") 300L else 100L
perturb <- if (length(args) >= 4) as.logical(args[[4]]) else TRUE
seed <- if (length(args) >= 5) as.integer(args[[5]]) else 42L

shapes <- if (length(args) >= 3) parse_csv_arg(args[[3]]) else c("sparse", "repeat_origins", "repeat_destinations", "volume_ranked")
graph <- make_benchmark_graph(dataset, perturb = perturb)
t_prepare <- system.time(cch <- cpp_cch_prepare(graph))
t_customize <- system.time(metric <- cpp_cch_customize(cch))
cat("DATASET", dataset, "\n")
cat("pairs_n=", pairs_n, "perturb=", perturb, "\n")
cat("cch_prepare_elapsed_sec=", t_prepare[["elapsed"]], "\n")
cat("cch_customize_elapsed_sec=", t_customize[["elapsed"]], "\n")

nodes <- graph$dict$ref
for (i in seq_along(shapes)) {
  od <- make_shaped_od(shapes[[i]], nodes, pairs_n, seed + i)
  compare_shape(graph, metric, shapes[[i]], od)
}
