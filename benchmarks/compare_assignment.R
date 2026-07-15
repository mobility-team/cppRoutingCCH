source("benchmarks/bench_utils.R")

RcppParallel::setThreadOptions(numThreads = 1)

run_assignment <- function(graph, od, method, algorithm, max_gap, max_it, cch) {
  elapsed <- system.time({
    result <- try(assign_traffic(
      graph, od$from, od$to, od$demand,
      algorithm = algorithm,
      max_gap = max_gap,
      max_it = max_it,
      aon_method = cch_aon_method(method),
      cch = if (is_prepared_cch_method(method)) cch else NULL,
      verbose = FALSE
    ), silent = TRUE)
  })[["elapsed"]]

  list(elapsed = elapsed, result = result)
}

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[[1]] else "chicago"
pairs_n <- if (length(args) >= 2) as.integer(args[[2]]) else if (dataset == "roads") 300L else 100L
max_it <- if (length(args) >= 3) as.integer(args[[3]]) else 10L
max_gap <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.05
algorithm <- if (length(args) >= 5) args[[5]] else "cfw"
methods <- if (length(args) >= 6) parse_csv_arg(args[[6]]) else c("bi", "cch_prepared")
shapes <- if (length(args) >= 7) parse_csv_arg(args[[7]]) else c("sparse", "repeat_origins", "repeat_destinations")
capacity <- if (length(args) >= 8) as.numeric(args[[8]]) else 50
seed <- if (length(args) >= 9) as.integer(args[[9]]) else 42L

graph <- make_benchmark_graph(dataset, capacity = capacity, perturb = TRUE)
cch <- NULL
if (any(vapply(methods, is_prepared_cch_method, logical(1)))) {
  t_prepare <- system.time(cch <- cpp_cch_prepare(graph))
  cat("cch_prepare_elapsed_sec=", t_prepare[["elapsed"]], "\n")
}

cat("DATASET", dataset, "\n")
cat("algorithm=", algorithm, "max_it=", max_it, "max_gap=", max_gap,
    "pairs_n=", pairs_n, "capacity=", capacity, "\n")
cat("methods=", paste(methods, collapse = ","),
    "shapes=", paste(shapes, collapse = ","), "\n")

for (shape_index in seq_along(shapes)) {
  shape <- shapes[[shape_index]]
  od <- make_shaped_od(shape, graph$dict$ref, pairs_n, seed + shape_index)
  timings <- data.table()
  results <- list()

  for (method in methods) {
    gc()
    one <- run_assignment(graph, od, method, algorithm, max_gap, max_it, cch)
    if (inherits(one$result, "try-error")) {
      timings <- rbind(timings, data.table(
        method, elapsed_sec = one$elapsed, elapsed_per_iteration = NA_real_,
        gap = NA_real_, iteration = NA_integer_, error = as.character(one$result)
      ))
    } else {
      results[[method]] <- one$result
      timings <- rbind(timings, data.table(
        method, elapsed_sec = one$elapsed,
        elapsed_per_iteration = one$elapsed / one$result$iteration,
        gap = one$result$gap, iteration = one$result$iteration,
        error = NA_character_
      ))
    }
  }

  cat("SHAPE", shape, "rows=", nrow(od),
      "unique_from=", uniqueN(od$from), "unique_to=", uniqueN(od$to), "\n")
  print(timings)
  if ("bi" %in% names(results)) {
    checks <- rbindlist(lapply(names(results), function(method) {
      check <- flow_check(results[["bi"]]$data, results[[method]]$data)
      data.table(
        method,
        flow_mismatches = unname(check[["flow_mismatches"]]),
        max_flow_diff = unname(check[["max_flow_diff"]]),
        max_cost_diff = unname(check[["max_cost_diff"]])
      )
    }))
    print(checks)
  }
}
