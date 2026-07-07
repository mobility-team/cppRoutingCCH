source("benchmarks/bench_utils.R")

RcppParallel::setThreadOptions(numThreads = 1)

is_prepared_cch_method <- function(method) {
  method == "cch_prepared"
}

cch_aon_method <- function(method) {
  switch(
    method,
    cch_prepared = "cch",
    method
  )
}

run_one <- function(graph, od, method, algorithm, max_gap, max_it, prepared_cch) {
  aon_method <- cch_aon_method(method)
  cch_arg <- if (is_prepared_cch_method(method)) prepared_cch else NULL

  elapsed <- system.time({
    result <- try(
      assign_traffic(
        graph,
        from = od$from,
        to = od$to,
        demand = od$demand,
        algorithm = algorithm,
        max_gap = max_gap,
        max_it = max_it,
        aon_method = aon_method,
        cch = cch_arg,
        verbose = FALSE
      ),
      silent = TRUE
    )
  })[["elapsed"]]

  if (inherits(result, "try-error")) {
    list(
      timing = data.table(
        method = method,
        elapsed_sec = elapsed,
        gap = NA_real_,
        iteration = NA_integer_,
        total_flow = NA_real_,
        error = as.character(result)
      ),
      result = NULL
    )
  } else {
    list(
      timing = data.table(
        method = method,
        elapsed_sec = elapsed,
        gap = result$gap,
        iteration = result$iteration,
        total_flow = sum(result$data$flow),
        error = NA_character_
      ),
      result = result
    )
  }
}

compare_shape <- function(graph, shape, od, methods, algorithm, max_gap, max_it, prepared_cch) {
  cat("SHAPE", shape, "\n")
  cat("od_rows=", nrow(od),
      "unique_from=", length(unique(od$from)),
      "unique_to=", length(unique(od$to)),
      "total_demand=", sum(od$demand), "\n")

  timings <- data.table()
  results <- list()

  for (method in methods) {
    gc()
    one <- run_one(graph, od, method, algorithm, max_gap, max_it, prepared_cch)
    timings <- rbind(timings, one$timing, fill = TRUE)
    if (!is.null(one$result)) results[[method]] <- one$result
  }

  timings[, elapsed_per_iteration := elapsed_sec / iteration]
  print(timings)

  if ("bi" %in% names(results)) {
    ref <- results[["bi"]]$data
    setorder(ref, from, to)
    checks <- data.table()
    for (method in names(results)) {
      candidate <- results[[method]]$data
      setorder(candidate, from, to)
      checks <- rbind(
        checks,
        data.table(
          method = method,
          max_flow_diff_vs_bi = max(abs(ref$flow - candidate$flow)),
          max_cost_diff_vs_bi = max(abs(ref$cost - candidate$cost)),
          flow_mismatches_vs_bi = sum(abs(ref$flow - candidate$flow) > 1e-7)
        )
      )
    }
    print(checks)
  }
}

compare_od_shapes <- function(dataset, pairs_n, max_it, max_gap,
                              algorithm, methods, shapes, capacity, perturb, seed) {
  cat("DATASET", dataset, "\n")
  cat("algorithm=", algorithm,
      "max_it=", max_it,
      "max_gap=", max_gap,
      "pairs_n=", pairs_n,
      "capacity=", capacity,
      "perturb=", perturb, "\n")
  cat("methods=", paste(methods, collapse = ","), "\n")
  cat("shapes=", paste(shapes, collapse = ","), "\n")

  graph <- make_benchmark_graph(dataset, capacity = capacity, perturb = perturb)

  prepared_cch <- NULL
  if (any(vapply(methods, is_prepared_cch_method, logical(1)))) {
    t_prepare <- system.time(prepared_cch <- cpp_cch_prepare(graph))
    cat("cch_prepare_elapsed_sec=", t_prepare[["elapsed"]], "\n")
  }

  nodes <- graph$dict$ref
  for (i in seq_along(shapes)) {
    od <- make_shaped_od(shapes[[i]], nodes, pairs_n, seed + i)
    compare_shape(graph, shapes[[i]], od, methods, algorithm, max_gap, max_it, prepared_cch)
  }
}

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[[1]] else "chicago"
pairs_n <- if (length(args) >= 2) as.integer(args[[2]]) else if (dataset == "roads") 300L else 100L
max_it <- if (length(args) >= 3) as.integer(args[[3]]) else 10L
max_gap <- if (length(args) >= 4) as.numeric(args[[4]]) else 0.05
algorithm <- if (length(args) >= 5) args[[5]] else "cfw"
methods <- if (length(args) >= 6) parse_csv_arg(args[[6]]) else c("bi", "d", "cbi", "cphast", "cch_prepared")
shapes <- if (length(args) >= 7) parse_csv_arg(args[[7]]) else c("sparse", "repeat_origins", "repeat_destinations", "volume_ranked")
capacity <- if (length(args) >= 8) as.numeric(args[[8]]) else 50
perturb <- if (length(args) >= 9) as.logical(args[[9]]) else TRUE
seed <- if (length(args) >= 10) as.integer(args[[10]]) else 42L

if (dataset == "chicago") {
  compare_od_shapes(
    "chicago", pairs_n, max_it, max_gap, algorithm, methods, shapes, capacity, perturb, seed
  )
} else if (dataset == "roads") {
  compare_od_shapes(
    "roads", pairs_n, max_it, max_gap, algorithm, methods, shapes, capacity, perturb, seed
  )
} else {
  stop("unknown dataset: ", dataset)
}
