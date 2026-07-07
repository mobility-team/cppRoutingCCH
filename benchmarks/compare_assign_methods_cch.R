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

compare_assignment_methods <- function(dataset, pairs_n, max_it, algorithm,
                                       methods, capacity, perturb = TRUE) {
  cat("DATASET", dataset, "\n")
  cat("algorithm=", algorithm, "max_it=", max_it, "pairs=", pairs_n,
      "capacity=", capacity, "perturb=", perturb, "\n")
  cat("methods=", paste(methods, collapse = ","), "\n")

  graph <- make_benchmark_graph(dataset, capacity = capacity, perturb = perturb)
  od <- make_random_od(graph$dict$ref, pairs_n, demand = TRUE)

  results <- list()
  timings <- data.frame()
  prepared_cch <- NULL
  if (any(vapply(methods, is_prepared_cch_method, logical(1)))) {
    t_prepare <- system.time(prepared_cch <- cpp_cch_prepare(graph))
    cat("cch_prepare_elapsed_sec=", t_prepare[["elapsed"]], "\n")
  }

  for (method in methods) {
    gc()
    aon_method <- cch_aon_method(method)
    cch_arg <- if (is_prepared_cch_method(method)) prepared_cch else NULL
    elapsed <- system.time({
      result <- try(
        assign_traffic(
          graph, od$from, od$to, od$demand,
          algorithm = algorithm,
          max_gap = 0,
          max_it = max_it,
          aon_method = aon_method,
          cch = cch_arg,
          verbose = FALSE
        ),
        silent = TRUE
      )
    })[["elapsed"]]

    if (inherits(result, "try-error")) {
      timings <- rbind(
        timings,
        data.frame(
          method = method,
          elapsed_sec = elapsed,
          elapsed_per_reported_iteration = NA_real_,
          gap = NA_real_,
          iteration = NA_integer_,
          total_flow = NA_real_,
          error = as.character(result)
        )
      )
    } else {
      results[[method]] <- result
      timings <- rbind(
        timings,
        data.frame(
          method = method,
          elapsed_sec = elapsed,
          elapsed_per_reported_iteration = elapsed / result$iteration,
          gap = result$gap,
          iteration = result$iteration,
          total_flow = sum(result$data$flow),
          error = NA_character_
        )
      )
    }
  }

  print(timings)

  if ("bi" %in% names(results)) {
    ref <- sort_flow_table(results[["bi"]]$data)
    checks <- data.frame()
    for (method in names(results)) {
      candidate <- sort_flow_table(results[[method]]$data)
      checks <- rbind(
        checks,
        data.frame(
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

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args)) args[[1]] else "chicago"
pairs_n <- if (length(args) >= 2) as.integer(args[[2]]) else if (dataset == "roads") 50L else 100L
max_it <- if (length(args) >= 3) as.integer(args[[3]]) else 3L
algorithm <- if (length(args) >= 4) args[[4]] else "msa"
methods <- if (length(args) >= 5) parse_csv_arg(args[[5]]) else c("bi", "cbi", "cphast", "cch")
capacity <- if (length(args) >= 6) as.numeric(args[[6]]) else 50
perturb <- if (length(args) >= 7) as.logical(args[[7]]) else TRUE

if (dataset == "chicago") {
  compare_assignment_methods(
    "chicago", pairs_n, max_it, algorithm, methods, capacity, perturb
  )
} else if (dataset == "roads") {
  compare_assignment_methods(
    "roads", pairs_n, max_it, algorithm, methods, capacity, perturb
  )
} else {
  stop("unknown dataset: ", dataset)
}
