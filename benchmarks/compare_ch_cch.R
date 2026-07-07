source("benchmarks/bench_utils.R")

RcppParallel::setThreadOptions(numThreads = 1)

compare_one <- function(dataset, pairs_n = 1000, order_source = "default") {
  cat("DATASET", dataset, "\n")

  edges <- read_benchmark_edges(dataset)
  cat(
    "edges=", nrow(edges),
    "nodes=", length(unique(c(edges$from, edges$to))),
    "pairs=", pairs_n, "\n"
  )

  t_graph <- system.time(graph <- makegraph(edges, directed = TRUE))
  od <- make_random_od(graph$dict$ref, pairs_n)

  t_dijkstra <- system.time(dijkstra <- get_distance_pair(graph, od$from, od$to, algorithm = "Dijkstra"))

  t_ch_pre <- system.time(ch <- cpp_contract(graph, silent = TRUE))
  t_ch_query <- system.time(ch_dist <- get_distance_pair(ch, od$from, od$to))

  if (order_source == "ch") {
    stop("CH-derived order is not enabled: it caused excessive CCH topology growth in the current CCH preparation path")
  } else if (order_source == "degree") {
    stop("Static degree order is not enabled: it caused excessive CCH topology growth in the current CCH preparation path")
  } else if (order_source != "default") {
    stop("unknown order_source: ", order_source)
  }

  t_cch_prepare <- system.time(cch <- cpp_cch_prepare(graph))
  t_cch_customize <- system.time(metric <- cpp_cch_customize(cch))
  t_cch_query <- system.time(cch_dist <- get_distance_pair(metric, od$from, od$to))

  print(data.frame(
    phase = c(
      "makegraph",
      "ch_preprocess",
      "cch_prepare",
      "cch_customize",
      "dijkstra_query",
      "ch_query",
      "cch_query"
    ),
    elapsed_sec = c(
      t_graph[["elapsed"]],
      t_ch_pre[["elapsed"]],
      t_cch_prepare[["elapsed"]],
      t_cch_customize[["elapsed"]],
      t_dijkstra[["elapsed"]],
      t_ch_query[["elapsed"]],
      t_cch_query[["elapsed"]]
    )
  ))

  print(rbind(
    ch_vs_dijkstra = compare_dist(dijkstra, ch_dist),
    cch_vs_dijkstra = compare_dist(dijkstra, cch_dist)
  ))

  cat("sample distances\n")
  print(head(data.frame(from = od$from, to = od$to, dijkstra = dijkstra, ch = ch_dist, cch = cch_dist), 10))
}

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args)) args[[1]] else "chicago"
order_source <- if (length(args) >= 2) args[[2]] else "default"

if (dataset == "chicago") {
  compare_one("chicago", pairs_n = 1000, order_source = order_source)
} else if (dataset == "roads") {
  compare_one("roads", pairs_n = 100, order_source = order_source)
} else {
  stop("unknown dataset: ", dataset)
}
