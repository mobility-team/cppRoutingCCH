source("benchmarks/bench_utils.R")

RcppParallel::setThreadOptions(numThreads = 1)

args <- commandArgs(trailingOnly = TRUE)
dataset <- if (length(args) >= 1) args[[1]] else "chicago"
pairs_n <- if (length(args) >= 2) as.integer(args[[2]]) else if (dataset == "roads") 100L else 1000L
value_columns <- if (length(args) >= 3) as.integer(args[[3]]) else 3L
seed <- if (length(args) >= 4) as.integer(args[[4]]) else 42L

graph <- make_benchmark_graph(dataset, perturb = TRUE)
od <- make_random_od(graph$dict$ref, pairs_n, seed = seed)
edge_id <- seq_len(nrow(graph$data))
values <- as.data.frame(lapply(seq_len(value_columns), function(column) {
  (edge_id %% (97L + column)) + column / 10
}))
names(values) <- paste0("value", seq_len(value_columns))

t_cch_prepare <- system.time(cch <- cpp_cch_prepare(graph))
t_customize <- system.time(metric <- cpp_cch_customize(cch))
t_ch_prepare <- system.time(ch <- cpp_contract(graph, silent = TRUE))
t_cch <- system.time(cch_values <- get_path_values_pair(metric, od$from, od$to, values))
t_ch <- system.time(ch_values <- get_path_values_pair(ch, od$from, od$to, values))

t_repeated <- system.time({
  reference <- lapply(values, function(edge_values) {
    graph_value <- graph
    graph_value$attrib$aux <- edge_values
    get_distance_pair(
      graph_value, od$from, od$to,
      algorithm = "Dijkstra", aggregate_aux = TRUE
    )
  })
})
reference <- as.data.frame(reference)

checks <- rbind(
  cch_cost = compare_dist(
    get_distance_pair(graph, od$from, od$to, algorithm = "Dijkstra"),
    cch_values$cost
  ),
  ch_cost = compare_dist(cch_values$cost, ch_values$cost)
)
for (column in names(values)) {
  checks <- rbind(
    checks,
    setNames(compare_dist(reference[[column]], cch_values[[column]]), colnames(checks)),
    setNames(compare_dist(reference[[column]], ch_values[[column]]), colnames(checks))
  )
  rownames(checks)[(nrow(checks) - 1):nrow(checks)] <- paste(c("cch", "ch"), column, sep = "_")
}

cat("DATASET", dataset, "pairs_n=", pairs_n,
    "value_columns=", value_columns, "\n")
print(data.frame(
  phase = c("cch_prepare", "cch_customize", "ch_prepare", "repeated_dijkstra", "cch_path_values", "ch_path_values"),
  elapsed_sec = c(
    t_cch_prepare[["elapsed"]], t_customize[["elapsed"]], t_ch_prepare[["elapsed"]],
    t_repeated[["elapsed"]], t_cch[["elapsed"]], t_ch[["elapsed"]]
  )
))
print(checks)
