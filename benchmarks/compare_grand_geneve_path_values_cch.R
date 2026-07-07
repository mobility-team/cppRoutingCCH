source("benchmarks/bench_utils.R")
library(arrow)

time_step <- function(label, expr) {
  gc()
  elapsed <- system.time({
    value <- force(expr)
  })[["elapsed"]]
  cat(label, "elapsed_sec=", elapsed, "\n")
  list(value = value, elapsed = elapsed)
}

build_vertex_pairs <- function(od_flows_fp, transport_zones_fp, vertices,
                               od_vertex_map, congestion_flows_scaling_factor,
                               target_max_vehicles_per_od_endpoint,
                               retained_volume_share, max_pairs) {
  od_flows <- as.data.table(read_parquet(od_flows_fp))
  if (any(is.na(od_flows$vehicle_volume))) {
    stop("Cannot build vertex pairs, some OD flow volumes are NA.")
  }

  buildings_fp <- file.path(
    dirname(transport_zones_fp),
    paste0(
      gsub("-transport_zones.gpkg", "", basename(transport_zones_fp)),
      "-transport_zones_buildings.parquet"
    )
  )

  buildings <- as.data.table(read_parquet(buildings_fp))
  buildings[, building_id := 1:.N]
  buildings <- merge(
    buildings,
    od_vertex_map[, list(building_id, vertex_id)],
    by = "building_id",
    all.x = TRUE,
    sort = FALSE
  )
  if (any(is.na(buildings$vertex_id))) {
    stop("OD vertex map is incomplete for the current buildings sample.")
  }
  buildings <- merge(
    buildings[, list(building_id, transport_zone_id, n_clusters, weight, vertex_id)],
    vertices,
    by = "vertex_id"
  )

  vertex_pairs <- tz_pairs_to_vertex_pairs(
    tz_id_from = od_flows$from,
    tz_id_to = od_flows$to,
    buildings = buildings,
    vehicle_volume = od_flows$vehicle_volume,
    congestion_flows_scaling_factor = congestion_flows_scaling_factor,
    target_max_vehicles_per_od_endpoint = target_max_vehicles_per_od_endpoint
  )

  od_flows <- merge(
    od_flows,
    vertex_pairs,
    by.x = c("from", "to"),
    by.y = c("tz_id_from", "tz_id_to")
  )
  od_flows[, vehicle_volume := vehicle_volume * weight]
  expanded_rows <- nrow(od_flows)

  od_flows <- od_flows[order(-vehicle_volume)]
  od_flows[, cum_share := cumsum(vehicle_volume) / sum(vehicle_volume)]
  od_flows[, flow_rank := .I]
  od_flows <- od_flows[cum_share <= retained_volume_share | flow_rank == 1]
  od_flows[, vehicle_volume := vehicle_volume / retained_volume_share * congestion_flows_scaling_factor]
  od_flows[, flow_rank := NULL]

  if (!is.na(max_pairs) && max_pairs > 0 && nrow(od_flows) > max_pairs) {
    od_flows <- od_flows[1:max_pairs]
  }

  list(
    pairs = unique(od_flows[, list(from = vertex_id_from, to = vertex_id_to)]),
    zone_rows = nrow(read_parquet(od_flows_fp)),
    expanded_rows = expanded_rows,
    retained_rows = nrow(od_flows)
  )
}

build_value_table <- function(graph) {
  values <- data.frame(time = graph$data$dist)

  attrib <- graph$attrib
  for (name in names(attrib)) {
    value <- attrib[[name]]
    if (is.numeric(value) && length(value) == nrow(graph$data)) {
      values[[name]] <- value
    }
  }

  if (!"aux" %in% names(values)) {
    stop("Graph has no edge-length aux column to benchmark.")
  }

  values
}

run_repeated_aux <- function(graph, from, to, values) {
  result <- data.frame(
    from = as.character(from),
    to = as.character(to),
    cost = get_distance_pair(graph, from, to, aggregate_aux = FALSE),
    stringsAsFactors = FALSE
  )

  for (name in names(values)) {
    graph$attrib$aux <- values[[name]]
    result[[name]] <- get_distance_pair(graph, from, to, aggregate_aux = TRUE)
  }

  result
}

run_cch_path_values <- function(graph, cch, from, to, values) {
  metric <- cpp_cch_customize(cch, weights = graph$data$dist)
  get_path_values_pair(metric, from, to, values)
}

summarise_differences <- function(reference, candidate, value_names) {
  checks <- data.frame(
    column = c("cost", value_names),
    max_abs_diff = NA_real_,
    mismatches_gt_1e_7 = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(checks))) {
    column <- checks$column[[i]]
    diff <- abs(reference[[column]] - candidate[[column]])
    checks$max_abs_diff[[i]] <- max(diff, na.rm = TRUE)
    checks$mismatches_gt_1e_7[[i]] <- sum(diff > 1e-7, na.rm = TRUE)
  }

  checks
}

args <- commandArgs(trailingOnly = TRUE)

package_fp <- if (length(args) >= 1) args[[1]] else "D:/dev/mobility/mobility"
graph_fp <- if (length(args) >= 2) args[[2]] else "d:/data/mobility/projects/grand-geneve/path_graph_car/modified/afe031259465406e5d58d115b7674005-car-modified-path-graph"
transport_zones_fp <- if (length(args) >= 3) args[[3]] else "d:/data/mobility/projects/grand-geneve/b8f7a3c54bdeef48b24ee5ecd4dc1518-transport_zones.gpkg"
od_flows_fp <- if (length(args) >= 4) args[[4]] else "d:/data/mobility/projects/grand-geneve/od_flows/4e7bc12ff9ddad07bc037d705ca69491-vehicle_od_flows_road.parquet"
cch_fp <- if (length(args) >= 5) args[[5]] else "d:/data/mobility/projects/grand-geneve/path_graph_car/cch/07cf4e5aafb3d3f7eff6d5666da70f00-cch.rds"
max_pairs <- if (length(args) >= 6) as.integer(args[[6]]) else 5000L
value_names <- if (length(args) >= 7) parse_csv_arg(args[[7]]) else c("time", "aux")
congestion_flows_scaling_factor <- if (length(args) >= 8) as.numeric(args[[8]]) else 0.5
target_max_vehicles_per_od_endpoint <- if (length(args) >= 9) as.numeric(args[[9]]) else 3000.0
retained_volume_share <- if (length(args) >= 10) as.numeric(args[[10]]) else 0.9
num_threads <- if (length(args) >= 11) {
  as.integer(args[[11]])
} else {
  as.integer(Sys.getenv("CPPROUTING_BENCH_THREADS", RcppParallel::defaultNumThreads()))
}
if (is.na(num_threads) || num_threads < 1L) {
  stop("num_threads must be a positive integer.")
}
RcppParallel::setThreadOptions(numThreads = num_threads)

source(file.path(package_fp, "transport", "graphs", "core", "cpprouting_io.R"))
source(file.path(package_fp, "transport", "graphs", "congested", "tz_pairs_to_vertex_pairs.R"))

cat("package_fp=", package_fp, "\n")
cat("graph_fp=", graph_fp, "\n")
cat("transport_zones_fp=", transport_zones_fp, "\n")
cat("od_flows_fp=", od_flows_fp, "\n")
cat("cch_fp=", cch_fp, "\n")
cat("max_pairs=", max_pairs, "\n")
cat("value_names=", paste(value_names, collapse = ","), "\n")
cat("num_threads=", num_threads, "\n")

graph_hash <- strsplit(basename(graph_fp), "-", fixed = TRUE)[[1]][1]
graph_loaded <- time_step("read_graph", read_cppr_graph(dirname(graph_fp), graph_hash))
graph <- graph_loaded$value

vertices <- read_parquet(file.path(dirname(dirname(graph_fp)), paste0(graph_hash, "-vertices.parquet")))
od_vertex_map <- as.data.table(read_parquet(file.path(dirname(dirname(graph_fp)), paste0(graph_hash, "-od-vertex-map.parquet"))))

pairs_step <- time_step(
  "build_vertex_pairs",
  build_vertex_pairs(
    od_flows_fp = od_flows_fp,
    transport_zones_fp = transport_zones_fp,
    vertices = vertices,
    od_vertex_map = od_vertex_map,
    congestion_flows_scaling_factor = congestion_flows_scaling_factor,
    target_max_vehicles_per_od_endpoint = target_max_vehicles_per_od_endpoint,
    retained_volume_share = retained_volume_share,
    max_pairs = max_pairs
  )
)
pairs_info <- pairs_step$value
pairs <- pairs_info$pairs

cat(
  "od_shape zone_rows=", pairs_info$zone_rows,
  "expanded_rows=", pairs_info$expanded_rows,
  "retained_rows=", pairs_info$retained_rows,
  "bench_pairs=", nrow(pairs),
  "unique_from=", length(unique(pairs$from)),
  "unique_to=", length(unique(pairs$to)),
  "\n"
)

values <- build_value_table(graph)
missing_values <- setdiff(value_names, names(values))
if (length(missing_values) > 0) {
  stop("Missing graph value columns: ", paste(missing_values, collapse = ", "))
}
values <- values[value_names]

cch <- NULL
if (nzchar(cch_fp) && file.exists(cch_fp)) {
  cch_loaded <- time_step("read_cch", readRDS(cch_fp))
  cch <- cch_loaded$value
} else {
  cch_prepared <- time_step("prepare_cch", cpp_cch_prepare(graph))
  cch <- cch_prepared$value
}

repeated_step <- time_step(
  "repeated_get_distance_pair",
  run_repeated_aux(graph, pairs$from, pairs$to, values)
)
repeated <- repeated_step$value

cch_step <- time_step(
  "cch_get_path_values_pair",
  run_cch_path_values(graph, cch, pairs$from, pairs$to, values)
)
cch_values <- cch_step$value

timings <- data.frame(
  method = c("repeated_get_distance_pair", "cch_get_path_values_pair"),
  elapsed_sec = c(repeated_step$elapsed, cch_step$elapsed),
  rows = nrow(pairs),
  value_columns = ncol(values),
  stringsAsFactors = FALSE
)

cat("\nTIMINGS\n")
print(timings)

cat("\nDIFFS_VS_REPEATED\n")
print(summarise_differences(repeated, cch_values, names(values)))
