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

build_vertex_flows <- function(od_flows_fp, graph_fp, transport_zones_fp, vertices,
                               od_vertex_map, congestion_flows_scaling_factor,
                               target_max_vehicles_per_od_endpoint,
                               retained_volume_share) {
  od_flows <- as.data.table(read_parquet(od_flows_fp))
  if (any(is.na(od_flows$vehicle_volume))) {
    stop("Cannot assign traffic, some OD flow volumes are NA.")
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

  list(
    flows = od_flows,
    zone_rows = nrow(read_parquet(od_flows_fp)),
    expanded_rows = expanded_rows,
    retained_rows = nrow(od_flows)
  )
}

run_method <- function(graph, vertex_flows, method, algorithm, max_gap, max_it, cch) {
  cch_arg <- if (method == "cch") cch else NULL
  args <- list(
    Graph = graph,
    from = vertex_flows$vertex_id_from,
    to = vertex_flows$vertex_id_to,
    demand = vertex_flows$vehicle_volume,
    algorithm = algorithm,
    aon_method = method,
    max_gap = max_gap,
    max_it = max_it,
    verbose = FALSE
  )
  if (!is.null(cch_arg)) args$cch <- cch_arg
  do.call(assign_traffic, args)
}

args <- commandArgs(trailingOnly = TRUE)

package_fp <- if (length(args) >= 1) args[[1]] else "D:/dev/mobility/mobility"
graph_fp <- if (length(args) >= 2) args[[2]] else "d:/data/mobility/projects/grand-geneve/path_graph_car/modified/afe031259465406e5d58d115b7674005-car-modified-path-graph"
transport_zones_fp <- if (length(args) >= 3) args[[3]] else "d:/data/mobility/projects/grand-geneve/b8f7a3c54bdeef48b24ee5ecd4dc1518-transport_zones.gpkg"
od_flows_fps <- if (length(args) >= 4) parse_csv_arg(args[[4]]) else c(
  "d:/data/mobility/projects/grand-geneve/od_flows/4e7bc12ff9ddad07bc037d705ca69491-vehicle_od_flows_road.parquet",
  "d:/data/mobility/projects/grand-geneve/od_flows/fa344627b54057422cc3679846021e96-vehicle_od_flows_road.parquet"
)
cch_fp <- if (length(args) >= 5) args[[5]] else "d:/data/mobility/projects/grand-geneve/path_graph_car/cch/07cf4e5aafb3d3f7eff6d5666da70f00-cch.rds"
methods <- if (length(args) >= 6) parse_csv_arg(args[[6]]) else c("cbi", "cch")
algorithm <- if (length(args) >= 7) args[[7]] else "cfw"
max_it <- if (length(args) >= 8) as.integer(args[[8]]) else 3L
max_gap <- if (length(args) >= 9) as.numeric(args[[9]]) else 0.15
congestion_flows_scaling_factor <- if (length(args) >= 10) as.numeric(args[[10]]) else 0.5
target_max_vehicles_per_od_endpoint <- if (length(args) >= 11) as.numeric(args[[11]]) else 3000.0
retained_volume_share <- if (length(args) >= 12) as.numeric(args[[12]]) else 0.9
num_threads <- if (length(args) >= 13) {
  as.integer(args[[13]])
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
cat("cch_fp=", cch_fp, "\n")
cat("methods=", paste(methods, collapse = ","), "\n")
cat("algorithm=", algorithm, "max_it=", max_it, "max_gap=", max_gap, "\n")
cat("num_threads=", num_threads, "\n")

graph_hash <- strsplit(basename(graph_fp), "-", fixed = TRUE)[[1]][1]
graph_loaded <- time_step("read_graph", read_cppr_graph(dirname(graph_fp), graph_hash))
graph <- graph_loaded$value

vertices <- read_parquet(file.path(dirname(dirname(graph_fp)), paste0(graph_hash, "-vertices.parquet")))
od_vertex_map <- as.data.table(read_parquet(file.path(dirname(dirname(graph_fp)), paste0(graph_hash, "-od-vertex-map.parquet"))))

cch <- NULL
if ("cch" %in% methods) {
  if (nzchar(cch_fp) && file.exists(cch_fp)) {
    cch_loaded <- time_step("read_cch", readRDS(cch_fp))
    cch <- cch_loaded$value
  } else {
    cch_prepared <- time_step("prepare_cch", cpp_cch_prepare(graph))
    cch <- cch_prepared$value
  }
}

all_results <- list()
timings <- data.frame()

for (od_flows_fp in od_flows_fps) {
  cat("\nOD_FILE", od_flows_fp, "\n")
  vertex_flow_step <- time_step(
    "build_vertex_flows",
    build_vertex_flows(
      od_flows_fp = od_flows_fp,
      graph_fp = graph_fp,
      transport_zones_fp = transport_zones_fp,
      vertices = vertices,
      od_vertex_map = od_vertex_map,
      congestion_flows_scaling_factor = congestion_flows_scaling_factor,
      target_max_vehicles_per_od_endpoint = target_max_vehicles_per_od_endpoint,
      retained_volume_share = retained_volume_share
    )
  )
  vertex_info <- vertex_flow_step$value
  vertex_flows <- vertex_info$flows

  cat(
    "od_shape zone_rows=", vertex_info$zone_rows,
    "expanded_rows=", vertex_info$expanded_rows,
    "retained_rows=", vertex_info$retained_rows,
    "unique_from=", length(unique(vertex_flows$vertex_id_from)),
    "unique_to=", length(unique(vertex_flows$vertex_id_to)),
    "total_demand=", sum(vertex_flows$vehicle_volume),
    "\n"
  )

  od_results <- list()
  for (method in methods) {
    method_step <- time_step(
      paste0("assign_", method),
      run_method(graph, vertex_flows, method, algorithm, max_gap, max_it, cch)
    )
    result <- method_step$value
    od_results[[method]] <- result
    timings <- rbind(
      timings,
      data.frame(
        od_file = basename(od_flows_fp),
        method = method,
        elapsed_sec = method_step$elapsed,
        assignment_iteration = result$iteration,
        gap = result$gap,
        total_flow = sum(result$data$flow),
        stringsAsFactors = FALSE
      )
    )
  }
  all_results[[basename(od_flows_fp)]] <- od_results
}

cat("\nTIMINGS\n")
print(timings)

if ("cbi" %in% methods) {
  cat("\nFLOW_DIFFS_VS_CBI\n")
  checks <- data.frame()
  for (od_file in names(all_results)) {
    ref <- all_results[[od_file]][["cbi"]]$data
    ref <- ref[order(ref$from, ref$to), ]
    for (method in names(all_results[[od_file]])) {
      candidate <- all_results[[od_file]][[method]]$data
      candidate <- candidate[order(candidate$from, candidate$to), ]
      checks <- rbind(
        checks,
        data.frame(
          od_file = od_file,
          method = method,
          max_flow_diff_vs_cbi = max(abs(ref$flow - candidate$flow)),
          max_cost_diff_vs_cbi = max(abs(ref$cost - candidate$cost)),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  print(checks)
}
