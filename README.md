# cppRoutingCCH

`cppRoutingCCH` is a Mobility-oriented fork of
[`vlarmet/cppRouting`](https://github.com/vlarmet/cppRouting).

The R package name remains `cppRouting`, so it can be used as a drop-in
dependency where Mobility needs the original cppRouting API plus the CCH
extensions described below.

For the full upstream package presentation, examples, and original routing
documentation, see the original cppRouting README:

https://github.com/vlarmet/cppRouting#readme

## What This Fork Adds

This fork adds a customizable contraction hierarchy (CCH) backend for
transport models where the graph topology is stable but edge costs change
repeatedly, for example during congestion assignment.

The main additions are:

- `cpp_cch_prepare()` prepares reusable CCH topology from a `makegraph()` graph.
- `cpp_cch_customize()` applies current edge weights to a prepared CCH.
- `get_distance_pair()`, `get_distance_matrix()`, and `get_aon()` can query a
  customized CCH metric.
- `get_path_values_pair()` routes on one cost and accumulates one or more edge
  value columns, such as distance, along the selected paths.
- `assign_traffic(..., aon_method = "cch", cch = cch)` uses the prepared CCH
  topology in link-based traffic assignment, so repeated cost updates avoid
  rebuilding a full contraction hierarchy.

## Installation

From GitHub:

```r
install.packages("remotes")
remotes::install_github("mobility-team/cppRoutingCCH")
```

After publication on R-universe, Mobility should install this fork from the
Mobility R-universe repository before falling back to CRAN.

## Basic CCH Usage

```r
library(cppRouting)

edges <- data.frame(
  from = c("a", "b", "a"),
  to = c("b", "c", "c"),
  time = c(1, 2, 5),
  distance = c(10, 20, 30)
)

graph <- makegraph(edges[, c("from", "to", "time")], directed = TRUE)

# Prepare once for this graph topology.
cch <- cpp_cch_prepare(graph)

# Customize whenever edge costs change.
metric <- cpp_cch_customize(cch, weights = graph$data$dist)

get_distance_pair(metric, from = "a", to = "c")

get_path_values_pair(
  metric,
  from = "a",
  to = "c",
  values = data.frame(distance = edges$distance)
)
```

## Traffic Assignment

For repeated congestion iterations, prepare the CCH once and pass it to
`assign_traffic()`:

```r
trips <- data.frame(
  from = c("a", "a"),
  to = c("c", "b"),
  demand = c(100, 50)
)

traffic <- assign_traffic(
  Graph = graph,
  from = trips$from,
  to = trips$to,
  demand = trips$demand,
  algorithm = "cfw",
  aon_method = "cch",
  cch = cch
)
```

## Versioning

This fork currently publishes `cppRouting` version `3.2.1`, based on upstream
cppRouting `3.2` plus the CCH additions.

## License

This fork follows the upstream cppRouting license: GPL (>= 2).
