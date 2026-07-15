# Benchmarks

These scripts use the public `chicago` and `roads` datasets and an installed
`cppRoutingCCH` package. Run them from the repository root.

## Preparation and queries

Compares CH and CCH preparation, repeated cost updates, pair queries, and CCH
distance matrices against Dijkstra. It also reports CH edge and CCH arc counts,
which are the first indicators to watch when evaluating a new CCH ordering.

```powershell
mamba run -n cpprouting-cch Rscript benchmarks/compare_prepare_query.R chicago 3 1000 25
mamba run -n cpprouting-cch Rscript benchmarks/compare_prepare_query.R roads 3 100 25
```

The repository includes the public-domain DIMACS Florida travel-time graph as
a 14 MB gzip file. Benchmark scripts decompress it to a temporary `.gr` file
and remove that file immediately after loading the graph:

```powershell
mamba run -n cpprouting-cch Rscript benchmarks/compare_prepare_query.R benchmarks/data/USA-road-t.FLA.gr.gz 1 100 25
```

The same dataset path works with every consolidated benchmark:

```powershell
mamba run -n cpprouting-cch Rscript benchmarks/compare_aon.R benchmarks/data/USA-road-t.FLA.gr.gz 1 300
mamba run -n cpprouting-cch Rscript benchmarks/compare_assignment.R benchmarks/data/USA-road-t.FLA.gr.gz 300 3 0.05 cfw 'bi,cch_prepared' 'sparse,repeat_origins,repeat_destinations'
mamba run -n cpprouting-cch Rscript benchmarks/compare_path_values.R benchmarks/data/USA-road-t.FLA.gr.gz 100 3
```

On the reference single-thread run, CCH preparation took 5.1 seconds versus
24.6 seconds for CH. Prepared-CCH congested assignment was 11.1x to 13.4x
faster than bidirectional Dijkstra across the three OD shapes, with identical
flows. These are baseline measurements, not fixed performance guarantees.

## All-or-Nothing assignment

Compares the public CCH elimination-tree query with Dijkstra across sparse,
repeated-origin, and repeated-destination OD shapes, including repeated cost
customization.

```powershell
mamba run -n cpprouting-cch Rscript benchmarks/compare_aon.R chicago 3 100
mamba run -n cpprouting-cch Rscript benchmarks/compare_aon.R roads 3 300
```

## Congested assignment

Compares prepared CCH against the selected assignment methods and OD shapes.
Quote comma-separated method and shape lists in PowerShell.

```powershell
mamba run -n cpprouting-cch Rscript benchmarks/compare_assignment.R chicago 100 3 0.05 cfw 'bi,cch_prepared' 'sparse,repeat_origins,repeat_destinations'
mamba run -n cpprouting-cch Rscript benchmarks/compare_assignment.R roads 300 3 0.05 cfw 'bi,cch_prepared' 'sparse,repeat_origins'
```

## Path values

Compares repeated auxiliary-value extraction with one multi-column CH or CCH
path-values query.

```powershell
mamba run -n cpprouting-cch Rscript benchmarks/compare_path_values.R chicago 1000 3
mamba run -n cpprouting-cch Rscript benchmarks/compare_path_values.R roads 100 3
```
