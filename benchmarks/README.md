# Bounded correctness and performance observations

These scripts separate prediction from summary I/O. They do not refit the
completed ecological case, reduce its trees/folds, or measure whole-workflow
acceleration. Results are observations on one machine, not general guarantees.

## Reproduce

Install castSDM >= 0.27.0 and the engines used in the saved fit. Set the project
library with `R_LIBS_USER` if needed. Run the benchmarks sequentially on an
otherwise idle machine; do not run package checks or other benchmarks alongside
them. `CAST_SOURCE_DIR` optionally uses `pkgload::load_all()` for development.

```sh
Rscript benchmarks/threshold_scan.R existing_comparison.tif benchmarks/output/io
Rscript benchmarks/paired_prediction.R completed_case_directory benchmarks/output/prediction
```

The first script also runs without a source argument using a synthetic raster.
Both scripts save all trial timings, verified numerical results and session
metadata. Temporary slice rasters stay in the chosen output directory. The
reference reader is retained only in the benchmark as a reproducible baseline;
it is not a production path. The prediction benchmark expects the versioned
single-species example's saved case, pool and scenario manifests.

## Observations: 7 November 2026

Windows 10 x64 build 19045, R 4.4.3, castSDM 0.27.0; Intel Xeon w7-3455
(24 physical / 48 logical cores). Prediction used four RF threads. Three warm
filesystem trials alternated method order, without concurrent R workloads.
Raw timing/count tables are in [results/2026-11-07](results/2026-11-07).

### Threshold-only I/O

A bounded 256 x 4861 cell strip (1,244,416 cells), beginning at source row 1600,
was copied from the saved late-SSP585 four-corner output. The source has 21
bands. The strip contains 234,497 valid four-corner cells and 130,988 supported
cells. Identical Float64 values, LZW compression and tiled geometry were used
for pixel- and band-interleaved copies. Three thresholds (0.3, 0.5, 0.7), all
eight contrasts and both full/supported domains had exactly equal counts.

| Method | Median elapsed, seconds | Range, seconds | Cumulative R allocation, MB |
|---|---:|---:|---:|
| 0.26-style short-strip `values()`, pixel layout | 7.45 | 7.43–7.56 | 335.81 |
| Persistent tile-aligned reader, same pixel layout | 1.04 | 0.83–1.05 | 309.48 |
| Persistent reader, band layout | 0.54 | 0.52–0.60 | 309.48 |

The same-file reader comparison is about 7.2 times faster here; reader plus
layout is about 13.8 times faster. The band file was 52.66 MB versus 50.78 MB
for pixel layout; faster selective reads did not imply a smaller artifact.
R allocation is a separate profiled run's cumulative allocation, excludes
GDAL/native allocation and is **not peak RSS**. Dataset creation, initial
input hashing and cold disk cache are not included in the timed scans.

Requested thresholds are now counted inside prediction blocks, so their later
summary uses saved counts rather than this scan at all. Tests verify that
additional diagnostic thresholds do not change raster hashes or call models.

### Paired predictions

The saved RF/MaxEnt/GAM fits and the same 4,000 eligible cells, weights,
Float64 calculations and clamp were used. The reference already shares current
predictions and uses exact additive GAM terms; it evaluates all three future
corners on all rows. The new path also reports sign agreement and support.

| Method | Median elapsed, seconds | Range, seconds |
|---|---:|---:|
| Reference paired block | 1.96 | 1.95–1.96 |
| Reuse unchanged corners, with added diagnostics | 1.93 | 1.92–1.94 |

Climate changed at all 4,000 cells; land use changed at 3,509. Required rows
per nonadditive engine fell from 12,000 to 11,018 (8.18%). Runtime was similar;
these three trials do not establish a material prediction speedup. Maximum
absolute discrepancies were 2.22e-16 for corner/contrast values and 3.85e-16
for paired SD. No model-performance improvement was tested.

## Evidence still needed for a manuscript benchmark

Use increasing grids, predictors, engines and scenario counts; include the
complete workload, setup/serialization/hash costs, peak process RSS with a
declared instrument, cold and resumed execution, and a second hardware/OS
setup. Compare relevant tools on identical fits, precision, domains and
outputs. The present strips and samples cannot support a claim that the
complete nine-scenario workflow is 7 or 14 times faster.
