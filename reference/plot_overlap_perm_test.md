# Plot overlap permutation test results as violin plots

Renders the core convert-weighted overlap figure: per trace type, a
violin of the permuted (null) pooled overlap fractions with an inner
boxplot, the observed pooled fraction as a red diamond, and the observed
counts (`n_overlap`/`n_converts`) labeled above each violin. Violins are
filled by overlap group (co-occurrence vs sequential).

Only the core figure is drawn. Cross-cluster/sequence-type aggregation,
significance annotations, and trace-type-specific dividers or headers
are left to the caller to add to the returned plot.

## Usage

``` r
plot_overlap_perm_test(
  perm_df,
  obs_df,
  title = "Overlap Permutation Test Results",
  subtitle = NULL,
  trace_labels = c(facility = "Facility", floor = "Floor", room = "Room", seq_facility =
    "Facility", seq_floor = "Floor", seq_room = "Room"),
  x_label = "Trace type",
  y_label = "Fraction of converts with overlap (convert-weighted)"
)
```

## Arguments

- perm_df:

  Long null-distribution data frame, one row per (trace type,
  permutation), with columns `trace_type` (factor; levels set x-axis
  order) and `overlap_fraction`.

- obs_df:

  Observed values, one row per trace type, with columns `trace_type`
  (same levels as `perm_df`), `overlap_fraction`, `n_overlap`,
  `n_converts`.

- title:

  Plot title.

- subtitle:

  Plot subtitle, or `NULL` for none.

- trace_labels:

  Named character vector mapping `trace_type` levels to x-axis labels.

- x_label, y_label:

  Axis titles.

## Value

A ggplot object with violin plots for each overlap type.
