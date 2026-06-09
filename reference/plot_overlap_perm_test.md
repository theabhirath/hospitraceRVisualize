# Plot overlap permutation test results as violin plots

Renders the core convert-weighted overlap figure: for each trace type a
violin of the permuted (null) pooled overlap fractions with a boxplot
inside, the observed pooled fraction as a red diamond, and the observed
counts (`n_overlap`/`n_converts`) labelled above each violin. Violins
are filled by overlap group (co-occurrence vs sequential).

This function draws only the core figure. Aggregation across clusters /
sequence types, significance/p-value annotations, and any group dividers
or headers that depend on the specific trace-type layout are left to the
caller to compute and add to the returned plot.

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

  A long data frame of the permuted (null) distribution with columns
  `trace_type` (a factor whose levels define the x-axis order) and
  `overlap_fraction`. One row per (trace type, permutation).

- obs_df:

  A data frame of observed values with one row per trace type: columns
  `trace_type` (a factor with the same levels as `perm_df`),
  `overlap_fraction`, `n_overlap` and `n_converts`.

- title:

  Plot title.

- subtitle:

  Plot subtitle. Defaults to `NULL` (no subtitle).

- trace_labels:

  A named character vector mapping `trace_type` levels to x-axis labels.
  The default maps facility/floor/room (and their `seq_` counterparts)
  to "Facility"/"Floor"/"Room".

- x_label, y_label:

  Axis titles.

## Value

A ggplot object with violin plots for each overlap type.
