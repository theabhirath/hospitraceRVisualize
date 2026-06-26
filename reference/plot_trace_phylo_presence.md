# Plot a phylogeny with a presence/absence trace heatmap and surveillance overlay

Creates a combined visualization of a phylogenetic tree with a patient
presence/absence trace heatmap. Surveillance data is encoded directly
into cell colors (not dots): dark gray for presence, white for absence,
and surveillance colors for results.

## Usage

``` r
plot_trace_phylo_presence(
  tree,
  isolate_lookup,
  trace_data,
  surv_df,
  cluster_filter = NULL,
  presence_color = "gray85",
  surv_colors = c(neg = "blue", pos = "red")
)
```

## Arguments

- tree:

  A phylogenetic tree of class `phylo`.

- isolate_lookup:

  Data frame from `get_isolate_lookup()` with columns isolate_id,
  patient_id, date, cluster, adm_pos, prev_surv.

- trace_data:

  Binary matrix (0/1) with patient IDs as row names and dates as column
  names.

- surv_df:

  Surveillance data frame with columns patient_id, genome_id, surv_date,
  result (0/1 for negative/positive).

- cluster_filter:

  Cluster IDs to include, or `NULL` for all.

- presence_color:

  Color for presence cells.

- surv_colors:

  Named colors for the negative and positive surveillance types.

## Value

A ggplot object combining the tree, heatmap, and surveillance overlay.
