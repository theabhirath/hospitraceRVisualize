# Plots a phylogeny with presence/absence trace heatmap and surveillance cell overlay

Creates a combined visualization of a phylogenetic tree with a patient
presence/absence trace heatmap. Surveillance data is encoded directly
into cell colors (not dots). Use dark gray for presence, white for
absence, and surveillance colors for results.

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

  A phylogenetic tree object of class `phylo`.

- isolate_lookup:

  A data frame from `get_isolate_lookup()` containing columns:
  isolate_id, patient_id, date, cluster, adm_pos, prev_surv.

- trace_data:

  A matrix with patient IDs as row names and dates (numeric) as column
  names. Values should be binary (0=absent, 1=present).

- surv_df:

  A data frame with surveillance data containing columns: patient_id,
  genome_id, surv_date, result (0/1 for negative/positive).

- cluster_filter:

  Numeric vector of cluster IDs to include, or NULL for all clusters.

- presence_color:

  Color for presence cells. Default "gray85" (light gray).

- surv_colors:

  Named character vector of colors for surveillance types: neg
  (negative) and pos (positive).

## Value

A ggplot object with the combined tree, heatmap, and surveillance
overlay.
