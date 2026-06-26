# Plot a phylogeny with a location trace heatmap

Creates a combined visualization of a phylogenetic tree with a patient
location trace heatmap. Uses ggtree for the phylogeny with continuous
line segments for the trace data (not discrete cells). Surveillance data
is overlaid as dots.

## Usage

``` r
plot_trace_phylo_tree(
  tree,
  isolate_lookup,
  trace_data,
  surv_df,
  cluster_filter = NULL,
  trace_colors = trace_color_palette(),
  surv_colors = surv_dot_colors(),
  surv_halos = surv_dot_halos(),
  surv_shapes = surv_dot_shapes(),
  clust_patient_categories = NULL,
  label_colors = c(index = "red", `multiply-colonized-index` = "darkred", `weak-index` =
    "orange", convert = "black", `adm-pos` = "forestgreen", `adm-pos-convert` = "blue",
    `secondary-convert` = "gray", `ambiguous-adm-pos` = "turquoise3", `ambiguous-convert`
    = "deeppink3", other = "purple"),
  inches_per_row = 0.2,
  row_thickness = 0.55,
  max_tree_width = 0.08,
  show_legend = TRUE
)
```

## Arguments

- tree:

  A phylogenetic tree of class `phylo`.

- isolate_lookup:

  Data frame from `get_isolate_lookup()` with columns isolate_id,
  patient_id, date, cluster, adm_pos, prev_surv.

- trace_data:

  Matrix with patient IDs as row names and numeric dates as column
  names; values are location categories (0 = absent, 1+ = locations).

- surv_df:

  Surveillance data frame with columns patient_id, genome_id, surv_date,
  result (0/1 for negative/positive).

- cluster_filter:

  Cluster IDs to include, or `NULL` for all.

- trace_colors:

  Color palette for trace location values.

- surv_colors:

  Named fill colors for surveillance types (neg, pos_clust,
  pos_nonclust).

- surv_halos:

  Named outline (halo) colors, keyed like `surv_colors`.

- surv_shapes:

  Named integer plotting shapes (filled, 21-25), keyed like
  `surv_colors`.

- clust_patient_categories:

  Named list from `cluster_patient_categorization()` mapping patient_id
  to category per cluster; `NULL` labels all patients "other".

- label_colors:

  Named colors for patient category labels; "other" is the fallback for
  any tip with no assigned category.

- inches_per_row:

  Height in inches per patient row when saving.

- row_thickness:

  Proportion of row spacing used for trace bar height.

- max_tree_width:

  Maximum tree width as a proportion of the time range; branch lengths
  are scaled to fit.

- show_legend:

  Logical. Attach the figure's own legends (location, patient status,
  surveillance); `FALSE` selects the separate-legend path.

## Value

A ggplot object with the following attributes for consistent sizing:

- `recommended_height`: Suggested height in inches

- `recommended_width`: Suggested width in inches

- `n_patients`: Number of patients in the plot

- `inches_per_row`: Height allocated per patient row

Use these when saving with `ggsave()` to maintain consistent row
heights. When `show_legend = TRUE` the figure carries its own native
legends and the returned object is a single ggplot/ggtree;
`recommended_width` is widened to fit the legend.
