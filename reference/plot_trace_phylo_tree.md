# Plots a phylogeny with trace heatmap using ggtree with consistent sizing

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
  trace_colors = paletteer_d("ggthemes::Classic_Cyclic"),
  surv_colors = c(neg = "blue", pos_nonclust = "red", pos_clust = "gold"),
  clust_patient_categories = NULL,
  label_colors = c(index = "red", `multiply-colonized-index` = "darkred", `weak-index` =
    "orange", convert = "black", `adm-pos` = "forestgreen", `adm-pos-convert` = "blue",
    `secondary-convert` = "grey", other = "purple"),
  inches_per_row = 0.2,
  row_thickness = 0.4,
  max_tree_width = 0.08
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
  names. Values are numeric categories (0=absent, 1+=location
  categories).

- surv_df:

  A data frame with surveillance data containing columns: patient_id,
  genome_id, surv_date, result (0/1 for negative/positive).

- cluster_filter:

  Numeric vector of cluster IDs to include, or NULL for all clusters.

- trace_colors:

  Color palette for trace location values.

- surv_colors:

  Named character vector of colors for surveillance types: neg
  (negative), pos_clust (positive cluster), pos_nonclust (positive
  non-cluster).

- clust_patient_categories:

  Named list from `cluster_patient_categorization()`. Each element is a
  named vector mapping patient_id to category. If NULL, all patients are
  labeled as "other".

- label_colors:

  Named character vector of colors for patient category labels.
  Categories from cluster_patient_categorization: index,
  multiply-colonized-index, weak-index, convert, adm-pos,
  adm-pos-convert, secondary-convert, other.

- inches_per_row:

  Height in inches allocated per patient row when saving. Default 0.2.
  Use this to control the visual thickness of rows in the output.

- row_thickness:

  Proportion of row spacing used for trace bar height. Default 0.4.

- max_tree_width:

  Maximum width for tree as proportion of time range. Default 0.08. Tree
  branch lengths are scaled to fit within this width.

## Value

A ggplot object with the following attributes for consistent sizing:

- `recommended_height`: Suggested height in inches

- `recommended_width`: Suggested width in inches

- `n_patients`: Number of patients in the plot

- `inches_per_row`: Height allocated per patient row

Use these when saving with `ggsave()` to maintain consistent row
heights.
