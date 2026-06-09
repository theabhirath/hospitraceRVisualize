# Create standalone legend for trace phylogeny plots

Generates a combined legend plot for trace locations, surveillance
results, and patient categories that can be saved separately from the
main plot.

## Usage

``` r
plot_trace_phylo_legend(
  trace_colors = paletteer_d("ggthemes::Classic_Cyclic"),
  trace_labels = NULL,
  surv_colors = c(neg = "blue", pos_nonclust = "red", pos_clust = "gold"),
  label_colors = c(index = "red", `multiply-colonized-index` = "darkred", `weak-index` =
    "orange", convert = "black", `adm-pos` = "forestgreen", `adm-pos-convert` = "blue",
    `secondary-convert` = "grey", other = "purple"),
  include_trace = TRUE,
  include_surv = TRUE,
  include_patient = TRUE,
  title_size = 10,
  text_size = 8,
  key_size = 0.5,
  layout = c("horizontal", "vertical")
)
```

## Arguments

- trace_colors:

  Color palette for trace location values. Should match the colors used
  in
  [`plot_trace_phylo_tree()`](https://theabhirath.github.io/hospitraceRVisualize/reference/plot_trace_phylo_tree.md).

- trace_labels:

  Character vector of labels for each trace location. Length should
  match `trace_colors`.

- surv_colors:

  Named character vector of colors for surveillance types. Default
  matches
  [`plot_trace_phylo_tree()`](https://theabhirath.github.io/hospitraceRVisualize/reference/plot_trace_phylo_tree.md)
  defaults.

- label_colors:

  Named character vector of colors for patient category labels. Default
  matches
  [`plot_trace_phylo_tree()`](https://theabhirath.github.io/hospitraceRVisualize/reference/plot_trace_phylo_tree.md)
  defaults.

- include_trace:

  Logical. Include trace location legend. Default TRUE.

- include_surv:

  Logical. Include surveillance legend. Default TRUE.

- include_patient:

  Logical. Include patient category legend. Default TRUE.

- title_size:

  Numeric. Font size for legend titles. Default 10.

- text_size:

  Numeric. Font size for legend text. Default 8.

- key_size:

  Numeric. Size of legend keys in cm. Default 0.5.

- layout:

  Character. Layout direction for combining legends. `"horizontal"`
  places legends side-by-side (default), `"vertical"` stacks them one
  below the other.

## Value

A ggplot object containing only the legend(s), suitable for saving with
`ggsave()`.
