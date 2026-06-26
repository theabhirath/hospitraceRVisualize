# Create a standalone legend for trace phylogeny plots

Generates a combined legend plot for trace locations, surveillance
results, and patient categories that can be saved separately from the
main plot.

## Usage

``` r
plot_trace_phylo_legend(
  trace_colors = trace_color_palette(),
  trace_labels = NULL,
  surv_colors = surv_dot_colors(),
  surv_shapes = surv_dot_shapes(),
  label_colors = c(index = "red", `multiply-colonized-index` = "darkred", `weak-index` =
    "orange", convert = "black", `adm-pos` = "forestgreen", `adm-pos-convert` = "blue",
    `secondary-convert` = "gray", `ambiguous-adm-pos` = "turquoise3", `ambiguous-convert`
    = "deeppink3"),
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

  Color palette for trace location values; should match
  [`plot_trace_phylo_tree()`](https://theabhirath.github.io/hospitraceRVisualize/reference/plot_trace_phylo_tree.md).

- trace_labels:

  Character vector of labels per trace location; length should match
  `trace_colors`.

- surv_colors:

  Named colors for surveillance types; should match
  [`plot_trace_phylo_tree()`](https://theabhirath.github.io/hospitraceRVisualize/reference/plot_trace_phylo_tree.md).

- surv_shapes:

  Named integer plotting shapes for surveillance types; should match
  [`plot_trace_phylo_tree()`](https://theabhirath.github.io/hospitraceRVisualize/reference/plot_trace_phylo_tree.md).

- label_colors:

  Named colors for patient category labels; should match
  [`plot_trace_phylo_tree()`](https://theabhirath.github.io/hospitraceRVisualize/reference/plot_trace_phylo_tree.md).

- include_trace:

  Logical. Include the trace location legend.

- include_surv:

  Logical. Include the surveillance legend.

- include_patient:

  Logical. Include the patient category legend.

- title_size:

  Numeric. Font size for legend titles.

- text_size:

  Numeric. Font size for legend text.

- key_size:

  Numeric. Size of legend keys in cm.

- layout:

  Character. `"horizontal"` places legends side-by-side, `"vertical"`
  stacks them.

## Value

A ggplot object containing only the legend(s), suitable for saving with
`ggsave()`.
