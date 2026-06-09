# Plots a heatmap of the Jaccard similarity between clusters created with two different methods.

Plots a heatmap of the Jaccard similarity between clusters created with
two different methods.

## Usage

``` r
plot_jaccard_similarity_heatmap(clusters1, clusters2, width = 10, height = 10)
```

## Arguments

- clusters1:

  A vector named by sequence IDs with values being subtrees defining the
  cluster.

- clusters2:

  A vector named by the same sequence IDs as clusters1 with values being
  subtrees defining the cluster.

- width:

  The width of the heatmap plot.

- height:

  The height of the heatmap plot.

## Value

A ggplot plot object showing the overlap between clusters.
