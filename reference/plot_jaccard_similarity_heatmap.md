# Plot a heatmap of Jaccard similarity between two clusterings

Plots the pairwise Jaccard overlap between the clusters of two
clusterings of the same sequences, with the overall fraction of
exactly-matching clusters in the title.

## Usage

``` r
plot_jaccard_similarity_heatmap(clusters1, clusters2, width = 10, height = 10)
```

## Arguments

- clusters1:

  Vector named by sequence IDs giving each sequence's cluster.

- clusters2:

  Vector named by the same sequence IDs as `clusters1` giving each
  sequence's cluster.

- width:

  Width of the heatmap.

- height:

  Height of the heatmap.

## Value

A ggplot object showing the overlap between clusters.
