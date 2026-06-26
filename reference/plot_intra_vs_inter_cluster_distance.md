# Plot max intra-cluster vs. min inter-cluster genetic distance

Plots each cluster's maximum within-cluster distance against the minimum
distance to an isolate in another cluster. Points below `y = x` are
clusters whose internal diversity exceeds their separation from the
nearest other cluster.

## Usage

``` r
plot_intra_vs_inter_cluster_distance(clusters, snp_dist)
```

## Arguments

- clusters:

  Vector named by sequence IDs giving each sequence's cluster;
  singletons dropped but their sequences still count as unclustered
  isolates.

- snp_dist:

  Matrix of SNP distances between isolates; row and column names define
  the full isolate universe, including unclustered isolates.

## Value

A ggplot object, suitable for saving with `ggsave()`.
