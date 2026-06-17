# Plot max intra-cluster vs. min inter-cluster genetic distance.

Visualizes, for each cluster, the maximum genetic distance within the
cluster against the minimum genetic distance to an isolate in another
cluster. Points below the `y = x` line are clusters whose internal
diversity exceeds their separation from the nearest other cluster.

## Usage

``` r
plot_intra_vs_inter_cluster_distance(clusters, snp_dist)
```

## Arguments

- clusters:

  A vector named by sequence IDs giving the cluster each sequence
  belongs to. Singletons are dropped; their sequences still count as
  unclustered isolates.

- snp_dist:

  A matrix of SNP distances between isolates. Its row/column names
  define the full universe of isolates, including those not assigned to
  any cluster.

## Value

A ggplot object, suitable for saving with `ggsave()`.
