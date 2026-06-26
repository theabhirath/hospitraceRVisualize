# Plot the per-isolate genetic distance needed to cluster each isolate, by cluster

For every non-singleton cluster, plots the genetic distance an
SNV-threshold clustering would need to place each isolate in the cluster
(clusters on x, distance on y, one jittered point per isolate). Single
linkage uses each isolate's nearest cluster-mate; `"complete"` uses the
farthest. A boxplot summarizes each cluster.

## Usage

``` r
plot_genetic_distance_by_cluster(
  clusters,
  snp_dist,
  linkage = c("single", "complete")
)
```

## Arguments

- clusters:

  Vector named by sequence IDs giving each sequence's cluster;
  singletons dropped.

- snp_dist:

  Matrix of SNP distances between isolates.

- linkage:

  Linkage rule, `"single"` or `"complete"`.

## Value

A ggplot object, suitable for saving with `ggsave()`.
