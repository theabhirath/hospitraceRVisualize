# Per-cluster maximum intra-cluster SNP distance.

Drops singletons and computes, for each remaining cluster, the maximum
within-cluster distance via
[`hospitraceR::cluster_pairwise_distances()`](https://theabhirath.github.io/hospitraceR/reference/cluster_pairwise_distances.html).
Clusters are visited in the same `sort(unique(...))` order as
[`hospitraceR::cluster_inter_distances()`](https://theabhirath.github.io/hospitraceR/reference/cluster_inter_distances.html),
so the result lines up row-for-row with that matrix and the two can be
plotted directly against each other.

## Usage

``` r
cluster_max_intra_distances(clusters, snp_dist)
```

## Arguments

- clusters:

  A vector named by sequence IDs giving the cluster each sequence
  belongs to.

- snp_dist:

  A matrix of SNP distances between isolates.

## Value

A named numeric vector of max intra-cluster distances, one per
non-singleton cluster.
