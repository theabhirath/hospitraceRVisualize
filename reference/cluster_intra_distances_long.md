# Compute the per-isolate genetic distance needed to cluster each isolate

Drops singletons and returns, per remaining cluster, one row per isolate
giving the genetic distance an SNV-threshold clustering would need to
place that isolate in the cluster: its distances to the other members
aggregated by the linkage rule (nearest member for single, farthest for
complete). Clusters are visited in `sort(unique(...))` order so factor
levels come out sorted.

## Usage

``` r
cluster_intra_distances_long(
  clusters,
  snp_dist,
  linkage = c("single", "complete")
)
```

## Arguments

- clusters:

  Vector named by sequence IDs giving each sequence's cluster.

- snp_dist:

  Matrix of SNP distances between isolates.

- linkage:

  Linkage rule, `"single"` or `"complete"`.

## Value

A data frame with a `cluster` factor column and a numeric `distance`
column.
