# Plot max intra-cluster vs. min inter-isolate genetic distance

Plots each cluster's maximum within-cluster distance against the minimum
distance to any other isolate (including unclustered ones). Points below
`y = x` are clusters whose internal diversity exceeds their separation
from the nearest neighboring isolate.

## Usage

``` r
plot_intra_vs_inter_isolate_distance(clusters, snp_dist)
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
