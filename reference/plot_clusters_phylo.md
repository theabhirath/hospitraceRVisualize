# Map clusters onto a phylogenetic tree

Plots a phylogenetic tree with tips colored by cluster. The minimal base
plot: takes only the tree and cluster assignments, plus an optional
sequence-to-patient mapping for tip labels.

## Usage

``` r
plot_clusters_phylo(tree, clusters, seq2pt = NULL)
```

## Arguments

- tree:

  A phylogenetic tree of class `phylo`.

- clusters:

  Vector named by sequence IDs giving each sequence's cluster.

- seq2pt:

  Optional named vector mapping sequence IDs to patient IDs; when given,
  adds patient tip labels.

## Value

A `ggtree` object with clusters visualized on the tree.
