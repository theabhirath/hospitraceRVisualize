# Maps clusters onto a phylogenetic tree and visualizes them

Maps clusters onto a phylogenetic tree and visualizes them

## Usage

``` r
plot_clusters_phylo(
  tree,
  clusters,
  seq2pt = NULL,
  dna_var = NULL,
  patient_label = FALSE,
  convert_status = FALSE,
  ip_seqs = NULL,
  dates = NULL,
  pt_trace = NULL
)
```

## Arguments

- tree:

  A phylogenetic tree object of class `phylo`.

- clusters:

  A vector named by sequence IDs with values indicating the cluster
  (subtree) each sequence belongs to.

- seq2pt:

  OPTIONAL. A named vector mapping sequence IDs to patient IDs. Default
  is NULL – no patient labels will be added to the tree.

- dna_var:

  OPTIONAL. A DNA alignment object of class `DNAbin`. Default is NULL –
  no variants will be plotted.

- patient_label:

  OPTIONAL. A logical value indicating whether to add patient labels to
  the tree. Default is FALSE – no patient labels will be added to the
  tree.

- convert_status:

  OPTIONAL. A logical value indicating whether to add convert status to
  the tree. Default is FALSE – no convert status will be added to the
  tree.

- ip_seqs:

  OPTIONAL. A vector of sequence IDs which correspond to admission
  positive sequences. Default is NULL – no admission positive sequences
  will be plotted.

- dates:

  OPTIONAL. A named vector mapping sequence IDs to dates. Default is
  NULL – no dates will be plotted.

- pt_trace:

  OPTIONAL. A matrix with patient IDs as row names and dates as column
  names. Default is NULL – no patient trace will be plotted.

## Value

A `ggtree` object with clusters visualized on the tree.
