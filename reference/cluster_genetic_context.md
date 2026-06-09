# Produces summary plots regarding genetic distances within and between transmission clusters.

Produces summary plots regarding genetic distances within and between
transmission clusters.

## Usage

``` r
cluster_genetic_context(clusters, seq2pt, ip_seqs, snp_dist, prefix)
```

## Arguments

- clusters:

  A vector named by sequence IDs with values indicating the cluster
  (subtree) each sequence belongs to.

- seq2pt:

  A named vector mapping sequence IDs to patient IDs.

- ip_seqs:

  A vector of sequence IDs corresponding to intake positive patients.

- snp_dist:

  A matrix of SNP distances between isolates.

- prefix:

  A descriptor indicating how the clusters were generated (used to name
  figure outputs).

## Value

A matrix containing intra-cluster and inter-cluster genetic distances.
