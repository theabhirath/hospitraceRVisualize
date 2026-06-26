#' Plot a heatmap of Jaccard similarity between two clusterings
#'
#' Plots the pairwise Jaccard overlap between the clusters of two clusterings of the same
#' sequences, with the overall fraction of exactly-matching clusters in the title.
#'
#' @param clusters1 Vector named by sequence IDs giving each sequence's cluster.
#' @param clusters2 Vector named by the same sequence IDs as `clusters1` giving each sequence's
#'   cluster.
#' @param width Width of the heatmap.
#' @param height Height of the heatmap.
#'
#' @return A ggplot object showing the overlap between clusters.
#'
#' @importFrom ggalign ggheatmap
#' @importFrom ggplot2 ggtitle
#' @importFrom hospitraceR cluster_contingency_table
#' @export
plot_jaccard_similarity_heatmap <- function(clusters1, clusters2, width = 10, height = 10) {
    # Shared-isolate counts: rows are clusters1 clusters, columns are clusters2 clusters
    cont_table <- cluster_contingency_table(clusters1, clusters2)

    # Jaccard overlap: intersection is the contingency count n_ij and the union is
    # |cluster_i| + |cluster_j| - n_ij, which the row and column marginals give
    row_sizes <- rowSums(cont_table)
    col_sizes <- colSums(cont_table)
    union_sizes <- outer(row_sizes, col_sizes, `+`) - cont_table
    cluster_overlap <- cont_table / union_sizes

    # calculate overall overlap percentage i.e. the fraction of isolates
    # that are in clusters which are exactly the same in both clusterings
    one_to_one_mapping <- rowSums(cluster_overlap == 1) == 1
    isolates_in_one_to_one_mapping <- sum(row_sizes[one_to_one_mapping])
    overlap_percentage <- isolates_in_one_to_one_mapping / sum(row_sizes)

    ggheatmap(cluster_overlap, width = width, height = height) +
        ggtitle(paste0(
            "Comparison of clusters. Isolate overlap: ",
            round(overlap_percentage * 100, 2),
            "%"
        ))
}
