#' Plots a heatmap of the Jaccard similarity between clusters created with two different methods.
#'
#' @param clusters1 A vector named by sequence IDs with values being subtrees defining the cluster.
#' @param clusters2 A vector named by the same sequence IDs as clusters1 with values being subtrees
#' defining the cluster.
#' @param width The width of the heatmap plot.
#' @param height The height of the heatmap plot.
#'
#' @return A ggplot plot object showing the overlap between clusters.
#'
#' @importFrom ggalign ggheatmap
#' @importFrom ggplot2 ggtitle
#' @importFrom hospitraceR cluster_contingency_table
#' @export
plot_jaccard_similarity_heatmap <- function(clusters1, clusters2, width = 10, height = 10) {
    # Count shared isolates between every pair of clusters. Rows are clusters
    # from clusters1, columns are clusters from clusters2.
    cont_table <- cluster_contingency_table(clusters1, clusters2)

    # Convert the shared-isolate counts to a Jaccard overlap. For clusters i and
    # j the intersection is the contingency count n_ij and the union is
    # |cluster_i| + |cluster_j| - n_ij, which the row and column marginals give.
    row_sizes <- rowSums(cont_table)
    col_sizes <- colSums(cont_table)
    union_sizes <- outer(row_sizes, col_sizes, `+`) - cont_table
    cluster_overlap <- cont_table / union_sizes

    # A clusters1 cluster maps 1-to-1 when exactly one clusters2 cluster has an
    # identical membership (overlap == 1). Report the fraction of isolates that
    # live in such cleanly-mapped clusters.
    one_to_one_mapping <- rowSums(cluster_overlap == 1) == 1
    isolates_in_one_to_one_mapping <- sum(row_sizes[one_to_one_mapping])
    overlap_percentage <- isolates_in_one_to_one_mapping / sum(row_sizes)

    # Return the ggheatmap plot
    ggheatmap(cluster_overlap, width = width, height = height) +
        ggtitle(paste0(
            "Comparison of clusters. Isolate overlap: ",
            round(overlap_percentage * 100, 2),
            "%"
        ))
}
