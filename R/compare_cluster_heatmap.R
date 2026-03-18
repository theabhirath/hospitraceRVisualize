#' Compares the content of clusters created with two different methods.
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
#' @importFrom ggplot2 theme element_text scale_fill_gradient
#' @export
compare_clusters <- function(clusters1, clusters2, width = 10, height = 10) {
    # Unique cluster labels
    clusters1_unique_labels <- as.character(sort(unique(clusters1)))
    clusters2_unique_labels <- as.character(sort(unique(clusters2)))

    # Create a matrix to store the overlap between clusters
    cluster_overlap <- matrix(
        nrow = length(clusters1_unique_labels),
        ncol = length(clusters2_unique_labels),
        dimnames = list(clusters1_unique_labels, clusters2_unique_labels)
    )

    # Compute cluster overlap
    for (cluster1 in clusters1_unique_labels) {
        for (cluster2 in clusters2_unique_labels) {
            cluster_overlap[cluster1, cluster2] <- length(intersect(
                names(clusters1)[clusters1 == cluster1],
                names(clusters2)[clusters2 == cluster2]
            )) /
                length(union(
                    names(clusters1)[clusters1 == cluster1],
                    names(clusters2)[clusters2 == cluster2]
                ))
        }
    }

    # get clusters that map exactly 1-1
    one_to_one_mapping <- rowSums(cluster_overlap == 1) == 1
    # get isolates that are in one-to-one mapping
    isolates_in_one_to_one_mapping <- names(clusters1)[one_to_one_mapping]
    # percentage of isolates that are in one-to-one mapping
    overlap_percentage <- length(isolates_in_one_to_one_mapping) /
        length(clusters1)

    # Return the ggheatmap plot
    ggheatmap(cluster_overlap, width = width, height = height) +
        ggtitle(paste0(
            "Comparison of clusters. Isolate overlap: ",
            round(overlap_percentage * 100, 2),
            "%"
        ))
}
