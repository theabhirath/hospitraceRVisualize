#' Jittered scatter of two per-cluster genetic-distance measures with a y = x reference line.
#'
#' @param x,y Numeric vectors of the per-cluster distances to plot on each axis.
#' @param xlab,ylab Axis labels.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_jitter annotate coord_cartesian labs theme_minimal
#' @importFrom rlang .data
#' @keywords internal
plot_genetic_context_scatter <- function(x, y, xlab, ylab) {
    # y = x reference line, drawn only up to the smaller of the two maxima
    diag_max <- min(max(x, na.rm = TRUE), max(y, na.rm = TRUE))

    ggplot(data.frame(x = x, y = y), aes(x = .data$x, y = .data$y)) +
        annotate("segment", x = 0, y = 0, xend = diag_max, yend = diag_max) +
        geom_jitter(width = 0.2, height = 0.2, shape = 1) +
        coord_cartesian(
            xlim = c(0, max(x, na.rm = TRUE) + 1),
            ylim = c(0, max(y, na.rm = TRUE) + 1)
        ) +
        labs(x = xlab, y = ylab) +
        theme_minimal()
}

#' Per-cluster maximum intra-cluster SNP distance.
#'
#' Drops singletons and computes, for each remaining cluster, the maximum within-cluster
#' distance via [hospitraceR::cluster_pairwise_distances()]. Clusters are visited in the same
#' `sort(unique(...))` order as [hospitraceR::cluster_inter_distances()], so the result lines up
#' row-for-row with that matrix and the two can be plotted directly against each other.
#'
#' @param clusters A vector named by sequence IDs giving the cluster each sequence belongs to.
#' @param snp_dist A matrix of SNP distances between isolates.
#'
#' @return A named numeric vector of max intra-cluster distances, one per non-singleton cluster.
#'
#' @importFrom hospitraceR remove_singleton_clusters cluster_pairwise_distances
#' @keywords internal
cluster_max_intra_distances <- function(clusters, snp_dist) {
    clusters <- remove_singleton_clusters(clusters)
    vapply(
        sort(unique(clusters)),
        function(cl) {
            cluster_pairwise_distances(
                names(clusters[clusters == cl]),
                snp_dist
            )[["max_genetic_distance"]]
        },
        numeric(1)
    )
}

#' Plot max intra-cluster vs. min inter-cluster genetic distance.
#'
#' Visualizes, for each cluster, the maximum genetic distance within the cluster against the
#' minimum genetic distance to an isolate in another cluster. Points below the `y = x` line are
#' clusters whose internal diversity exceeds their separation from the nearest other cluster.
#'
#' @param clusters A vector named by sequence IDs giving the cluster each sequence belongs to.
#'   Singletons are dropped; their sequences still count as unclustered isolates.
#' @param snp_dist A matrix of SNP distances between isolates. Its row/column names define the full
#'   universe of isolates, including those not assigned to any cluster.
#'
#' @return A ggplot object, suitable for saving with `ggsave()`.
#'
#' @importFrom hospitraceR cluster_inter_distances
#' @export
plot_intra_vs_inter_cluster_distance <- function(clusters, snp_dist) {
    plot_genetic_context_scatter(
        cluster_max_intra_distances(clusters, snp_dist),
        cluster_inter_distances(clusters, snp_dist)[, "min_inter_cluster"],
        "Max genetic distance within cluster",
        "Min genetic distance to another cluster"
    )
}

#' Plot max intra-cluster vs. min inter-isolate genetic distance.
#'
#' Visualizes, for each cluster, the maximum genetic distance within the cluster against the
#' minimum genetic distance to any other isolate (including isolates not assigned to a cluster).
#' Points below the `y = x` line are clusters whose internal diversity exceeds their separation
#' from the nearest neighbouring isolate.
#'
#' @param clusters A vector named by sequence IDs giving the cluster each sequence belongs to.
#'   Singletons are dropped; their sequences still count as unclustered isolates.
#' @param snp_dist A matrix of SNP distances between isolates. Its row/column names define the full
#'   universe of isolates, including those not assigned to any cluster.
#'
#' @return A ggplot object, suitable for saving with `ggsave()`.
#'
#' @importFrom hospitraceR cluster_inter_distances
#' @export
plot_intra_vs_inter_isolate_distance <- function(clusters, snp_dist) {
    plot_genetic_context_scatter(
        cluster_max_intra_distances(clusters, snp_dist),
        cluster_inter_distances(clusters, snp_dist)[, "min_inter_isolate"],
        "Max genetic distance within cluster",
        "Min genetic distance to another isolate"
    )
}
