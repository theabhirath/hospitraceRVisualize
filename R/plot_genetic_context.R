#' Plot a jittered scatter of two per-cluster genetic-distance measures with a y = x reference line
#'
#' @param x,y Numeric vectors of per-cluster distances for each axis.
#' @param xlab,ylab Axis labels.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_jitter annotate coord_cartesian labs theme_minimal
#' @importFrom rlang .data
#' @keywords internal
plot_genetic_context_scatter <- function(x, y, xlab, ylab) {
    # y = x line drawn only up to the smaller of the two maxima
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

#' Compute the per-isolate genetic distance needed to cluster each isolate
#'
#' Drops singletons and returns, per remaining cluster, one row per isolate giving the genetic
#' distance an SNV-threshold clustering would need to place that isolate in the cluster: its
#' distances to the other members aggregated by the linkage rule (nearest member for single,
#' farthest for complete). Clusters are visited in `sort(unique(...))` order so factor levels
#' come out sorted.
#'
#' @param clusters Vector named by sequence IDs giving each sequence's cluster.
#' @param snp_dist Matrix of SNP distances between isolates.
#' @param linkage Linkage rule, `"single"` or `"complete"`.
#'
#' @return A data frame with a `cluster` factor column and a numeric `distance` column.
#'
#' @importFrom hospitraceR remove_singleton_clusters
#' @keywords internal
cluster_intra_distances_long <- function(
    clusters,
    snp_dist,
    linkage = c("single", "complete")
) {
    linkage <- match.arg(linkage)
    aggregate_distance <- switch(linkage, single = min, complete = max)

    clusters <- remove_singleton_clusters(clusters)
    cluster_ids <- sort(unique(clusters))

    rows <- lapply(cluster_ids, function(cl) {
        cluster_seqs <- names(clusters[clusters == cl])
        intra <- snp_dist[cluster_seqs, cluster_seqs, drop = FALSE]
        distances <- vapply(
            cluster_seqs,
            function(s) aggregate_distance(intra[s, setdiff(cluster_seqs, s)]),
            numeric(1)
        )
        data.frame(cluster = cl, distance = distances)
    })

    out <- do.call(rbind, rows)
    out$cluster <- factor(out$cluster, levels = cluster_ids)
    out
}

#' Plot the per-isolate genetic distance needed to cluster each isolate, by cluster
#'
#' For every non-singleton cluster, plots the genetic distance an SNV-threshold clustering would
#' need to place each isolate in the cluster (clusters on x, distance on y, one jittered point per
#' isolate). Single linkage uses each isolate's nearest cluster-mate; `"complete"` uses the
#' farthest. A boxplot summarizes each cluster.
#'
#' @param clusters Vector named by sequence IDs giving each sequence's cluster; singletons dropped.
#' @param snp_dist Matrix of SNP distances between isolates.
#' @param linkage Linkage rule, `"single"` or `"complete"`.
#'
#' @return A ggplot object, suitable for saving with `ggsave()`.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter labs theme_minimal
#' @importFrom rlang .data
#' @export
plot_genetic_distance_by_cluster <- function(
    clusters,
    snp_dist,
    linkage = c("single", "complete")
) {
    linkage <- match.arg(linkage)
    dist_df <- cluster_intra_distances_long(clusters, snp_dist, linkage)

    ggplot(dist_df, aes(x = .data$cluster, y = .data$distance)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, height = 0, shape = 1) +
        labs(x = "Cluster", y = "Genetic distance needed to cluster isolate") +
        theme_minimal()
}

#' Compute each cluster's maximum intra-cluster SNP distance
#'
#' Drops singletons and computes each remaining cluster's maximum within-cluster distance via
#' [hospitraceR::cluster_pairwise_distances()]. Clusters are visited in the same `sort(unique(...))`
#' order as [hospitraceR::cluster_inter_distances()], so the result lines up row-for-row with that
#' matrix and the two can be plotted directly against each other.
#'
#' @param clusters Vector named by sequence IDs giving each sequence's cluster.
#' @param snp_dist Matrix of SNP distances between isolates.
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

#' Plot max intra-cluster vs. min inter-cluster genetic distance
#'
#' Plots each cluster's maximum within-cluster distance against the minimum distance to an isolate
#' in another cluster. Points below `y = x` are clusters whose internal diversity exceeds their
#' separation from the nearest other cluster.
#'
#' @param clusters Vector named by sequence IDs giving each sequence's cluster; singletons dropped
#'   but their sequences still count as unclustered isolates.
#' @param snp_dist Matrix of SNP distances between isolates; row and column names define the full
#'   isolate universe, including unclustered isolates.
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

#' Plot max intra-cluster vs. min inter-isolate genetic distance
#'
#' Plots each cluster's maximum within-cluster distance against the minimum distance to any other
#' isolate (including unclustered ones). Points below `y = x` are clusters whose internal diversity
#' exceeds their separation from the nearest neighboring isolate.
#'
#' @param clusters Vector named by sequence IDs giving each sequence's cluster; singletons dropped
#'   but their sequences still count as unclustered isolates.
#' @param snp_dist Matrix of SNP distances between isolates; row and column names define the full
#'   isolate universe, including unclustered isolates.
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
