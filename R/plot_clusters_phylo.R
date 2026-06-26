#' Map clusters onto a phylogenetic tree
#'
#' Plots a phylogenetic tree with tips colored by cluster. The minimal base plot:
#' takes only the tree and cluster assignments, plus an optional sequence-to-patient
#' mapping for tip labels.
#'
#' @param tree A phylogenetic tree of class `phylo`.
#' @param clusters Vector named by sequence IDs giving each sequence's cluster.
#' @param seq2pt Optional named vector mapping sequence IDs to patient IDs; when given,
#'               adds patient tip labels.
#'
#' @return A `ggtree` object with clusters visualized on the tree.
#'
#' @importFrom ggtree ggtree geom_tippoint geom_tiplab geom_rootedge
#' @importFrom hues iwanthue
#' @importFrom rlang .data
#' @importFrom dplyr left_join
#' @importFrom ggplot2 aes scale_color_manual theme element_text unit guides guide_legend margin
#' @importFrom stats setNames
#' @export
plot_clusters_phylo <- function(
    tree,
    clusters,
    seq2pt = NULL
) {
    parents <- unique(tree$edge[, 1])
    children <- unique(tree$edge[, 2])
    root_node <- setdiff(parents, children)

    # Shorten the root clade's branch so the rootedge stays proportionate
    max_branch_length <- max(tree$edge.length)
    tree$edge.length[tree$edge[, 1] == root_node] <- max_branch_length * 0.05

    tree <- ggtree(tree)

    cluster_df <- data.frame(
        isolate = names(clusters),
        clust_id = factor(clusters)
    )

    cluster_colors <- setNames(
        iwanthue(length(unique(cluster_df$clust_id))),
        levels(cluster_df$clust_id)
    )

    tree$data <- tree$data |> left_join(cluster_df, by = c("label" = "isolate"))

    p <- tree +
        geom_tippoint(aes(color = .data$clust_id), size = 2, alpha = 0.8) +
        geom_rootedge(rootedge = max_branch_length * 0.01)

    if (!is.null(seq2pt)) {
        p <- p +
            geom_tiplab(
                aes(label = paste0(label, " (", seq2pt[label], ")")),
                size = 2,
                offset = 0.5
            )
    }

    # Scale legend sizing down as the cluster count grows so it stays compact
    n_clusters <- length(cluster_colors)
    legend_ncol <- max(1L, ceiling(n_clusters / 20))
    tier <- if (n_clusters > 40) {
        1L
    } else if (n_clusters > 20) {
        2L
    } else {
        3L
    }
    legend_text_size <- c(6, 7, 9)[tier]
    legend_key_size <- c(0.25, 0.3, 0.4)[tier]
    legend_point_size <- c(2.5, 3, 4)[tier]

    p +
        scale_color_manual(values = cluster_colors) +
        guides(
            color = guide_legend(
                title = "Cluster ID",
                ncol = legend_ncol,
                byrow = TRUE,
                override.aes = list(size = legend_point_size, shape = 15)
            )
        ) +
        theme(
            plot.title.position = "plot",
            plot.title = element_text(hjust = 0.5),
            legend.position = "right",
            legend.key.height = unit(legend_key_size, "cm"),
            legend.key.width = unit(legend_key_size, "cm"),
            legend.text = element_text(size = legend_text_size),
            legend.title = element_text(size = legend_text_size + 2, face = "bold"),
            legend.spacing.x = unit(0.15, "cm"),
            legend.spacing.y = unit(0.05, "cm"),
            legend.margin = margin(2, 2, 2, 2)
        )
}
