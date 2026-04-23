#' Maps clusters onto a phylogenetic tree and visualizes them
#'
#'
#' @param tree A phylogenetic tree object of class `phylo`.
#' @param clusters A vector named by sequence IDs with values indicating the cluster (subtree) each sequence belongs to.
#' @param seq2pt OPTIONAL. A named vector mapping sequence IDs to patient IDs. Default is NULL – no patient labels will
#'               be added to the tree.
#' @param patient_label OPTIONAL. A logical value indicating whether to add patient labels to the tree. Default is FALSE
#'                      – no patient labels will be added to the tree.
#' @param dna_var OPTIONAL. A DNA alignment object of class `DNAbin`. Default is NULL – no variants will be plotted.
#' @param convert_status OPTIONAL. A logical value indicating whether to add convert status to the tree. Default is
#'                       FALSE – no convert status will be added to the tree.
#' @param ip_seqs OPTIONAL. A vector of sequence IDs which correspond to admission positive sequences. Default is
#'                NULL – no admission positive sequences will be plotted.
#' @param dates OPTIONAL. A named vector mapping sequence IDs to dates. Default is NULL – no dates will be plotted.
#' @param pt_trace OPTIONAL. A matrix with patient IDs as row names and dates as column names. Default is NULL –
#'                 no patient trace will be plotted.
#'
#' @return A `ggtree` object with clusters visualized on the tree.
#'
#' @importFrom ggtree ggtree geom_tippoint geom_facet geom_tiplab geom_rootedge scaleClade
#' @importFrom hues iwanthue
#' @importFrom rlang .data
#' @importFrom dplyr left_join
#' @importFrom ggplot2 aes scale_color_manual theme element_text unit guides guide_legend margin
#' @importFrom stats setNames
#' @importFrom ape dist.dna as.DNAbin
#' @export
plot_clusters_phylo <- function(
    tree,
    clusters,
    seq2pt = NULL,
    dna_var = NULL,
    patient_label = FALSE,
    convert_status = FALSE,
    ip_seqs = NULL,
    dates = NULL,
    pt_trace = NULL
) {
    # ensure that epi data is provided if convert status is TRUE
    if (convert_status) {
        if (is.null(pt_trace) || is.null(ip_seqs) || is.null(dates) || is.null(seq2pt)) {
            stop(paste(
                "pt_trace, ip_seqs, dates, and seq2pt must be provided if convert_status",
                "is TRUE, since convert status requires epi data to be provided."
            ))
        }
    }
    parents <- unique(tree$edge[, 1])
    children <- unique(tree$edge[, 2])
    root_node <- setdiff(parents, children)

    # find max branch length
    max_branch_length <- max(tree$edge.length)
    # reduce branch length of root node clade to 50% of max branch length
    tree$edge.length[tree$edge[, 1] == root_node] <- max_branch_length * 0.05

    # convert phylo object to ggtree object
    tree <- ggtree(tree)

    # Format the clusters into a dataframe
    cluster_df <- data.frame(
        isolate = names(clusters),
        clust_id = factor(clusters)
    )

    # Use iwanthue to generate colors for clusters
    cluster_colors <- setNames(
        iwanthue(length(unique(cluster_df$clust_id))),
        levels(cluster_df$clust_id)
    )

    # identify convert isolates
    if (convert_status) {
        is_convert_isolate <- sapply(names(clusters), function(id) {
            if (is.na(seq2pt[id])) {
                print(paste0("Sequence ", id, " has no patient label."))
                return(FALSE)
            }
            row_vals <- pt_trace[as.character(seq2pt[id]), ]
            i_negative <- which(row_vals %in% c(1.25))
            i_positive <- which(row_vals %in% c(1.5))
            if (length(i_negative) == 0 || length(i_positive) == 0) {
                return(FALSE)
            }
            is_convert <- min(i_negative) < min(i_positive)
            trace_date <- as.numeric(colnames(pt_trace)[min(i_positive)])
            iso_date <- dates[id]
            (iso_date - trace_date < 7) && is_convert && !(id %in% ip_seqs)
        })
        # index isolates are those which are in ip_seqs
        is_index_isolate <- names(clusters) %in% ip_seqs
        # add convert class information to cluster_df
        convert_class <- rep("Convert patient", nrow(cluster_df))
        convert_class[is_index_isolate] <- "Index patient"
        cluster_df$convert_class <- convert_class
    }

    # add cluster information to the tree
    tree$data <- tree$data |> left_join(cluster_df, by = c("label" = "isolate"))

    # build the plot
    # if convert_status is true, use shape of nodes to indicate convert class
    if (convert_status) {
        p <- tree +
            geom_tippoint(
                aes(color = .data$clust_id, shape = .data$convert_class),
                size = 2,
                alpha = 0.8
            ) +
            geom_rootedge(rootedge = max_branch_length * 0.01)
    }

    # if dna_var is provided, plot the variants
    if (!is.null(dna_var)) {
        # coerce DNAbin to character matrix properly using ape's method
        dna_mat <- as.character(as.matrix(dna_var))
        # ensure column names (position IDs) exist for consistent reordering
        if (is.null(colnames(dna_mat))) {
            colnames(dna_mat) <- paste0("Pos", seq_len(ncol(dna_mat)))
        }
        # hierarchical clustering of variant columns using DNA distance
        # rows = positions, columns = samples for distance over positions
        dna_for_clust <- t(dna_mat)
        # ensure row names (position identifiers) mirror position column names
        rownames(dna_for_clust) <- colnames(dna_mat)
        # filter to informative positions, ignoring unknowns ('N' and '-')
        unknowns <- c("-", "N", "n")
        keep_rows <- apply(dna_for_clust, 1, function(row_vals) {
            row_no_unknown <- row_vals[!(row_vals %in% unknowns)]
            uniq <- unique(row_vals)
            uniq_no_unknown <- setdiff(uniq, unknowns)
            if (length(uniq_no_unknown) == 0) {
                FALSE # all N/- → drop
            } else if (length(uniq_no_unknown) == 1 && length(uniq) == 1) {
                FALSE # all same → drop
            } else if (sum(table(row_no_unknown) >= 2) < 2) {
                FALSE # the column has at least two variants but minor variants only appeared once → drop
            } else {
                TRUE # informative site
            }
        })
        dna_filtered <- dna_for_clust[keep_rows, , drop = FALSE]
        if (nrow(dna_filtered) >= 2) {
            clustered_pos_order <- tryCatch(
                {
                    dna_dist <- dist.dna(
                        as.DNAbin(dna_filtered),
                        pairwise.deletion = TRUE,
                        model = "N"
                    )
                    hclust(dna_dist, method = "average")$order
                },
                error = function(e) {
                    seq_len(nrow(dna_filtered))
                }
            )
            ordered_pos_labels <- rownames(dna_filtered)[clustered_pos_order]
            # reorder original columns to clustered order; append any dropped/filtered cols at the end
            remaining_cols <- setdiff(colnames(dna_mat), ordered_pos_labels)
            cols_ordered <- c(ordered_pos_labels, remaining_cols)
            dna_mat <- dna_mat[, cols_ordered, drop = FALSE]
        }
        out_group <- 1
        # calculate 0/1/2 variant codes relative to the specified outgroup row in dna_var
        # 0: no variant, 1: variant, 2: isolate base unknown ("-","N","n")
        out_bases <- dna_mat[out_group, , drop = TRUE]
        # start with 0s
        variant_code_mat <- matrix(
            0L,
            nrow = nrow(dna_mat),
            ncol = ncol(dna_mat),
            dimnames = list(rownames(dna_mat), colnames(dna_mat))
        )
        # mark variants as 1 where isolate base differs from outgroup base
        diff_mat <- dna_mat !=
            matrix(
                out_bases,
                nrow = nrow(dna_mat),
                ncol = ncol(dna_mat),
                byrow = TRUE
            )
        variant_code_mat[diff_mat] <- 1L
        # override to 2 where the isolate base is unknown (do not consider outgroup unknowns)
        is_unknown_isolate <- dna_mat %in% unknowns
        variant_code_mat[is_unknown_isolate] <- 2L
        # filter positions: if all variant, or all unknown, or all non-variant, drop
        # check this for all non-reference genomes
        keep_cols <- apply(variant_code_mat[-1, ], 2, function(col_vals) {
            if (all(col_vals %in% unknowns)) {
                FALSE
            } else if (all(col_vals == col_vals[1])) {
                FALSE
            } else {
                TRUE
            }
        })
        variant_code_mat <- variant_code_mat[, keep_cols, drop = FALSE]

        # build long-format data directly from matrix (avoid data.frame conversion issues)
        row_labels <- if (!is.null(rownames(variant_code_mat))) {
            rownames(variant_code_mat)
        } else {
            as.character(seq_len(nrow(variant_code_mat)))
        }
        col_labels <- colnames(variant_code_mat)
        long_df <- data.frame(
            label = rep(row_labels, times = length(col_labels)),
            variant = rep(col_labels, each = length(row_labels)),
            value = factor(as.vector(variant_code_mat), levels = c(0L, 1L, 2L)),
            stringsAsFactors = FALSE
        )
        # derive tip order strictly by y-position to match plotted order
        tip_order <- with(tree$data[tree$data$isTip, ], label[order(y)])
        # filter to only labels present in tree tips
        long_df <- long_df[long_df$label %in% tip_order, ]
        long_df$label <- factor(long_df$label, levels = tip_order)
        # numeric y aligned to tip order to satisfy continuous scale inside facet panel
        label_index_map <- setNames(seq_along(tip_order), tip_order)
        long_df$label_index <- label_index_map[as.character(long_df$label)]
        # map variant IDs to a numeric index to avoid continuous scale errors in facet panel
        variant_levels <- unique(col_labels)
        variant_index_map <- setNames(seq_along(variant_levels), variant_levels)
        long_df$variant_index <- as.numeric(variant_index_map[long_df$variant])

        # this has to be plotted as a heatmap aligned to the tree tips
        if (convert_status) {
            p <- tree +
                geom_tippoint(
                    aes(color = .data$clust_id, shape = .data$convert_class),
                    size = 2,
                    alpha = 0.8
                ) +
                geom_rootedge(rootedge = max_branch_length * 0.01) +
                geom_facet(
                    panel = "Variants",
                    data = long_df,
                    mapping = aes(x = variant_index, y = label_index, fill = value),
                    geom = geom_tile
                ) +
                scale_fill_manual(
                    name = "Variant",
                    values = c("0" = "lightgray", "1" = "red", "2" = "black"),
                    drop = FALSE
                )
        } else {
            p <- tree +
                geom_tippoint(aes(color = .data$clust_id), size = 2, alpha = 0.8) +
                geom_rootedge(rootedge = max_branch_length * 0.01) +
                geom_facet(
                    panel = "Variants",
                    data = long_df,
                    mapping = aes(x = variant_index, y = label_index, fill = value),
                    geom = geom_tile
                ) +
                scale_fill_manual(
                    name = "Variant",
                    values = c("0" = "lightgray", "1" = "red", "2" = "black"),
                    drop = FALSE
                )
        }
    } else {
        p <- tree +
            geom_tippoint(aes(color = .data$clust_id), size = 2, alpha = 0.8) +
            geom_rootedge(rootedge = max_branch_length * 0.01)
    }

    # Add patient labels if requested
    if (!is.null(seq2pt) && patient_label) {
        p <- p +
            geom_tiplab(
                aes(label = paste0(label, " (", seq2pt[label], ")")),
                size = 2,
                offset = 0.5
            )
    }

    # Adapt legend layout to number of clusters so it stays compact
    n_clusters <- length(cluster_colors)
    legend_ncol <- max(1L, ceiling(n_clusters / 20))
    legend_text_size <- if (n_clusters > 40) 6 else if (n_clusters > 20) 7 else 9
    legend_key_size <- if (n_clusters > 40) 0.25 else if (n_clusters > 20) 0.3 else 0.4
    legend_point_size <- if (n_clusters > 40) 2.5 else if (n_clusters > 20) 3 else 4

    # Return the final plot
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
