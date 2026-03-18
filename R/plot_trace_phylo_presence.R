#' Plots a phylogeny with presence/absence trace heatmap and surveillance cell overlay
#'
#' Creates a combined visualization of a phylogenetic tree with a patient presence/absence
#' trace heatmap. Surveillance data is encoded directly into cell colors (not dots).
#' Use dark gray for presence, white for absence, and surveillance colors for results.
#'
#' @param tree A phylogenetic tree object of class `phylo`.
#' @param isolate_lookup A data frame from `get_isolate_lookup()` containing columns:
#'   isolate_id, patient_id, date, cluster, adm_pos, prev_surv.
#' @param trace_data A matrix with patient IDs as row names and dates (numeric) as
#'   column names. Values should be binary (0=absent, 1=present).
#' @param surv_df A data frame with surveillance data containing columns:
#'   patient_id, genome_id, surv_date, result (0/1 for negative/positive).
#' @param cluster_filter Numeric vector of cluster IDs to include, or NULL for all clusters.
#' @param presence_color Color for presence cells. Default "gray85" (light gray).
#' @param surv_colors Named character vector of colors for surveillance types:
#'   neg (negative) and pos (positive).
#'
#' @return A ggplot object with the combined tree, heatmap, and surveillance overlay.
#'
#' @importFrom ggtree ggtree geom_tiplab gheatmap vexpand
#' @importFrom ape keep.tip
#' @importFrom ggtreeExtra geom_fruit
#' @importFrom ggnewscale new_scale_fill
#' @importFrom ggplot2 aes geom_tile scale_fill_manual scale_fill_gradientn theme
#'   element_text unit
#' @importFrom hues iwanthue
#' @importFrom dplyr left_join filter mutate case_when group_by arrange slice_min
#'   ungroup pull
#' @importFrom rlang .data
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames
#' @export
plot_trace_phylo_presence <- function(
    tree,
    isolate_lookup,
    trace_data,
    surv_df,
    cluster_filter = NULL,
    presence_color = "gray85",
    surv_colors = c(neg = "blue", pos = "red")
) {
    # ========================================================================
    # 1. Filter data by cluster if specified
    # ========================================================================
    if (!is.null(cluster_filter)) {
        isolate_lookup <- isolate_lookup[
            isolate_lookup$cluster %in% cluster_filter,
        ]
    }

    # ========================================================================
    # 2. Prepare isolate ordering
    # ========================================================================
    cluster_order <- sort(unique(isolate_lookup$cluster))

    isolate_order <- isolate_lookup |>
        group_by(.data$patient_id) |>
        mutate(cluster_factor = factor(.data$cluster, levels = cluster_order)) |>
        arrange(.data$cluster_factor, .data$date) |>
        slice_min(.data$date, with_ties = FALSE) |>
        ungroup() |>
        pull(.data$isolate_id)

    pt_in_trace <- isolate_lookup$patient_id %in% row.names(trace_data)
    valid_isolates <- isolate_lookup$isolate_id[pt_in_trace]
    keep_isolates <- intersect(isolate_order, valid_isolates)

    if (length(keep_isolates) == 0) {
        stop("No isolates found with matching patients in trace_data.")
    }

    # ========================================================================
    # 3. Prepare trace data with surveillance overlay
    # ========================================================================
    isolate_patient_map <- setNames(
        isolate_lookup$patient_id,
        isolate_lookup$isolate_id
    )
    trace_sub <- trace_data[
        as.character(isolate_patient_map[keep_isolates]),
        ,
        drop = FALSE
    ]
    row.names(trace_sub) <- keep_isolates

    # Get visible patients
    visible_patients <- unique(isolate_patient_map[keep_isolates])

    # Categorize surveillance results
    surv_categorized <- surv_df |>
        filter(.data$patient_id %in% visible_patients) |>
        mutate(
            surv_type = case_when(
                .data$result == 0 ~ "neg",
                .data$result == 1 ~ "pos",
                TRUE ~ NA_character_
            )
        ) |>
        filter(!is.na(.data$surv_type))

    # Encode surveillance into trace matrix using special values:
    # 0 = absent, 1 = present, 1.25 = neg surv, 1.5 = pos surv
    trace_with_surv <- as.matrix(trace_sub)

    if (nrow(surv_categorized) > 0) {
        # Create reverse map: patient_id -> isolate_id (for isolates in trace)
        patient_to_iso <- setNames(
            keep_isolates,
            isolate_patient_map[keep_isolates]
        )

        # Filter to surveillance records with valid row/column mappings
        surv_to_encode <- surv_categorized |>
            mutate(
                iso_id = patient_to_iso[as.character(.data$patient_id)],
                col_name = as.character(.data$surv_date),
                surv_val = ifelse(.data$surv_type == "neg", 1.25, 1.5)
            ) |>
            filter(
                !is.na(.data$iso_id),
                .data$col_name %in% colnames(trace_with_surv)
            )

        # Apply all surveillance values at once using matrix indexing
        if (nrow(surv_to_encode) > 0) {
            idx <- cbind(
                match(surv_to_encode$iso_id, rownames(trace_with_surv)),
                match(surv_to_encode$col_name, colnames(trace_with_surv))
            )
            trace_with_surv[idx] <- surv_to_encode$surv_val
        }
    }

    # ========================================================================
    # 4. Subset and prepare tree
    # ========================================================================
    tree_sub <- keep.tip(tree, keep_isolates)

    tip_label_map <- setNames(
        paste0(
            isolate_patient_map[tree_sub$tip.label],
            " (",
            tree_sub$tip.label,
            ")"
        ),
        tree_sub$tip.label
    )
    tree_sub$tip.label <- tip_label_map

    # ========================================================================
    # 5. Dynamic scaling calculations
    # ========================================================================
    n_tips <- length(tree_sub$tip.label)
    n_cols <- ncol(trace_with_surv)

    font_size <- max(0.5, min(2, 20 / sqrt(n_tips)))

    # Create tree plot and extract depth for scaling
    tree_plot <- ggtree(tree_sub)
    tip_data <- tree_plot$data[tree_plot$data$isTip, ]
    tree_depth <- max(tip_data$x)

    if (is.na(tree_depth) || tree_depth <= 0) {
        tree_depth <- n_tips
    }

    heatmap_width <- max(
        tree_depth * 0.5,
        min(tree_depth * 4, n_cols * tree_depth * 0.02)
    )

    # ========================================================================
    # 6. Prepare heatmap data and color mapping
    # ========================================================================
    custom_breaks <- c(0, 1, 1.25, 1.5)
    color_map <- c(
        "0" = "white",
        "1" = presence_color,
        "1.25" = unname(surv_colors["neg"]),
        "1.5" = unname(surv_colors["pos"])
    )

    trace_df <- as.data.frame(trace_with_surv)
    trace_df[] <- lapply(trace_df, function(x) {
        factor(as.numeric(as.character(x)), levels = custom_breaks)
    })

    row.names(trace_df) <- tip_label_map[row.names(trace_df)]
    trace_df <- trace_df[tree_sub$tip.label, , drop = FALSE]

    col_lab <- colnames(trace_df)
    label_interval <- max(1, floor(n_cols / 20))
    idx_keep_labels <- seq(1, length(col_lab), label_interval)
    col_lab[setdiff(seq_along(col_lab), idx_keep_labels)] <- ""

    # ========================================================================
    # 7. Build base plot
    # ========================================================================
    tree_plot <- tree_plot +
        geom_tiplab(size = font_size, align = TRUE, offset = 0.5, linetype = NULL)

    base_plot <- gheatmap(
        tree_plot,
        trace_df,
        offset = 3.5,
        width = heatmap_width / tree_depth,
        color = NULL,
        font.size = font_size * 0.8,
        custom_column_labels = col_lab,
        colnames_angle = 90,
        colnames_offset_y = -n_tips * 0.02
    ) +
        vexpand(0.1, -1) +
        scale_fill_manual(
            name = "Trace",
            values = color_map,
            breaks = custom_breaks,
            labels = c("absent", "present", "negative surv.", "positive surv."),
            na.value = "white",
            drop = FALSE
        )

    # ========================================================================
    # 8. Add annotations (cluster)
    # ========================================================================
    # Only include clusters that are actually present in the filtered data
    clusters_present <- unique(isolate_lookup$cluster[
        isolate_lookup$isolate_id %in% keep_isolates
    ])
    cluster_order_filtered <- cluster_order[cluster_order %in% clusters_present]

    annotation_row <- data.frame(
        id = factor(tip_label_map[keep_isolates], levels = tree_sub$tip.label),
        Cluster = factor(
            as.character(isolate_lookup$cluster[match(
                keep_isolates,
                isolate_lookup$isolate_id
            )]),
            levels = as.character(cluster_order_filtered)
        )
    )

    cluster_colors <- setNames(
        iwanthue(length(cluster_order_filtered)),
        as.character(cluster_order_filtered)
    )

    base_plot <- base_plot +
        new_scale_fill() +
        geom_fruit(
            data = annotation_row,
            geom = geom_tile,
            mapping = aes(y = id, x = 1, fill = Cluster),
            offset = -0.2,
            width = 0.3
        ) +
        scale_fill_manual(
            name = "Cluster",
            values = cluster_colors,
            limits = as.character(cluster_order_filtered),
            drop = FALSE
        )

    # ========================================================================
    # 9. Final theme adjustments
    # ========================================================================
    base_plot +
        theme(
            legend.title = element_text(size = 11),
            legend.text = element_text(size = 9),
            legend.position = "right",
            legend.key.height = unit(0.8, "cm"),
            legend.key.width = unit(0.5, "cm")
        )
}
