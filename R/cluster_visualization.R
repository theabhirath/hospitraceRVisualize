#' Plots a phylogeny with trace heatmap using ggtree with consistent sizing
#'
#' Creates a combined visualization of a phylogenetic tree with a patient location
#' trace heatmap. Uses ggtree for the phylogeny with continuous line segments for
#' the trace data (not discrete cells). Surveillance data is overlaid as dots.
#'
#' @param tree A phylogenetic tree object of class `phylo`.
#' @param isolate_lookup A data frame from `get_isolate_lookup()` containing columns:
#'   isolate_id, patient_id, date, cluster, adm_pos, prev_surv.
#' @param trace_data A matrix with patient IDs as row names and dates (numeric) as
#'   column names. Values are numeric categories (0=absent, 1+=location categories).
#' @param surv_df A data frame with surveillance data containing columns:
#'   patient_id, genome_id, surv_date, result (0/1 for negative/positive).
#' @param cluster_filter Numeric vector of cluster IDs to include, or NULL for all clusters.
#' @param trace_colors Color palette for trace location values.
#' @param surv_colors Named character vector of colors for surveillance types:
#'   neg (negative), pos_clust (positive cluster), pos_nonclust (positive non-cluster).
#' @param clust_patient_categories Named list from `cluster_patient_categorization()`.
#'   Each element is a named vector mapping patient_id to category. If NULL,
#'   all patients are labeled as "other".
#' @param label_colors Named character vector of colors for patient category labels.
#'   Categories from cluster_patient_categorization: index, multiply-colonized-index,
#'   weak-index, convert, adm-pos, adm-pos-convert, secondary-convert, other.
#' @param inches_per_row Height in inches allocated per patient row when saving.
#'   Default 0.2. Use this to control the visual thickness of rows in the output.
#' @param row_thickness Proportion of row spacing used for trace bar height. Default 0.4.
#' @param max_tree_width Maximum width for tree as proportion of time range. Default 0.08.
#'   Tree branch lengths are scaled to fit within this width.
#'
#' @return A ggplot object with the following attributes for consistent sizing:
#'   \itemize{
#'     \item \code{recommended_height}: Suggested height in inches
#'     \item \code{recommended_width}: Suggested width in inches
#'     \item \code{n_patients}: Number of patients in the plot
#'     \item \code{inches_per_row}: Height allocated per patient row
#'   }
#'   Use these when saving with \code{ggsave()} to maintain consistent row heights.
#'
#' @importFrom ggtree ggtree geom_tiplab
#' @importFrom ape keep.tip
#' @importFrom ggplot2 ggplot aes geom_rect geom_point scale_fill_manual
#'   scale_color_manual scale_x_continuous scale_y_continuous theme_minimal theme
#'   element_blank element_text element_rect coord_cartesian labs annotation_custom
#' @importFrom dplyr left_join filter mutate case_when group_by arrange slice_min
#'   ungroup pull summarise inner_join anti_join slice_head bind_rows select distinct
#' @importFrom paletteer paletteer_d
#' @importFrom ggnewscale new_scale_color
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @export
plot_trace_phylo_tree <- function(
    tree,
    isolate_lookup,
    trace_data,
    surv_df,
    cluster_filter = NULL,
    trace_colors = paletteer_d("ggthemes::Classic_Cyclic"),
    surv_colors = c(
        neg = "blue",
        pos_nonclust = "red",
        pos_clust = "gold"
    ),
    clust_patient_categories = NULL,
    label_colors = c(
        "index" = "red",
        "multiply-colonized-index" = "darkred",
        "weak-index" = "orange",
        "convert" = "black",
        "adm-pos" = "forestgreen",
        "adm-pos-convert" = "blue",
        "secondary-convert" = "grey",
        "other" = "purple"
    ),
    inches_per_row = 0.2,
    row_thickness = 0.4,
    max_tree_width = 0.08
) {
    # ========================================================================
    # 1. Filter data by cluster if specified
    # ========================================================================
    if (!is.null(cluster_filter)) {
        isolate_lookup <- isolate_lookup[
            isolate_lookup$cluster %in% cluster_filter,
        ]
    }
    target_clusters <- unique(isolate_lookup$cluster)

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

    tree_sub <- keep.tip(tree, keep_isolates)
    n_tips <- length(tree_sub$tip.label)
    n_cols <- ncol(trace_sub)

    # ========================================================================
    # 3. Prepare tree and calculate scaling
    # ========================================================================
    time_points <- as.numeric(colnames(trace_sub))
    time_range <- max(time_points) - min(time_points)
    target_tree_width <- time_range * max_tree_width
    time_step <- if (n_cols > 1) median(diff(time_points)) else 1

    # Check if tree has branch lengths (phylogram) or not (cladogram/parsimony)
    is_cladogram <- is.null(tree_sub$edge.length) ||
        length(tree_sub$edge.length) == 0

    # Create tip label map before modifying tree
    tip_label_map <- setNames(
        paste0(
            isolate_patient_map[tree_sub$tip.label],
            " (",
            tree_sub$tip.label,
            ")"
        ),
        tree_sub$tip.label
    )

    # Get y-coordinates from tree before changing labels (y-coords depend on topology, not branch lengths)
    # Also get tree depth for scaling calculations
    initial_tree <- ggtree(
        tree_sub,
        branch.length = if (is_cladogram) "none" else NULL
    )
    tip_data <- initial_tree$data[initial_tree$data$isTip, ]
    tip_to_y <- setNames(tip_data$y, tip_data$label)
    initial_tree_depth <- max(tip_data$x)

    # Scale branch lengths for phylograms
    if (!is_cladogram && !is.na(initial_tree_depth) && initial_tree_depth > 0) {
        tree_sub$edge.length <- tree_sub$edge.length *
            (target_tree_width / initial_tree_depth)
    }

    # Update tree labels for display and create final plot
    tree_sub$tip.label <- tip_label_map[tree_sub$tip.label]
    tree_plot <- ggtree(
        tree_sub,
        branch.length = if (is_cladogram) "none" else NULL
    )
    tree_depth <- max(tree_plot$data$x[tree_plot$data$isTip])

    # Calculate effective tree depth (with cladogram scaling if needed)
    cladogram_scale <- if (is_cladogram && tree_depth > 0) {
        target_tree_width / tree_depth
    } else {
        1
    }
    effective_tree_depth <- tree_depth * cladogram_scale

    # ========================================================================
    # 4. Calculate heatmap position
    # ========================================================================
    # Estimate label width based on longest label
    max_label_chars <- max(nchar(tip_label_map))
    estimated_label_width <- max_label_chars * time_range * 0.004

    # Heatmap starts after tree + gap + label width + buffer
    heatmap_offset <- effective_tree_depth *
        1.1 +
        estimated_label_width +
        time_range * 0.02
    time_to_x <- function(t) heatmap_offset + (t - min(time_points))

    # ========================================================================
    # 5. Create continuous segments for trace data
    # ========================================================================
    heatmap_list <- lapply(names(tip_to_y), function(iso_id) {
        y_pos <- tip_to_y[iso_id]
        values <- as.numeric(trace_sub[iso_id, ])

        # Find runs of consecutive identical values
        run_starts <- c(1, which(diff(values) != 0) + 1)
        run_ends <- c(which(diff(values) != 0), n_cols)
        n_runs <- length(run_starts)

        data.frame(
            xmin = time_to_x(time_points[run_starts] - time_step / 2),
            xmax = time_to_x(time_points[run_ends] + time_step / 2),
            y = rep(unname(y_pos), n_runs),
            value = as.character(values[run_starts])
        )
    })
    heatmap_df <- do.call(rbind, heatmap_list)

    # ========================================================================
    # 6. Prepare color mapping
    # ========================================================================
    observed_breaks <- sort(unique(as.numeric(as.matrix(trace_sub))))
    custom_breaks <- sort(unique(c(0, observed_breaks[observed_breaks > 0])))

    color_map <- setNames(
        c(
            "transparent",
            trace_colors[seq_len(min(
                length(custom_breaks) - 1,
                length(trace_colors)
            ))]
        ),
        as.character(custom_breaks[seq_len(min(
            length(custom_breaks),
            length(trace_colors) + 1
        ))])
    )

    # ========================================================================
    # 7. Categorize surveillance results and calculate positions
    # ========================================================================
    visible_patients <- unique(isolate_patient_map[keep_isolates])

    # Get all isolate_ids from the full lookup
    all_lookup_isolates <- unique(isolate_lookup$isolate_id)

    surv_filtered <- surv_df |>
        filter(.data$patient_id %in% visible_patients) |>
        mutate(
            genome_id = as.character(.data$genome_id),
            in_lookup = .data$genome_id %in% all_lookup_isolates
        ) |>
        left_join(
            isolate_lookup[, c("isolate_id", "cluster")],
            by = c("genome_id" = "isolate_id")
        ) |>
        mutate(
            is_pos_clust = .data$result == 1 &
                !is.na(.data$cluster) &
                .data$cluster %in% target_clusters
        )

    # Find patient+date combinations that have a pos_clust
    dates_with_pos_clust <- surv_filtered |>
        filter(.data$is_pos_clust) |>
        select("patient_id", "surv_date") |>
        distinct()

    # Categorize surveillances:
    # - Discard positives not in lookup ONLY if there's a pos_clust on same date
    # - Otherwise keep them as pos_nonclust
    surv_categorized <- surv_filtered |>
        mutate(
            has_pos_clust_same_date = paste(.data$patient_id, .data$surv_date) %in%
                paste(dates_with_pos_clust$patient_id, dates_with_pos_clust$surv_date),
            surv_type = case_when(
                .data$result == 0 ~ "neg",
                .data$is_pos_clust ~ "pos_clust",
                .data$result == 1 & .data$in_lookup ~ "pos_nonclust",
                .data$result == 1 & !.data$has_pos_clust_same_date ~ "pos_nonclust",
                TRUE ~ NA_character_
            )
        ) |>
        filter(!is.na(.data$surv_type))

    # ========================================================================
    # 7b. Get patient categories from clust_patient_categories parameter
    # ========================================================================
    # Extract patient categories for the filtered cluster(s)
    patient_status <- setNames(
        rep("other", length(visible_patients)),
        visible_patients
    )

    if (!is.null(clust_patient_categories) && !is.null(cluster_filter)) {
        for (cl in as.character(cluster_filter)) {
            if (cl %in% names(clust_patient_categories)) {
                cluster_cats <- clust_patient_categories[[cl]]
                for (pt in names(cluster_cats)) {
                    if (pt %in% names(patient_status)) {
                        patient_status[pt] <- cluster_cats[pt]
                    }
                }
            }
        }
    }

    # Create patient to y-coordinate lookup
    # Use keep_isolates directly since it has one isolate per patient (from slice_min)
    # and these are the exact isolates in the tree
    patient_to_y <- setNames(
        tip_to_y[keep_isolates],
        isolate_patient_map[keep_isolates]
    )

    # Calculate surveillance dot positions
    surv_plot_df <- if (nrow(surv_categorized) > 0) {
        surv_with_coords <- surv_categorized |>
            mutate(
                y = patient_to_y[as.character(.data$patient_id)],
                surv_time = as.numeric(.data$surv_date)
            ) |>
            filter(.data$surv_time %in% time_points, !is.na(.data$y))

        if (nrow(surv_with_coords) > 0) {
            data.frame(
                x = time_to_x(surv_with_coords$surv_time),
                y = surv_with_coords$y,
                surv_type = surv_with_coords$surv_type
            )
        } else {
            data.frame(x = numeric(), y = numeric(), surv_type = character())
        }
    } else {
        data.frame(x = numeric(), y = numeric(), surv_type = character())
    }

    # ========================================================================
    # 8. Build combined plot
    # ========================================================================
    font_size <- max(1.5, min(3, 30 / sqrt(n_tips)))
    dot_size <- max(0.5, min(2, 20 / sqrt(n_tips)))

    # Scale cladogram tree coordinates
    if (is_cladogram) {
        tree_plot$data$x <- tree_plot$data$x * cladogram_scale
        if ("branch" %in% names(tree_plot$data)) {
            tree_plot$data$branch <- tree_plot$data$branch * cladogram_scale
        }
    }

    # Add patient status to tree data for colored tip labels
    # Map from new tip labels back to patient_id, then to status
    tree_plot$data$patient_status <- sapply(tree_plot$data$label, function(lbl) {
        if (is.na(lbl) || lbl == "") {
            return(NA_character_)
        }
        # Extract patient_id from label format "patient_id (isolate_id)"
        pt_id <- sub(" \\(.*\\)$", "", lbl)
        if (pt_id %in% names(patient_status)) patient_status[pt_id] else "other"
    })

    # X-axis breaks at 14-day intervals
    x_axis_breaks <- seq(
        ceiling(min(time_points) / 14) * 14,
        floor(max(time_points) / 14) * 14,
        by = 14
    )

    # Dotted grid data - dots at 7-day intervals
    grid_df <- expand.grid(
        time = seq(
            ceiling(min(time_points) / 7) * 7,
            floor(max(time_points) / 7) * 7,
            by = 7
        ),
        y = seq_len(n_tips)
    )
    grid_df$x <- time_to_x(grid_df$time)

    # Build plot
    p <- tree_plot +
        geom_point(
            data = grid_df,
            aes(x = x, y = y),
            color = "gray80",
            size = 0.3,
            shape = 16,
            inherit.aes = FALSE
        ) +
        geom_rect(
            data = heatmap_df,
            aes(
                xmin = xmin,
                xmax = xmax,
                ymin = y - row_thickness / 2,
                ymax = y + row_thickness / 2,
                fill = value
            ),
            inherit.aes = FALSE
        ) +
        scale_fill_manual(
            name = "Location",
            values = color_map,
            na.value = "transparent"
        ) +
        geom_tiplab(
            aes(color = patient_status),
            size = font_size,
            offset = effective_tree_depth * 0.05,
            align = FALSE
        ) +
        scale_color_manual(
            name = "Patient Status",
            values = label_colors,
            na.value = "black"
        ) +
        scale_x_continuous(
            name = "Time (days)",
            breaks = time_to_x(x_axis_breaks),
            labels = x_axis_breaks
        ) +
        theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(size = 8),
            axis.title.x = element_text(size = 10),
            legend.position = "none",
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA)
        ) +
        coord_cartesian(clip = "off")

    # Add surveillance dots if present (use new_scale_color for separate color scale)
    if (nrow(surv_plot_df) > 0) {
        p <- p +
            new_scale_color() +
            geom_point(
                data = surv_plot_df,
                aes(x = x, y = y, color = surv_type),
                size = dot_size,
                inherit.aes = FALSE
            ) +
            scale_color_manual(
                name = "Surveillance",
                values = surv_colors,
                labels = c(
                    neg = "Negative",
                    pos_nonclust = "Positive (non-cluster)",
                    pos_clust = "Positive (cluster)"
                )
            )
    }

    # ========================================================================
    # 9. Calculate recommended dimensions and return
    # ========================================================================
    base_margin <- 0.8
    recommended_height <- n_tips * inches_per_row + base_margin
    recommended_width <- max(8, min(14, time_range / 8 + 4))

    # Attach dimensions as attributes
    attr(p, "recommended_height") <- recommended_height
    attr(p, "recommended_width") <- recommended_width
    attr(p, "n_patients") <- n_tips
    attr(p, "inches_per_row") <- inches_per_row

    p
}

#' Create standalone legend for trace phylogeny plots
#'
#' Generates a combined legend plot for trace locations, surveillance results,
#' and patient categories that can be saved separately from the main plot.
#'
#' @param trace_colors Color palette for trace location values. Should match
#'   the colors used in `plot_trace_phylo_tree()`.
#' @param trace_labels Character vector of labels for each trace location.
#'   Length should match `trace_colors`.
#' @param surv_colors Named character vector of colors for surveillance types.
#'   Default matches `plot_trace_phylo_tree()` defaults.
#' @param label_colors Named character vector of colors for patient category labels.
#'   Default matches `plot_trace_phylo_tree()` defaults.
#' @param include_trace Logical. Include trace location legend. Default TRUE.
#' @param include_surv Logical. Include surveillance legend. Default TRUE.
#' @param include_patient Logical. Include patient category legend. Default TRUE.
#' @param title_size Numeric. Font size for legend titles. Default 10.
#' @param text_size Numeric. Font size for legend text. Default 8.
#' @param key_size Numeric. Size of legend keys in cm. Default 0.5.
#' @param layout Character. Layout direction for combining legends.
#'   `"horizontal"` places legends side-by-side (default),
#'   `"vertical"` stacks them one below the other.
#'
#' @return A ggplot object containing only the legend(s), suitable for saving
#'   with `ggsave()`.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_tile scale_fill_manual
#'   scale_color_manual theme_void theme element_text unit guides guide_legend
#' @importFrom cowplot get_legend plot_grid
#' @importFrom paletteer paletteer_d
#' @importFrom stats setNames
#' @export
plot_trace_phylo_legend <- function(
    trace_colors = paletteer_d("ggthemes::Classic_Cyclic"),
    trace_labels = NULL,
    surv_colors = c(
        neg = "blue",
        pos_nonclust = "red",
        pos_clust = "gold"
    ),
    label_colors = c(
        "index" = "red",
        "multiply-colonized-index" = "darkred",
        "weak-index" = "orange",
        "convert" = "black",
        "adm-pos" = "forestgreen",
        "adm-pos-convert" = "blue",
        "secondary-convert" = "grey",
        "other" = "purple"
    ),
    include_trace = TRUE,
    include_surv = TRUE,
    include_patient = TRUE,
    title_size = 10,
    text_size = 8,
    key_size = 0.5,
    layout = c("horizontal", "vertical")
) {
    layout <- match.arg(layout)
    legend_list <- list()

    # Common theme for extracting legends
    legend_theme <- theme_void() +
        theme(
            legend.title = element_text(size = title_size, face = "bold"),
            legend.text = element_text(size = text_size),
            legend.key.size = unit(key_size, "cm")
        )

    # ========================================================================
    # 1. Trace location legend
    # ========================================================================
    if (include_trace && length(trace_colors) > 0) {
        # Use trace_labels length if provided, otherwise use trace_colors length
        if (is.null(trace_labels)) {
            n_items <- length(trace_colors)
            trace_labels <- paste("Location", seq_len(n_items))
        } else {
            n_items <- length(trace_labels)
        }

        trace_df <- data.frame(
            x = seq_len(n_items),
            y = 1,
            location = factor(trace_labels, levels = trace_labels)
        )

        trace_plot <- ggplot(trace_df, aes(x = x, y = y, fill = location)) +
            geom_tile() +
            scale_fill_manual(
                name = "Location",
                values = setNames(
                    as.character(trace_colors[seq_len(n_items)]),
                    trace_labels
                )
            ) +
            guides(fill = guide_legend(ncol = 1)) +
            legend_theme

        legend_list$trace <- get_legend(trace_plot)
    }

    # ========================================================================
    # 2. Surveillance legend
    # ========================================================================
    if (include_surv && length(surv_colors) > 0) {
        surv_labels <- c(
            neg = "Negative",
            pos_nonclust = "Positive (non-cluster)",
            pos_clust = "Positive (cluster)"
        )

        surv_df <- data.frame(
            x = seq_along(surv_colors),
            y = 1,
            surv_type = factor(names(surv_colors), levels = names(surv_colors))
        )

        surv_plot <- ggplot(surv_df, aes(x = x, y = y, color = surv_type)) +
            geom_point(size = 3) +
            scale_color_manual(
                name = "Surveillance",
                values = surv_colors,
                labels = surv_labels[names(surv_colors)]
            ) +
            guides(color = guide_legend(ncol = 1)) +
            legend_theme

        legend_list$surv <- get_legend(surv_plot)
    }

    # ========================================================================
    # 3. Patient category legend
    # ========================================================================
    if (include_patient && length(label_colors) > 0) {
        patient_labels <- c(
            "index" = "Index",
            "multiply-colonized-index" = "Multiply-colonized index",
            "weak-index" = "Weak index",
            "convert" = "Convert",
            "adm-pos" = "Admission positive",
            "adm-pos-convert" = "Adm. positive convert",
            "secondary-convert" = "Secondary convert",
            "other" = "Other"
        )

        patient_df <- data.frame(
            x = seq_along(label_colors),
            y = 1,
            category = factor(names(label_colors), levels = names(label_colors))
        )

        patient_plot <- ggplot(patient_df, aes(x = x, y = y, color = category)) +
            geom_text(label = "P", size = 4, fontface = "bold") +
            scale_color_manual(
                name = "Patient Category",
                values = label_colors,
                labels = patient_labels[names(label_colors)]
            ) +
            guides(color = guide_legend(ncol = 1)) +
            legend_theme

        legend_list$patient <- get_legend(patient_plot)
    }

    # ========================================================================
    # 4. Combine legends
    # ========================================================================
    if (length(legend_list) == 0) {
        stop("At least one legend type must be included.")
    }

    if (layout == "horizontal") {
        plot_grid(plotlist = legend_list, ncol = length(legend_list), align = "h")
    } else {
        plot_grid(plotlist = legend_list, nrow = length(legend_list), align = "v")
    }
}
