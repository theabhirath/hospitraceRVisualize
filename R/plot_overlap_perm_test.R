#' Plot overlap permutation test results as violin plots
#'
#' @description
#' Renders the core convert-weighted overlap figure: per trace type, a violin of
#' the permuted (null) pooled overlap fractions with an inner boxplot, the
#' observed pooled fraction as a red diamond, and the observed counts
#' (`n_overlap`/`n_converts`) labeled above each violin. Violins are filled by
#' overlap group (co-occurrence vs sequential).
#'
#' Only the core figure is drawn. Cross-cluster/sequence-type aggregation,
#' significance annotations, and trace-type-specific dividers or headers are left
#' to the caller to add to the returned plot.
#'
#' @param perm_df Long null-distribution data frame, one row per (trace type,
#'   permutation), with columns `trace_type` (factor; levels set x-axis order)
#'   and `overlap_fraction`.
#' @param obs_df Observed values, one row per trace type, with columns
#'   `trace_type` (same levels as `perm_df`), `overlap_fraction`, `n_overlap`,
#'   `n_converts`.
#' @param title Plot title.
#' @param subtitle Plot subtitle, or `NULL` for none.
#' @param trace_labels Named character vector mapping `trace_type` levels to
#'   x-axis labels.
#' @param x_label,y_label Axis titles.
#'
#' @return A ggplot object with violin plots for each overlap type.
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_point geom_text
#' @importFrom ggplot2 scale_x_discrete scale_fill_manual scale_y_continuous
#' @importFrom ggplot2 coord_cartesian labs theme_minimal theme element_text margin
#' @importFrom rlang .data
#' @export
plot_overlap_perm_test <- function(
    perm_df,
    obs_df,
    title = "Overlap Permutation Test Results",
    subtitle = NULL,
    trace_labels = c(
        facility = "Facility",
        floor = "Floor",
        room = "Room",
        seq_facility = "Facility",
        seq_floor = "Floor",
        seq_room = "Room"
    ),
    x_label = "Trace type",
    y_label = "Fraction of converts with overlap (convert-weighted)"
) {
    trace_types <- levels(perm_df$trace_type)

    # seq_* trace types are sequential overlap; the rest are co-occurrence
    perm_df$group <- factor(
        ifelse(grepl("^seq_", as.character(perm_df$trace_type)), "Sequential", "Co-occurrence"),
        levels = c("Co-occurrence", "Sequential")
    )

    # Observed-count label ("n_overlap/n_converts (fraction)") above each violin
    obs_idx <- match(trace_types, as.character(obs_df$trace_type))
    n_overlap <- obs_df$n_overlap[obs_idx]
    n_converts <- obs_df$n_converts[obs_idx]
    obs_fraction <- obs_df$overlap_fraction[obs_idx]
    count_label <- ifelse(
        n_converts > 0,
        sprintf("%.0f/%.0f (%.2f)", n_overlap, n_converts, obs_fraction),
        sprintf("%.0f/%.0f", n_overlap, n_converts)
    )
    annot_df <- data.frame(
        trace_type = factor(trace_types, levels = trace_types),
        count_label = count_label,
        count_y = 1.0
    )

    group_fills <- c("Co-occurrence" = "#A6CEE3", "Sequential" = "#B2DF8A")

    ggplot(perm_df, aes(x = .data$trace_type, y = .data$overlap_fraction)) +
        geom_violin(aes(fill = .data$group), color = "black", alpha = 0.4) +
        geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
        geom_point(
            data = obs_df,
            aes(x = .data$trace_type, y = .data$overlap_fraction),
            shape = 23,
            size = 5,
            fill = "red",
            color = "black"
        ) +
        geom_text(
            data = annot_df,
            aes(x = .data$trace_type, y = .data$count_y, label = .data$count_label),
            size = 3.2,
            fontface = "bold",
            vjust = 0
        ) +
        scale_x_discrete(labels = trace_labels) +
        scale_fill_manual(values = group_fills, name = NULL, guide = "none") +
        labs(
            title = title,
            subtitle = subtitle,
            x = x_label,
            y = y_label
        ) +
        scale_y_continuous(breaks = seq(0, 1, 0.2)) +
        coord_cartesian(ylim = c(0, 1), clip = "off") +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 16, margin = margin(b = 6)),
            plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 30)),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
        )
}
