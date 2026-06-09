#' Plot overlap permutation test results as violin plots
#'
#' @description
#' Renders the core convert-weighted overlap figure: for each trace type a violin
#' of the permuted (null) pooled overlap fractions with a boxplot inside, the
#' observed pooled fraction as a red diamond, and the observed counts
#' (`n_overlap`/`n_converts`) labelled above each violin. Violins are filled by
#' overlap group (co-occurrence vs sequential).
#'
#' This function draws only the core figure. Aggregation across clusters /
#' sequence types, significance/p-value annotations, and any group dividers or
#' headers that depend on the specific trace-type layout are left to the caller
#' to compute and add to the returned plot.
#'
#' @param perm_df A long data frame of the permuted (null) distribution with
#'   columns `trace_type` (a factor whose levels define the x-axis order) and
#'   `overlap_fraction`. One row per (trace type, permutation).
#' @param obs_df A data frame of observed values with one row per trace type:
#'   columns `trace_type` (a factor with the same levels as `perm_df`),
#'   `overlap_fraction`, `n_overlap` and `n_converts`.
#' @param title Plot title.
#' @param subtitle Plot subtitle. Defaults to `NULL` (no subtitle).
#' @param trace_labels A named character vector mapping `trace_type` levels to
#'   x-axis labels. The default maps facility/floor/room (and their `seq_`
#'   counterparts) to "Facility"/"Floor"/"Room".
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

    # Co-occurrence (facility/floor/room) vs sequential (seq_*) grouping
    classify_group <- function(tt) {
        ifelse(grepl("^seq_", as.character(tt)), "Sequential", "Co-occurrence")
    }

    perm_df$group <- factor(
        classify_group(perm_df$trace_type),
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
