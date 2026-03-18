#' Plot overlap permutation test results as violin plots
#'
#' @description
#' Creates violin plots comparing observed overlap fractions against permuted
#' distributions for facility, floor, and room levels, for both concurrent and
#' sequential overlap.
#'
#' @param perm_result A list returned by [transclust::cluster_overlap_perm_test()], containing
#'   `observed`, `permuted`, and `valid_clusters`.
#'
#' @return A ggplot object with violin plots for each overlap type.
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_point position_jitter
#'   scale_color_manual labs theme_minimal theme element_text facet_wrap
#' @importFrom stats median
#' @export
plot_overlap_perm_test <- function(perm_result) {
    observed <- perm_result$observed
    perm_array <- perm_result$permuted
    valid_clusters <- perm_result$valid_clusters

    trace_types <- dimnames(perm_array)[[2]]

    # Build a long data frame from permuted results
    perm_rows <- lapply(trace_types, function(tt) {
        vals <- as.vector(perm_array[, tt, ])
        data.frame(
            trace_type = tt,
            value = vals[!is.na(vals)],
            source = "Permuted"
        )
    })
    perm_df <- do.call(rbind, perm_rows)

    # Build a long data frame from observed results
    obs_rows <- lapply(trace_types, function(tt) {
        vals <- observed[[tt]]
        vals <- vals[!is.na(vals)]
        data.frame(
            trace_type = tt,
            value = vals,
            source = "Observed"
        )
    })
    obs_df <- do.call(rbind, obs_rows)

    # Label mapping for display
    type_labels <- c(
        facility = "Facility",
        floor = "Floor",
        room = "Room",
        seq_facility = "Sequential Facility",
        seq_floor = "Sequential Floor",
        seq_room = "Sequential Room"
    )
    perm_df$trace_label <- factor(
        type_labels[perm_df$trace_type],
        levels = type_labels
    )
    obs_df$trace_label <- factor(
        type_labels[obs_df$trace_type],
        levels = type_labels
    )

    ggplot(perm_df, aes(x = trace_label, y = value)) +
        geom_violin(fill = "lightblue", alpha = 0.6, scale = "width") +
        geom_point(
            data = obs_df,
            aes(x = trace_label, y = value),
            color = "red",
            size = 2,
            alpha = 0.7,
            position = position_jitter(width = 0.05)
        ) +
        labs(
            title = "Overlap Permutation Test Results",
            x = "Overlap Type",
            y = "Fraction of Converts with Overlap"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)
        )
}
