#' Dot plot ranking LR pairs by number of significant hotspot bins
#'
#' @param LR_results Data frame returned in the \code{LR_results} slot of a BLISA
#'   result list. Rows are ligand-receptor pairs; must contain columns
#'   \code{sig_numbers} and \code{annotation}.
#' @param top Integer or \code{NULL}. Number of top LR pairs (by
#'   \code{sig_numbers}) to display. Default \code{30}.
#' @param pt.size Numeric. Point size passed to \code{geom_point}. Default 4.
#' @param flip Logical. When \code{TRUE}, LR pairs are placed on the x-axis
#'   and the hotspot count on the y-axis (vertical orientation). Default
#'   \code{FALSE} (LR pairs on y-axis, horizontal orientation).
#'
#' @return A \code{ggplot} object.
#' @export
plotLRrank <- function(LR_results, top = 30, pt.size = 4, flip = FALSE) {

  if (!is.null(top)) {
    LR_results <- LR_results[seq_len(min(top, nrow(LR_results))), ]
  }

  n_shown <- nrow(LR_results)
  title   <- paste0("Top ", n_shown, " LR pairs by hotspot count")

  LR_results$LR_pair <- rownames(LR_results)
  LR_results <- LR_results[order(-LR_results$sig_numbers), ]
  # flip=FALSE: levels low→high so highest appears at top of y-axis
  # flip=TRUE:  levels high→low so highest appears at left of x-axis
  LR_results$LR_pair <- factor(LR_results$LR_pair,
                            levels = if (flip) LR_results$LR_pair else rev(LR_results$LR_pair))

  color_scale <- scale_color_manual(
    values = c(
      "Secreted Signaling"    = "#90a955",
      "ECM-Receptor"          = "#219ebc",
      "Cell-Cell Contact"     = "#f7b801",
      "Non-protein Signaling" = "#9f86c0"
    )
  )

  if (flip) {
    p <- ggplot(LR_results, aes(x = LR_pair, y = sig_numbers, color = annotation)) +
      geom_point(size = pt.size) +
      color_scale +
      scale_y_continuous(expand = expansion(add = 100)) +
      labs(x = "Ligand–Receptor Pair", y = "Sig Spot Numbers",
           color = "Annotation", title = title) +
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title    = element_text(size = 12, face = "bold"),
        legend.text     = element_text(size = 10),
        panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.text.x     = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y     = element_text(angle = 0, hjust = 1)
      ) +
      coord_cartesian(clip = "off")
  } else {
    p <- ggplot(LR_results, aes(x = sig_numbers, y = LR_pair, color = annotation)) +
      geom_point(size = pt.size) +
      color_scale +
      scale_x_continuous(expand = expansion(add = 100)) +
      labs(x = "Sig Spot Numbers", y = "Ligand–Receptor Pair",
           color = "Annotation", title = title) +
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title    = element_text(size = 12, face = "bold"),
        legend.text     = element_text(size = 10),
        panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.text.y     = element_text(size = 12),
        axis.text.x     = element_text(angle = 45, hjust = 1)
      ) +
      coord_cartesian(clip = "off")
  }

  return(p)
}
