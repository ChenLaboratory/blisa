#' Dot plot ranking LR pairs by number of significant hotspot bins
#'
#' @param LR_out Data frame returned in the \code{LR_out} slot of a BLISA
#'   result list. Rows are ligand-receptor pairs; must contain columns
#'   \code{sig_numbers} and \code{annotation}.
#' @param top Integer or \code{NULL}. If set, only the top \code{top} rows
#'   (already ordered by \code{sig_numbers}) are plotted. Default \code{NULL}
#'   plots all pairs.
#' @param pt.size Numeric. Point size passed to \code{geom_point}. Default 4.
#'
#' @return A \code{ggplot} object.
#' @export
plotLRIsum <- function(LR_out, top = NULL, pt.size = 4) {

  if (!is.null(top)) {
    LR_out <- LR_out[seq_len(min(top, nrow(LR_out))), ]
  }

  LR_out$LR_pair <- rownames(LR_out)
  LR_out <- LR_out[order(-LR_out$sig_numbers), ]
  LR_out$LR_pair <- factor(LR_out$LR_pair, levels = rev(LR_out$LR_pair))

  p <- ggplot(LR_out, aes(x = `sig_numbers`, y = LR_pair, color = annotation)) +
    geom_point(size = pt.size)+
    scale_color_manual(
      values = c(
        "Secreted Signaling" = "#90a955",   # blue
        "ECM-Receptor"       = "#219ebc",   # green
        "Cell-Cell Contact"  = "#f7b801",   # red
        "Non-protein Signaling" = "#9f86c0" # purple
      )
    ) +
    scale_x_continuous(expand = expansion(add = 100)) +
    labs(x = "Sig Spot Numbers", y = "Ligand–Receptor Pair", color = "Annotation") +
    theme_minimal()+
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.text.y = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )+
    coord_cartesian(clip = "off")

  return(p)

}
