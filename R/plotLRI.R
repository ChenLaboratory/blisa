#' Spatial significance map for one ligand-receptor pair
#'
#' Plots each bin coloured by significance status: empty, non-significant,
#' or significant (continuous gradient of -log10 p-value or 1-p-value).
#'
#' @param bins The \code{bins} slot of a \code{blisa} object (i.e.
#'   \code{res$bins}).
#' @param LR_results The \code{LR_results} slot of a \code{blisa} object (i.e.
#'   \code{res$LR_results}).
#' @param index Integer. Row index into \code{LR_results} selecting the
#'   ligand-receptor pair to visualise.
#' @param log_pval Logical. If \code{TRUE} (default), colour significant bins
#'   by -log10(p-value). If \code{FALSE}, use 1 - p-value.
#'
#' @return A \code{ggplot} object.
#' @export
plotLRI.sf <- function(bins, LR_results, index, log_pval = TRUE) {

  interaction <- rownames(LR_results)[index]
  sig_indices <- LR_results$sig_index[[index]]
  p_values    <- LR_results$sig_pval[[index]]
  sig_num     <- LR_results$sig_numbers[index]

  if (log_pval) {
    plot_vals <- -log10(p_values)
    lgd_title <- "-log10(pval)"
  } else {
    plot_vals <- 1 - p_values
    lgd_title <- "1-pval"
  }

  # Pre-compute a fill colour for every bin
  fill_cols <- rep("#F0F0F0", nrow(bins))   # empty
  fill_cols[bins$n_cells > 0] <- "#D3D3D3" # non-significant

  if (length(sig_indices) > 0) {
    pval_range <- range(plot_vals)
    if (pval_range[1] == pval_range[2]) pval_range[2] <- pval_range[1] + 1
    pval_norm  <- (plot_vals - pval_range[1]) / diff(pval_range)
    pal        <- colorRampPalette(LRI_spatial_colors)(100)
    fill_cols[sig_indices] <- pal[as.integer(pval_norm * 99) + 1]
  }

  bins$fill_col <- fill_cols

  # Invisible dummy points to drive the continuous colorbar legend
  val_range <- if (length(plot_vals) > 0) range(plot_vals) else c(0, 1)
  bbox <- sf::st_bbox(bins)
  cx   <- unname((bbox["xmin"] + bbox["xmax"]) / 2)
  cy   <- unname((bbox["ymin"] + bbox["ymax"]) / 2)
  dummy <- data.frame(
    x   = cx,
    y   = cy,
    val = seq(val_range[1], val_range[2], length.out = 100)
  )

  p <- ggplot(bins) +
    geom_sf(aes(fill = fill_col), color = NA) +
    scale_fill_identity() +
    geom_point(data = dummy, aes(x = x, y = y, color = val), alpha = 0) +
    scale_color_gradientn(colours = LRI_spatial_colors, name = lgd_title) +
    guides(color = guide_colorbar(title.position = "top")) +
    labs(
      title    = paste0(interaction, " - ", sig_num, " hotspots"),
      subtitle = "Grey scale: Light (Empty bins) | Medium (Non-sig bins)"
    ) +
    theme_void() +
    theme(legend.position = "right")

  return(p)
}

