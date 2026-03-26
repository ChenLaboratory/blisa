library(ggplot2)
library(sf)


LRI_spatial_colors <- c("#FFFFCC", "#FFD700", "#FF7F00", "#D7301F")

plotLRI.sf <- function(hex_sf, LR_out, index, log_pval = TRUE) {

  # 1. Prepare Data
  interaction <- rownames(LR_out)[index]
  sig_indices <- LR_out$sig_index[[index]]
  p_values <- LR_out$sig_pval[[index]]
  sig_num <- LR_out$sig_numbers[index]

  # Initialize the status column
  # Default: "Empty" (no cells)
  hex_sf$bin_status <- "Empty"

  # Bins with cells but not significant
  hex_sf$bin_status[hex_sf$n_cells > 0] <- "Non-Significant"

  # Bins that are significant
  hex_sf$bin_status[sig_indices] <- "Significant"

  # Map the actual values for the color scale
  hex_sf$plot_val <- NA
  if (log_pval) {
    hex_sf$plot_val[sig_indices] <- -log10(p_values)
    lgd_title <- "-log10(pval)"
  } else {
    hex_sf$plot_val[sig_indices] <- 1 - p_values
    lgd_title <- "1-pval"
  }

  # 2. Plot
  p <- ggplot(hex_sf) +
    # Layer 1: The "Background" (Empty and Non-Significant bins)
    geom_sf(aes(fill = bin_status), color = NA) +

    # Layer 2: The "Significance" (Continuous gradient for hotspots)
    # We use 'new_scale_fill' if using ggnewscale,
    # but here we can just use a single fill scale with manual overrides
    scale_fill_manual(
      values = c(
        "Empty" = "#F0F0F0",          # Very light grey
        "Non-Significant" = "#D3D3D3" # Medium light grey
      ),
      guide = "none" # Hide the discrete legend
    ) +

    # Layer 3: Overlay the significant gradients
    # Note: To do two fill scales in one ggplot, we use the ggnewscale package
    # If you don't have it, we'll use a trick with color/fill or subsetting
    ggnewscale::new_scale_fill() +
    geom_sf(data = subset(hex_sf, bin_status == "Significant"),
            aes(fill = plot_val), color = NA) +
    scale_fill_gradientn(
      colours = LRI_spatial_colors,
      name = lgd_title,
      na.value = "transparent"
    ) +
    labs(
      title = paste0(interaction, " - ", sig_num, " hotspots"),
      subtitle = "Grey scale: Light (Empty bins) | Medium (Non-sig bins)"
    ) +
    theme_void()

  return(p)
}

# plotLRI <- function(spe, LRI_sum, index, STtype = c("iST", "sST"), size = 1, alpha = 0.8, img = FALSE, pol.border = FALSE, log_pval = T) {
#   STtype <- match.arg(STtype)
#
#   interaction <- rownames(LRI_sum)[index]
#   num <- LRI_sum$sig_spot_numbers[index]
#   plot_title <- paste0(interaction, " - ", num, " hotspots")
#
#   # prepare pval vector
#   lisa_pval_HH <- c(1:dim(spe)[2])
#   if (log_pval) {
#     lisa_pval_HH[LRI_sum$`sig_index`[[index]]] <- -log(LRI_sum$`sig_pval`[[index]], 10)
#     lgd_title <- "-log10(pval)"
#   } else {
#     lisa_pval_HH[LRI_sum$`sig_index`[[index]]] <- 1 - LRI_sum$`sig_pval`[[index]]
#     lgd_title <- "1-pval"
#   }
#
#   lisa_pval_HH[-LRI_sum$`sig_index`[[index]]] <- NA
#
#   if (STtype == "iST") {
#     spe_hex <- spe
#
#     spe_hex@metadata$grid_density$plotP <- lisa_pval_HH
#
#     # get grid_data
#     grid_data <- spe_hex@metadata$grid_density
#
#     # Build hex polygons with p-values
#     poly <- grid2df(
#       spe_hex,
#       grid_data$node_x,
#       grid_data$node_y,
#       group = grid_data$plotP
#     )
#
#     # Plot
#     p <- plotImage(spe_hex) +
#       geom_polygon(
#         data = poly,
#         aes(
#           x = X,
#           y = Y,
#           group = L2,
#           fill = group
#         ),
#         alpha = alpha,
#         color = if (pol.border) "black" else NA
#       ) +
#       coord_fixed() +
#       theme_minimal() +
#       labs(
#         title = plot_title,
#         x = NULL, y = NULL,
#         fill = lgd_title
#       ) +
#       scale_fill_gradientn(
#         colours = LRI_spatial_colors,
#         na.value = "grey70"
#       )+
#       theme(legend.position = "right",
#             legend.title = element_text(size = 12, face = "bold"),
#             legend.text = element_text(size = 10),
#             axis.title.x = element_blank(),
#             axis.title.y = element_blank())
#
#     # Update bounds (same as plotGrid)
#     p <- update_bound(
#       p,
#       x = spe_hex@metadata$grid_info$xlim + c(-1, 1) * spe_hex@metadata$grid_info$xstep / 2,
#       y = spe_hex@metadata$grid_info$ylim + c(-1, 1) * spe_hex@metadata$grid_info$ystep / 2
#     )
#
#   } else if (STtype == "sST") {
#     spe$plotP <- lisa_pval_HH
#
#     if(!img & !is.null(imgData(spe))) {
#       imgData(spe) <- NULL
#     }
#
#     p <- plotSpatial(spe, group.by = "plotP", pt.size = size, pt.alpha = alpha) +
#       coord_fixed() +
#       theme_minimal() +
#       labs(
#         title = plot_title,
#         x = NULL, y = NULL,
#         fill = lgd_title
#       ) +
#       scale_color_gradientn(
#         name = lgd_title,
#         colours = LRI_spatial_colors,
#         na.value = "grey70"
#       )+
#       theme(legend.position = "right",
#             legend.title = element_text(size = 12, face = "bold"),
#             legend.text = element_text(size = 10),
#             axis.title.x = element_blank(),
#             axis.title.y = element_blank())
#   }
#
#   p <- p + theme_void() # remove x/y axis
#
#   return(p)
#
# }
#
#
# plot_hex_lisa <- function(hex_sf, LR_out, index, log_pval = TRUE) {
#
#   # Get interaction info
#   interaction <- rownames(LR_out)[index]
#   sig_indices <- LR_out$sig_index[[index]]
#   p_values <- LR_out$sig_pval[[index]]
#
#   # Create a placeholder vector for all hexagons
#   plot_vals <- rep(NA, nrow(hex_sf))
#
#   # Calculate mapping values (-log10 pval or 1-pval)
#   if (log_pval) {
#     plot_vals[sig_indices] <- -log10(p_values)
#     lgd_title <- "-log10(pval)"
#   } else {
#     plot_vals[sig_indices] <- 1 - p_values
#     lgd_title <- "1-pval"
#   }
#
#   # Add the values to the sf object
#   hex_sf$plot_val <- plot_vals
#
#   # Plot using geom_sf
#   p <- ggplot(hex_sf) +
#     geom_sf(aes(fill = plot_val), color = NA) + # Set color = "black" for borders
#     scale_fill_gradientn(
#       colours = LRI_spatial_colors,
#       na.value = "grey90", # Light grey for non-significant/empty hexagons
#       name = lgd_title
#     ) +
#     labs(
#       title = paste0(interaction, " Hotspots"),
#       subtitle = paste0("Significant hexagons: ", length(sig_indices))
#     ) +
#     theme_void() +
#     theme(
#       legend.position = "right",
#       plot.title = element_text(hjust = 0.5, face = "bold"),
#       legend.title = element_text(size = 10, face = "bold")
#     )
#
#   return(p)
# }
#
#
# plot_hex_lisa_base <- function(hex_sf, LR_out, index, log_pval = TRUE) {
#
#   # 1. Extract Data
#   interaction <- rownames(LR_out)[index]
#   sig_indices <- LR_out$sig_index[[index]]
#   p_values <- LR_out$sig_pval[[index]]
#
#   # 2. Create the combined value column
#   # Default everything to -2 (Empty / Very Light Grey)
#   hex_sf$plot_val <- -2
#
#   # Set bins with cells to -1 (Non-sig / Medium Light Grey)
#   hex_sf$plot_val[hex_sf$n_cells > 0] <- -1
#
#   # Set significant bins to their actual positive values
#   if (log_pval) {
#     hex_sf$plot_val[sig_indices] <- -log10(p_values)
#     lgd_label <- "-log10(pval)"
#   } else {
#     hex_sf$plot_val[sig_indices] <- 1 - p_values
#     lgd_label <- "1-pval"
#   }
#
#   # 3. Define the custom gradient logic
#   # We need to calculate where the "0" point is in the scale relative to the min/max
#   val_min <- -2
#   val_max <- max(hex_sf$plot_val, na.rm = TRUE)
#
#   # This creates a normalized position for the colors
#   # 0 is Empty, slightly above 0 is Non-sig, then the rest is the palette
#   rescale_mid <- (0 - val_min) / (val_max - val_min)
#   rescale_low <- (-1 - val_min) / (val_max - val_min)
#
#   # 4. Plot
#   p <- ggplot(hex_sf) +
#     geom_sf(aes(fill = plot_val), color = NA) +
#     scale_fill_gradientn(
#       colors = c("#F0F0F0", "#D3D3D3", "#FFFFCC", "#FFD700", "#FF7F00", "#D7301F"),
#       values = c(0, rescale_low, rescale_mid,
#                  rescale_mid + 0.01, (rescale_mid + 1)/2, 1),
#       limits = c(-2, val_max),
#       name = lgd_label,
#       # Hide the negative numbers in the legend
#       breaks = seq(0, floor(val_max)),
#       labels = seq(0, floor(val_max))
#     ) +
#     labs(title = interaction) +
#     theme_void()
#
#   return(p)
# }
