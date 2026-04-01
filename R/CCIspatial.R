library(data.table)
library(dplyr)
library(ggplot2)
library(SummarizedExperiment)

cols <- c(
  # Highly distinct / saturated first
  "#E69F00", # orange
  "#56B4E9", # light blue
  "#009E73", # green
  "#D55E00", # vermillion
  "#CC79A7", # magenta
  "#117a77", # blue
  "#F0E442", # yellow
  "#0b2b5e", # strong blue
  "#33A02C", # strong green
  "#E31A1C", # red
  "#6A3D9A", # purple
  "#B15928", # brown

  # Softer but still distinct
  "#FB8072", # coral
  "#80B1D3", # sky blue
  "#FDB462", # orange pastel
  "#B3DE69", # lime pastel
  "#CAB2D6", # lavender
  "#FB9A99", # salmon
  "#a18e6a", # gold pastel
  "#FF7F00", # bright orange

  # Pastel / similar shades pushed later
  "#8DD3C7", # teal pastel
  "#FFFFB3", # pale yellow
  "#BEBADA", # pale lavender
  "#FCCDE5", # pink pastel
  "#BC80BD", # violet pastel
  "#CCEBC5", # mint
  "#FFED6F", # soft yellow
  "#A6CEE3", # pale blue
  "#B2DF8A", # pale green
  "#FFB3BA"  # baby pink
)

#' plot spatial plot of cell-cell interaction pattern for a specific LR pair
#'
#' @param spe cell-level spe object
#' @param BLISA_output result outout from runBLISA
#' @param index LR index
#' @export
#'
CCIspatial <- function(
    spe,
    BLISA_output,
    index,
    ct_group = "cell_type",
    top = 30,
    hex_size = 50,
    dmax = 250
) {
  # 1. Setup Data
  LRI_sum <- BLISA_output$LR_out
  hex_sf  <- BLISA_output$hex_sf # contains n_cells

  interaction <- unname(unlist(LRI_sum[index, c("ligand.symbol", "receptor.symbol")]))
  sigHH <- LRI_sum$sig_index[[index]]
  mode  <- LRI_sum$ccc_mode[index]

  # 2. Neighbors logic (spdep)
  centroids <- sf::st_centroid(hex_sf)
  coords <- sf::st_coordinates(centroids)

  if (mode == "nearby") {
    nb_list <- spdep::dnearneigh(centroids, 0, 1.2 * hex_size)
  } else {
    nb_list <- spdep::dnearneigh(coords, 0, dmax)
  }

  # 3. Efficient Mapping
  cell_to_hex <- get_cell_hex_mapping(spe, hex_sf)

  cell_data <- data.table::data.table(
    hex_id = as.integer(cell_to_hex),
    ct = as.character(SummarizedExperiment::colData(spe)[[ct_group]]),
    ligand_expr = as.numeric(counts(spe)[interaction[1], ]),
    receptor_expr = as.numeric(counts(spe)[interaction[2], ])
  )

  # 4. Interaction Scoring (Receiver in HH, Sender in Neighbors)
  rcpt_summary <- cell_data[hex_id %in% sigHH, .(r_sum = sum(receptor_expr)), by = .(hex_id, ct_r = ct)]

  hh_nb_map <- data.table::data.table(
    hh_hex = rep(sigHH, sapply(nb_list[sigHH], length)),
    nb_hex = unlist(nb_list[sigHH])
  )

  lig_summary <- merge(hh_nb_map, cell_data, by.x = "nb_hex", by.y = "hex_id", allow.cartesian = TRUE)
  lig_summary <- lig_summary[, .(l_sum = sum(ligand_expr)), by = .(hh_hex, ct_l = ct)]

  # Merge and find top pair per hotspot
  merged_scores <- merge(rcpt_summary, lig_summary, by.x = "hex_id", by.y = "hh_hex", allow.cartesian = TRUE)
  merged_scores[, product := (log10(r_sum + 1) + log10(l_sum + 1)) / 2]
  merged_scores[, cell_pair := paste(ct_l, ct_r, sep = " → ")]

  top_pairs <- merged_scores[merged_scores[, .I[which.max(product)], by = hex_id]$V1]

  # 5. Legend and Filtering Logic
  tbl <- sort(table(top_pairs$cell_pair), decreasing = TRUE)

  if (length(tbl) <= top) {
    filtered_pairs <- names(tbl)
    legend_title <- "All pairs (ordered)"
  } else {
    filtered_pairs <- names(tbl[1:top])
    legend_title <- paste0("Top ", top, " pairs (ordered)")
  }

  top_pairs <- top_pairs %>% mutate(cell_pair_plot = ifelse(cell_pair %in% filtered_pairs, cell_pair, "rare pairs"))

  # 6. Categorize ALL hexagons
  # Start by identifying bins with cells vs empty bins
  hex_sf$cell_pair_plot <- ifelse(hex_sf$n_cells > 0, "Non-Significant", "Empty")

  # Overwrite hotspots with their top interacting pair
  hex_sf$cell_pair_plot[top_pairs$hex_id] <- top_pairs$cell_pair_plot

  # Order levels for the legend
  legend_levels <- c(filtered_pairs, "rare pairs", "Non-Significant", "Empty")
  hex_sf$cell_pair_plot <- factor(hex_sf$cell_pair_plot, levels = legend_levels)

  # 7. Visualization
  p <- ggplot(hex_sf) +
    geom_sf(aes(fill = cell_pair_plot), color = NA) +
    scale_fill_manual(
      values = c(
        "Empty" = "#F0F0F0",
        "Non-Significant" = "#D3D3D3",
        "rare pairs" = "#818589",
        setNames(cols[1:length(filtered_pairs)], filtered_pairs)
      )
    ) +
    theme_void() +
    guides(fill = guide_legend(title = legend_title)) +
    labs(
      title = paste0(rownames(LRI_sum)[index], ": Interacting Hotspots"),
      subtitle = "Grey scale: Light (Empty) | Medium (Non-sig cells)"
    )

  return(p)
}
