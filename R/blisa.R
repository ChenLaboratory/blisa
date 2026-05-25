#' Run BLISA spatial cell-cell communication analysis
#'
#' Generic function for running BLISA (Bivariate Local Indicator of Spatial
#' Association). Dispatches on the class of \code{x}:
#' \itemize{
#'   \item \code{blisa.default} accepts a pre-binned gene-by-bin count
#'     matrix and a matching \code{bins} polygon object.
#'   \item \code{blisa.SpatialExperiment} accepts a cell-level
#'     \code{SpatialExperiment} object and bins cells into hexagonal tiles
#'     internally via \code{\link{hexBinCells}} before running the analysis.
#' }
#'
#' @param x A gene-by-bin count matrix (for \code{blisa.default}) or a
#'   cell-level \code{SpatialExperiment} object (for
#'   \code{blisa.SpatialExperiment}).
#' @param ... Additional arguments passed to the relevant method.
#'
#' @return A list; see individual method documentation for details.
#' @export
blisa <- function(x, ...) UseMethod("blisa")


#' @describeIn blisa Method for a gene-by-bin count matrix.
#'
#' @param bins An \code{sf} object of bin polygons. Row order must match the
#'   column order of \code{x}.
#' @param LR_df Data frame of ligand-receptor pairs with columns
#'   \code{ligand.symbol} and \code{receptor.symbol}. When \code{NULL},
#'   CellChatDB for the chosen \code{species} is downloaded automatically.
#' @param hex_size Numeric. Bin spacing used to define queen adjacency for
#'   nearby-mode interactions. Default \code{50}.
#' @param dmax Numeric. Maximum distance for diffuse-mode neighbours.
#'   Default \code{250}.
#' @param nsim Integer. Number of permutations for Moran's I significance.
#'   Default \code{999}.
#' @param p_cutoff Numeric. P-value threshold for High-High hotspots.
#'   Default \code{0.05}.
#' @param min_ligand Numeric. Minimum ligand count threshold. Default \code{10}.
#' @param min_receptor Numeric. Minimum receptor count threshold.
#'   Default \code{10}.
#' @param min_cells_per_bin Integer. Bins with fewer cells are excluded from
#'   Moran's I and assigned neutral statistics (\emph{p} = 1, LISA = 0).
#'   Ignored when \code{n_cells_col = NA}. Default \code{1}.
#' @param n_cells_col Character or \code{NA}. Column in \code{bins} holding
#'   per-bin cell counts used for \code{min_cells_per_bin} filtering.
#'   Set to \code{NA} to skip (default).
#' @param col Character. Column in \code{LR_df} specifying interaction
#'   category used for communication-mode assignment. Default
#'   \code{"annotation"}.
#' @param default_mode Character. CCC mode assigned to LR pairs whose
#'   annotation does not match \code{diffuse_category}. Default
#'   \code{"diffuse"}.
#' @param diffuse_category Character vector of annotation categories treated
#'   as diffuse signalling.
#' @param species Character. Which CellChatDB to download when
#'   \code{LR_df = NULL}. One of \code{"human"} (default) or \code{"mouse"}.
#' @param genes Character vector of gene names to consider when matching
#'   ligand-receptor pairs. Only LR pairs whose ligand and receptor genes are
#'   all present in \code{genes} and pass the expression thresholds will be
#'   tested. Defaults to \code{rownames(x)} (all genes in the count matrix).
#'   Genes specified here but absent from \code{rownames(x)} are silently
#'   ignored.
#' @param counts_by_group Named list of gene-by-bin count matrices, one per
#'   group level (e.g. cell type), as returned by \code{\link{hexBinCells}}
#'   when \code{group} is supplied. When provided, \code{\link{runCCI}} is
#'   called automatically after the BLISA loop and its output is included in
#'   the result as \code{CCI_scores}. Default \code{NULL}.
#'
#' @return An object of class \code{blisa} with four components:
#' \describe{
#'   \item{LR_results}{Data frame of BLISA results for each LR pair, including
#'     \code{ccc_mode}, \code{sig_numbers}, \code{sig_index}, \code{sig_pval},
#'     \code{all_pval}, \code{all_lisa}, \code{all_quadrant}, and original
#'     columns from \code{LR_df}. \code{all_quadrant} is a character vector
#'     of \code{spdep::hotspot} quadrant labels (\code{"High-High"},
#'     \code{"Low-Low"}, etc.) for every bin; non-tested bins are \code{NA}.}
#'   \item{bins}{Bin-level \code{sf} object of hexagonal polygons.}
#'   \item{spatial_weights}{Spatial weights list from \code{\link{computeSpatialWeights}}.}
#'   \item{CCI_scores}{\code{NULL} unless \code{counts_by_group} is supplied,
#'     in which case a wide data frame of interaction scores from
#'     \code{\link{runCCI}}: rows are \code{"Sender->Receiver"} group pairs,
#'     columns are LR pairs.}
#' }
#' @export
blisa.default <- function(
    x,
    bins,
    LR_df             = NULL,
    hex_size          = 50,
    dmax              = 250,
    nsim              = 999,
    p_cutoff          = 0.05,
    min_ligand        = 10,
    min_receptor      = 10,
    min_cells_per_bin = 1,
    n_cells_col       = NA,
    col               = "annotation",
    default_mode      = "diffuse",
    diffuse_category  = c("Secreted Signaling", "Non-protein Signaling"),
    species           = c("human", "mouse"),
    genes             = NULL,
    counts_by_group   = NULL,
    ...
) {
  sw             <- computeSpatialWeights(bins, hex_size, dmax, min_cells_per_bin, n_cells_col)
  queen_wt       <- sw$queen_wt
  dist_wt        <- sw$dist_wt
  keep_idx_queen <- sw$keep_idx_queen
  keep_idx_dist  <- sw$keep_idx_dist

  # Restrict gene panel for LR pair matching; default to all genes in x
  genes <- if (is.null(genes)) rownames(x) else intersect(genes, rownames(x))

  LR_df_filtered <- filterLRpairs(
    counts       = x[genes, , drop = FALSE],
    min_ligand   = min_ligand,
    min_receptor = min_receptor,
    LR_df        = LR_df,
    species      = species
  )

  LR_results <- LR_df_add_mode(LR_df_filtered, col, default_mode, diffuse_category)

  LR_results$sig_numbers   <- integer(nrow(LR_results))
  LR_results$sig_index     <- vector("list", nrow(LR_results))
  LR_results$sig_pval      <- vector("list", nrow(LR_results))
  LR_results$all_pval      <- vector("list", nrow(LR_results))
  LR_results$all_lisa      <- vector("list", nrow(LR_results))
  LR_results$all_quadrant  <- vector("list", nrow(LR_results))

  for (i in seq_len(nrow(LR_results))) {
    message(rownames(LR_results)[i])

    ligand   <- LR_results$ligand.symbol[i]
    receptor <- LR_results$receptor.symbol[i]
    mode     <- LR_results$ccc_mode[i]

    if (mode == "nearby") {
      wt       <- queen_wt
      keep_idx <- keep_idx_queen
    } else {
      wt       <- dist_wt
      keep_idx <- keep_idx_dist
    }
    message("ccc mode is ", mode)

    x_full <- get_min_expr(receptor, x)
    y_full <- get_min_expr(ligand, x)

    res_bv <- spdep::localmoran_bv(x_full[keep_idx], y_full[keep_idx], wt, nsim = nsim)

    full_pval <- rep(1, ncol(x))
    full_lisa <- rep(0, ncol(x))
    full_pval[keep_idx] <- res_bv[, "Pr(folded) Sim"]
    full_lisa[keep_idx] <- res_bv[, "Ibvi"]

    hs <- spdep::hotspot(
      res_bv,
      Prname        = "Pr(folded) Sim",
      cutoff        = p_cutoff,
      quadrant.type = "pysal",
      p.adjust      = "none"
    )

    hs_idx_hh <- !is.na(hs) & hs == "High-High"
    HH_idx    <- keep_idx[which(hs_idx_hh)]
    HH_pval   <- full_pval[HH_idx]

    full_quadrant <- rep(NA_character_, ncol(x))
    full_quadrant[keep_idx] <- as.character(hs)

    LR_results$sig_numbers[i]   <- length(HH_idx)
    LR_results$sig_index[[i]]   <- HH_idx
    LR_results$sig_pval[[i]]    <- HH_pval
    LR_results$all_pval[[i]]    <- full_pval
    LR_results$all_lisa[[i]]    <- full_lisa
    LR_results$all_quadrant[[i]] <- full_quadrant
    LR_results$ccc_mode[i]      <- mode
  }

  LR_results <- LR_results[order(-LR_results$sig_numbers), , drop = FALSE]
  front_cols <- c("ccc_mode", "sig_numbers", "sig_index", "sig_pval")
  LR_results <- LR_results[, c(front_cols, setdiff(colnames(LR_results), front_cols))]

  CCI_scores <- NULL
  if (!is.null(counts_by_group)) {
    message("Running CCI analysis...")
    tmp <- new_blisa(LR_results, bins, sw)
    CCI_scores <- runCCI.default(tmp, counts_by_group)
  }

  new_blisa(LR_results, bins, sw, CCI_scores)
}


#' @describeIn blisa Method for a cell-level SpatialExperiment object.
#'   Bins cells into hexagonal tiles via \code{\link{hexBinCells}} then
#'   delegates to \code{blisa.default}.
#'
#' @param bin_size Numeric. Width of each hexagonal bin in coordinate units
#'   (e.g. microns). Passed to \code{\link{hexBinCells}} and also used as
#'   \code{hex_size} for queen-adjacency computation. Default \code{50}.
#' @param group Character. Column name in \code{colData(x)} to use as the
#'   grouping variable (e.g. cell type) for per-group bin aggregation and
#'   downstream CCI analysis via \code{\link{runCCI}}. If the column is not
#'   found in \code{colData(x)}, a message is issued and CCI is skipped.
#'   Default \code{"cell_type"}.
#' @param genes Character vector of gene names to consider when matching
#'   ligand-receptor pairs. Defaults to \code{rownames(x)} (all genes in the
#'   \code{SpatialExperiment} object).
#'
#' @export
blisa.SpatialExperiment <- function(x, bin_size = 50, LR_df = NULL,
                                       group = "cell_type",
                                       genes = NULL, ...) {
  coords <- as.data.frame(SpatialExperiment::spatialCoords(x))

  # Resolve group vector from colData
  cd_cols <- colnames(SummarizedExperiment::colData(x))
  if (!is.null(group) && group %in% cd_cols) {
    group_vec <- SummarizedExperiment::colData(x)[[group]]
  } else {
    if (!is.null(group))
      message("Column '", group, "' not found in colData(x) — CCI analysis will be skipped.")
    group_vec <- NULL
  }

  if (is.null(genes)) genes <- rownames(x)

  binned <- hexBinCells(coords, SummarizedExperiment::assay(x, "counts"), bin_size = bin_size, group = group_vec)

  blisa.default(
    x               = binned$counts_matrix,
    bins          = binned$bins,
    LR_df           = LR_df,
    hex_size        = bin_size,
    n_cells_col     = "n_cells",
    genes           = genes,
    counts_by_group = binned$counts_by_group,
    ...
  )
}
