# Estimate the Visium spot center-to-center spacing from spot coordinates as the
# median nearest-neighbour distance (robust to tissue-edge spots). Because the
# physical Visium pitch is a known constant (100 um), the returned value doubles
# as a unit ruler: 100 / pitch is microns per coordinate unit.
.estimateSpotPitch <- function(xy) {
  xy  <- as.matrix(xy)
  nn  <- spdep::knearneigh(xy, k = 1)$nn[, 1]
  nnd <- sqrt(rowSums((xy - xy[nn, , drop = FALSE])^2))
  stats::median(nnd)
}

#' Build spot-level hexagonal bins from a Visium SpatialExperiment
#'
#' Turns a 10x Visium \code{SpatialExperiment} (as produced by
#' \code{scider::readVisium()}) into the pre-binned inputs expected by
#' \code{\link{blisa.default}}, treating each Visium spot as one bin. The spot
#' coordinates are rebuilt from the integer \code{array_col}/\code{array_row}
#' indices onto an exact, de-tilted hexagonal lattice with 100 um
#' (\code{spot_pitch}) spacing -- the standalone equivalent of
#' \code{scider::realignVisium()}.
#'
#' @details
#' On Visium, each spot is already a spatial bin, so no cell-to-bin aggregation
#' is needed (unlike \code{\link{hexBinCells}}). Because \code{array_col}
#' increments by 2 within a row and alternate rows shift \code{array_col} by 1,
#' the half-row stagger is already encoded in the indices; the lattice therefore
#' needs no odd/even phase correction. Only the bin centroids are used by
#' \code{\link{computeSpatialWeights}}; the hexagon polygons are for plotting
#' (they do tessellate the lattice exactly).
#'
#' When \code{array_col}/\code{array_row} are absent (e.g. a
#' \code{SpatialExperiment} not built by \code{scider::readVisium()}), the
#' function falls back to using \code{spatialCoords(spe)} directly: spots are
#' still one bin each, but the de-tilt/exact-lattice step is skipped. In that
#' case, if \code{spot_pitch} is \code{NULL}, the pitch is measured from the
#' coordinates as the median nearest-neighbour distance. Since the physical
#' Visium pitch is 100 um, \code{100 / pitch} is microns-per-unit; the
#' coordinates are \emph{rescaled to microns} by that factor (near-identity when
#' they are already microns), so the returned \code{pitch} is 100 and downstream
#' \code{bin_size}/\code{dmax} are in microns. If you instead supply
#' \code{spot_pitch}, no rescaling is done and everything stays in the
#' coordinates' native units.
#'
#' @param spe A Visium \code{SpatialExperiment} with a \code{"counts"} assay
#'   whose row names are gene symbols. \code{array_col}/\code{array_row} in
#'   \code{colData} are used to rebuild the exact lattice when present;
#'   otherwise \code{spatialCoords(spe)} is used as-is.
#' @param spot_pitch Numeric or \code{NULL}. Spot center-to-center spacing used
#'   as \code{bin_size}. \code{NULL} (default) means: use 100 um for the
#'   array-index lattice, or measure it from the coordinates in the fallback.
#'   Supply a number to force it (interpreted in coordinate units).
#' @param in_tissue_only Logical. Keep only spots with \code{in_tissue == 1}.
#'   Default \code{TRUE}.
#' @param hexagons Logical. If \code{TRUE} (default), build hexagonal bin
#'   polygons; if \code{FALSE}, use point geometry (faster, but
#'   \code{plotHotspots} then draws points).
#' @param min_total_counts Numeric. Drop spots whose total counts fall below
#'   this threshold. Default \code{10}. Set to \code{0} to keep all spots.
#' @param verbose Logical. Print progress messages. Default \code{FALSE}.
#'
#' @return A list with:
#' \describe{
#'   \item{counts_matrix}{Gene-by-spot count matrix (columns match \code{bins}).}
#'   \item{bins}{An \code{sf} object of per-spot bins with \code{bin_id},
#'     \code{n_cells} (always 1), \code{total_counts}, and \code{img_x}/
#'     \code{img_y} (the raw image-registered spot coordinates, for overlaying
#'     plots on the H&E image) columns. Row order matches the columns of
#'     \code{counts_matrix}.}
#'   \item{pitch}{The spot spacing used (= the recommended \code{bin_size}).}
#'   \item{coord_unit}{\code{"micron"}, \code{"pixel/other"}, or \code{NA}
#'     (when \code{spot_pitch} was supplied) -- the inferred coordinate unit.}
#' }
#' @seealso \code{\link{blisa}}, \code{\link{hexBinCells}}
#' @export
visiumSpotBins <- function(spe,
                           spot_pitch       = NULL,
                           in_tissue_only   = TRUE,
                           hexagons         = TRUE,
                           min_total_counts = 10,
                           verbose          = FALSE) {
  cd        <- SummarizedExperiment::colData(spe)
  has_array <- all(c("array_col", "array_row") %in% colnames(cd))

  # Restrict to in-tissue spots
  if (in_tissue_only) {
    if (!"in_tissue" %in% colnames(cd)) {
      if (verbose) message("No 'in_tissue' column found; using all spots.")
    } else {
      spe <- spe[, cd$in_tissue == 1]
      cd  <- SummarizedExperiment::colData(spe)
    }
  }

  counts <- SummarizedExperiment::assay(spe, "counts")

  # Keep the raw, image-registered spot coordinates so plots can overlay the
  # H&E image. The analysis lattice below may be de-tilted/rescaled and thus not
  # match the frame the image is aligned to; these do.
  img_xy <- SpatialExperiment::spatialCoords(spe)

  if (has_array) {
    # De-tilted, exactly-spaced hex lattice rebuilt from the integer array
    # indices (standalone equivalent of scider::realignVisium): array_col
    # carries the half-row stagger, so no odd/even correction is needed. The
    # rebuilt lattice is in microns, so the pitch is 100 um by construction.
    eff_pitch <- if (is.null(spot_pitch)) 100 else spot_pitch
    coord_unit <- "micron"
    x_ideal <- as.numeric(cd$array_col) * (eff_pitch / 2)
    y_ideal <- as.numeric(cd$array_row) * (eff_pitch / 2 * sqrt(3))

    # Anchor the lattice onto the original coordinates so the bins still overlay
    # the H&E image. Only safe when coordinates are in microns.
    if (identical(spe@metadata$coord_unit, "micron")) {
      xy0     <- SpatialExperiment::spatialCoords(spe)[1, ]
      x_ideal <- x_ideal + (xy0[1] - x_ideal[1])
      y_ideal <- y_ideal + (xy0[2] - y_ideal[1])
    } else if (verbose) {
      message("coord_unit is not 'micron'; lattice origin left at (0,0) ",
              "(bins will not overlay the image).")
    }
  } else {
    # No array indices (e.g. an SPE not built by scider::readVisium): fall back
    # to the stored spot coordinates directly. Spots are still one-bin-each; we
    # only skip the de-tilt/exact-lattice step, which blisa's adjacency radius
    # (1.2 * bin_size) tolerates anyway.
    xy <- SpatialExperiment::spatialCoords(spe)

    if (is.null(spot_pitch)) {
      # The physical Visium pitch is 100 um, so 100 / measured-pitch is
      # microns-per-unit. Rescale the coordinates to microns so that bin_size
      # (100) and dmax (e.g. 250) are physically meaningful regardless of the
      # original unit. When already in microns, um_per_unit ~ 1 (near identity).
      measured    <- .estimateSpotPitch(xy)
      um_per_unit <- 100 / measured
      inferred    <- if (um_per_unit > 0.8 && um_per_unit < 1.25)
        "micron" else "pixel/other"
      xy         <- xy * um_per_unit
      eff_pitch  <- 100
      coord_unit <- "micron"
      message(sprintf(
        paste0("array_col/array_row not found; using spatialCoords directly. ",
               "Measured pitch = %.1f units (inferred %s, %.4f um/unit); ",
               "rescaled to microns, so bin_size = 100 and dmax are in um."),
        measured, inferred, um_per_unit))
    } else {
      # User forced the pitch: work in the coordinates' native units, no rescale.
      eff_pitch  <- spot_pitch
      coord_unit <- NA_character_
      if (verbose)
        message("array_col/array_row not found; using spatialCoords directly ",
                "with supplied spot_pitch = ", eff_pitch, " (native units).")
    }

    x_ideal <- as.numeric(xy[, 1])
    y_ideal <- as.numeric(xy[, 2])
  }

  # Bin geometry. Pointy-top hexagons of circumradius eff_pitch/sqrt(3)
  # tessellate the lattice exactly; only their centroids affect blisa results.
  if (hexagons) {
    r       <- eff_pitch / sqrt(3)
    ang     <- (seq_len(6) - 1) * (pi / 3) + pi / 6
    cospart <- r * cos(ang)
    sinpart <- r * sin(ang)
    geom <- sf::st_sfc(lapply(seq_along(x_ideal), function(i) {
      vx <- x_ideal[i] + cospart
      vy <- y_ideal[i] + sinpart
      sf::st_polygon(list(cbind(c(vx, vx[1]), c(vy, vy[1]))))
    }))
  } else {
    geom <- sf::st_sfc(lapply(seq_along(x_ideal), function(i)
      sf::st_point(c(x_ideal[i], y_ideal[i]))))
  }

  total_counts <- as.numeric(Matrix::colSums(counts))
  bins <- sf::st_sf(bin_id       = seq_along(x_ideal),
                    n_cells      = 1L,
                    total_counts = total_counts,
                    img_x        = as.numeric(img_xy[, 1]),
                    img_y        = as.numeric(img_xy[, 2]),
                    geometry     = geom)

  # Drop empty spots if requested (keep counts and bins in lockstep)
  if (min_total_counts > 0) {
    keep <- total_counts >= min_total_counts
    if (verbose)
      message(sum(!keep), " spots dropped (< ", min_total_counts,
              " total counts).")
    bins   <- bins[keep, , drop = FALSE]
    counts <- counts[, keep, drop = FALSE]
  }

  list(counts_matrix = counts, bins = bins,
       pitch = eff_pitch, coord_unit = coord_unit)
}
