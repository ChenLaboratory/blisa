#' Compute Theoretical Visium Spot Centre-to-Centre Distance
#'
#' Calculates the theoretical centre-to-centre distance between adjacent Visium
#' spots in image pixel units using the scale factors stored in a Seurat object.
#' Formula: \code{spot * (100/55) * scale_factor}, where 100 µm is the physical
#' centre-to-centre distance and 55 µm is the spot diameter.
#'
#' @param seu A \code{Seurat} object containing a Visium image with scale factors.
#' @param scale Character. Which scale factor image to use. One of
#'   \code{"lowres"} (default) or \code{"hires"}.
#'
#' @return Numeric scalar: the centre-to-centre spot distance in image pixels.
#' @export
visium_spot_distance <- function(seu, scale = "lowres") {
  sf           <- seu@images[[1]]@scale.factors
  spot_distance <- sf$spot * (100 / 55) * sf[[scale]]
  return(spot_distance)
}
