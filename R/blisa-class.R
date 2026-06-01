#' Create a blisa object
#'
#' Constructor for the \code{blisa} S3 class, which stores the full output of
#' \code{\link{blisa}}.
#'
#' @param LR_results Data frame of BLISA results, one row per ligand-receptor
#'   pair.
#' @param bins An \code{sf} object of hexagonal bin polygons.
#' @param spatial_weights Spatial weights list from
#'   \code{\link{computeSpatialWeights}}.
#' @param CCI_scores Wide data frame of cell-cell interaction scores from
#'   \code{\link{runCCI}}, or \code{NULL} if CCI was not computed.
#'
#' @return An object of class \code{blisa}.
#' @keywords internal
new_blisa <- function(LR_results, bins, spatial_weights, CCI_scores = NULL) {
  structure(
    list(
      LR_results      = LR_results,
      bins            = bins,
      spatial_weights = spatial_weights,
      CCI_scores      = CCI_scores
    ),
    class = "blisa"
  )
}


#' @export
print.blisa <- function(x, ...) {
  cat("A blisa object\n")
  cat(" LR pairs tested  :", nrow(x$LR_results), "\n")
  cat(" Significant pairs:", sum(x$LR_results$sig_numbers > 0), "\n")
  cat(" Bins             :", nrow(x$bins), "\n")
  cat(" CCI computed     :", !is.null(x$CCI_scores), "\n")
  invisible(x)
}


#' Test if an object is a blisa object
#'
#' @param x Any R object.
#' @return Logical.
#' @examples
#' is.blisa(list())           # FALSE
#' is.blisa("not a blisa")    # FALSE
#' @export
is.blisa <- function(x) inherits(x, "blisa")
