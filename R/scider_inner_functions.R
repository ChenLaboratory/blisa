grid2df <- function(spe,
                    x=spe@metadata$grid_density$node_x,
                    y=spe@metadata$grid_density$node_y,
                    reverseY = FALSE,
                    ...) {
  if (is.null(spe@metadata$grid_info)) stop("Missing grid. Compute Density first")
  if (is.null(x) || is.null(y)) stop("Missing x or y")
  if (length(x) != length(x)) stop("x, y must be of equal length")

  x <- as.numeric(x)
  y <- as.numeric(y)

  nx <- spe@metadata$grid_info$dims[1]
  ny <- spe@metadata$grid_info$dims[2]

  if (spe@metadata$grid_info$grid_type=="hex") {
    dx <- spe@metadata$grid_info$xstep/2
    dy <- spe@metadata$grid_info$ystep/3
    offset <- c(spe@metadata$grid_info$xlim[1] - dx,
                spe@metadata$grid_info$ylim[1] - dy * 2)

    xc <- offset[1] + (0:(nx*2+1)) * dx
    yc <- offset[2] + (0:(ny*3+1)) * dy

    n_vertices = 7
    x_pattern = c(0,1,1,0,-1,-1,0)
    y_pattern = c(2,1,-1,-2,-1,1,2)

    x <- 2*x + !y%%2
    y <- 3*y
  } else { # square
    dx <- c(spe@metadata$grid_info$xstep,spe@metadata$grid_info$ystep)
    offset <- c(spe@metadata$grid_info$xlim[1],spe@metadata$grid_info$ylim[1])

    xc <- offset[1] + (0:nx) * dx[1]
    yc <- offset[2] + (0:ny) * dx[2]

    n_vertices = 5
    x_pattern = c(0,1,1,0,0)
    y_pattern = c(0,0,1,1,0)
  }
  x <- xc[rep(x, each = n_vertices) + x_pattern]
  y <- yc[rep(y, each = n_vertices) + y_pattern]
  L2 <- rep(seq_len(length(x)/n_vertices),each=n_vertices)
  res <- data.frame(X=x,Y=y,L2=L2)

  lst <- lapply(list(...),function(ii) rep(ii,each=n_vertices))
  if (length(lst)) res <- cbind(res,lst)
  return(res)
}

update_bound <- function(p,x=NULL,y=NULL) {
  # Whether reverseY has been applied to p
  reverseY <- any(sapply(p$scales$scales, function(s) {
    inherits(s, "ScaleContinuousPosition") &&
      "y" %in% s$aesthetics &&
      "reverse" %in% s$trans$name
  }))
  p$coordinates$limits$x = range(p$coordinates$limits$x,x)
  p$coordinates$limits$y = range(p$coordinates$limits$y,y)

  if (reverseY) p$coordinates$limits$y <- rev(p$coordinates$limits$y)
  return(p)
}

