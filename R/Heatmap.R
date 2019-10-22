
### creat a matrix of x,y
createMap <- function(smooth.data,
                      x.vec = NULL,
                      y.vec = NULL) {
  coord.long <- data.frame('x' = smooth.data$x,
                           'y' = smooth.data$y) %>%
    round(. , digits = 0) %>%
    dplyr::group_by(x, y) %>%
    dplyr::summarise(count = n())
  if (is.null(x.vec) & is.null(y.vec)) {
    x.vec <- seq(min(coord.long$x), max(coord.long$x), 1)
    y.vec <- seq(min(coord.long$y), max(coord.long$y), 1)
  }
  full.mat <- expand.grid('x' = x.vec, 'y' = y.vec)
  full.mat <- dplyr::left_join(full.mat, coord.long) %>% dplyr::mutate(count = replace(count, which(is.na(count)), 0))
  return(xtabs(count ~ x + y, data = full.mat))
}

smoothMap <- function(wide.coord.mat,
                      sigma.para   = 9,
                      n            = 63) {
  x.vec <- as.numeric(row.names(wide.coord.mat))
  y.vec <- as.numeric(colnames(wide.coord.mat))
  wide.heat.smooth <- smoothie:::kernel2dsmooth(wide.coord.mat, kernel.type ="gauss", nx = n, ny = n, sigma = sigma.para)
  colnames(wide.heat.smooth) <- y.vec
  rownames(wide.heat.smooth) <- x.vec
  return(wide.heat.smooth)
}


### Find maxima
FindMax <- function(wide.heat.smooth, quant) {
  y.dim <- ncol(wide.heat.smooth)
  x.dim <- nrow(wide.heat.smooth)
  wide.heat.smooth[wide.heat.smooth < quantile(wide.heat.smooth, quant)] <- 0
  M1 <- cbind(0,wide.heat.smooth[ ,-y.dim])
  M2 <- cbind(wide.heat.smooth[ ,-1], 0)
  M3 <- rbind(0, wide.heat.smooth[-x.dim, ])
  M4 <- rbind(wide.heat.smooth[-1, ], 0)
  M5 <- rbind(0,M1[-1, ])
  M6 <- rbind(0, M2[-1, ])
  M7 <- rbind(M1[-y.dim,],0)
  M8 <- rbind(M2[-1,], 0);
  max.point <- which(wide.heat.smooth >= M1 & wide.heat.smooth > M2 & wide.heat.smooth > M3 &
                       wide.heat.smooth > M4 & wide.heat.smooth > M5 & wide.heat.smooth > M6 &
                       wide.heat.smooth > M7 & wide.heat.smooth > M8, arr.ind = TRUE)
  y.coord.max <- as.numeric(colnames(wide.heat.smooth)[max.point[ ,2]])
  x.coord.max <- as.numeric(rownames(wide.heat.smooth)[max.point[ ,1]])
  return(cbind('x' = x.coord.max, 'y' = y.coord.max))
}

#' HeatMaker - Calculates smoothed dwell time at each coodinates.
#' @param smooth.data A data.frame contains x and y coordinates
#'  name of x coordinates vector should contain "x" or "X"
#'   column name of y coordinates vector should contain "y" or "Y"
#' @param x.vec The grid of x-axis for the heatmap. Defaults to seq(min(x), max(x), 1)
#' @param y.vec The grid of y-axis for the heatmap. Defaults to seq(min(y), max(y), 1)
#' @param sigma The standard deviation of the Gaussian smoothing kernel, controls the level of the
#' smoothing (larger, the heatmap will be smoother)
#' @param n The filter size is n * n
#' @param quant Quantile threshold, returns the > quant local maxima
#' @param to.plot Flag, if TRUE will also plot the heatmap
#' @return  The function returns a list with following objects:
  #'  \describe{
  #'  \item{smooth.dwell.time}{ a data.frame, x, y, smooth.dwell.time. For each cell coordinates the smoothed dwell time spent}
  #'  \item{maxpoint}{coordinates of local maxima above the predetermined threshold}}
#' @export

SmoothDwell <- function(smooth.data,
                        x.vec = NULL,
                        y.vec = NULL,
                        sigma = 9,
                        n     = 57,
                        quant = 0.95,
                        to.plot = FALSE) {
  smooth.heat.map <- smoothMap(createMap(smooth.data, x.vec, y.vec), n = n, sigma.para =  sigma)
  max.points      <- FindMax(smooth.heat.map, quant)
  smooth.data.dwell.time <- reshape2::melt(smooth.heat.map)
  colnames(smooth.data.dwell.time) <- c('x', 'y', 'smooth.dwell.time')
  if (to.plot) {
    print(DrawHeat(smooth.data.dwell.time))
  }
  return(list('smooth.dwell.time' = smooth.data.dwell.time,
              'local.maxima'      = max.points))
}

### Ilan favourite pallette
col.pal <- c(scales::seq_gradient_pal('white', 'red')(seq(0, 1, length.out = 20)),
             scales::seq_gradient_pal('red', 'yellow')(seq(0, 1, length.out = 30)),
             scales::seq_gradient_pal('yellow', 'yellow')(10))

### Create heatmap
DrawHeat <- function(smooth.dwell.df) {
  ggplot2::ggplot() +
    ggplot2::geom_tile(data =  smooth.dwell.df, ggplot2::aes(x = x, y = y, fill = smooth.dwell.time)) +
    ggplot2::geom_contour(data =  smooth.dwell.df, ggplot2::aes(x = x, y = y, z = smooth.dwell.time),
                 color = "darkgrey", alpha = 1, bins = 10) +
    ggplot2::scale_fill_gradientn(colors = col.pal, values = seq(0, 37, length.out = 50) / 37)  +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = 'none') +
    ggplot2::xlab('') +
    ggplot2::ylab('') +
    ggplot2::coord_fixed(ratio = 1)
}

