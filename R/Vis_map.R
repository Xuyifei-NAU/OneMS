#' To normalize data to a certin ranges for visualization
#'
#'A convenient function to modify the data distribution range
#'
#' @param x original matrix or vector or data.frame, must be numeric
#' @param range your expected range
#' @param from.range
#'
#' @return matrix or vector or data.frame after mapping to ranges
#' @export
#'
#' @examples
#' Vis_map(matrix(1:9,ncol = 3))
#' Vis_map(1:9)
#' Vis_map(data.frame(a=1:3,b = 7:9))
Vis_map <- function(x, range = c(0,1), from.range=NA) {
  if(any(is.na(from.range))) from.range <- range(x, na.rm=TRUE)

  ## check if all values are the same
  if(!diff(from.range)) return(
    matrix(mean(range), ncol=ncol(x), nrow=nrow(x),
           dimnames = dimnames(x)))

  ## map to [0,1]
  x <- (x-from.range[1])
  x <- x/diff(from.range)
  ## handle single values
  if(diff(from.range) == 0) x <- 0

  ## map from [0,1] to [range]
  if (range[1]>range[2]) x <- 1-x
  x <- x*(abs(diff(range))) + min(range)

  x[x<min(range) | x>max(range)] <- NA

  x
}
