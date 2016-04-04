# Script that returns a PSM weight
# based on the given distance-decay model.

distance.decay.threshold <- function(x, thresh) {
  if (x < thresh) {
    y <- 1
  } else {
    y <- 0
  }
  return(y)
}

distance.decay.linear.threshold <- function(x, thresh) {
  if (x < thresh) {
    y <- 1
  } else {
    y <- thresh/x
  }
  return(y)
}

distance.decay.spherical <- function(x, thresh) {
  if (x < thresh) {
    y <- 1 - ((3/2) * (x/thresh) - (1/2) * (x/thresh)^3)
  } else {
    y <- 0
  }
  return(y)
}

distance.decay.gaussian.semivariance <- function(x, thresh) {
  y <- exp(-(x^2 / thresh^2))
  return(y)
}

distance.decay.exponential.semivariance <- function(x, thresh) {
  y <- exp(-(x / thresh))
  return(y)
}

# # weights values based off of correlogram using distance
# # this might need additional arguments for correl/dist
# distance.decay.morans <- function(x, thresh, correl?, dist?) {
#
# }


run.distance.decay <- function(thresh, dist, func) {

  if (class(dist) == "matrix") {
    weights <- apply(dist, MARGIN=c(1, 2), FUN=func, thresh=thresh)
  } else {
    weights <- sapply(dist, FUN=func, thresh=thresh)
  }

  return(weights)
}
