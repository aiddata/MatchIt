# Script that returns a PSM weight
# based on the given distance-decay model.

distance.decay.threshold <- function(x, thresh) {
  if (x < thresh) {
    NA
  } else {
    1*x
  }
}

distance.decay.spherical <- function(x, thresh) {
  (3/2) * (x/thresh) - (1/2) * (x/thresh)^3
}

distance.decay.gaussian.semivariance <- function(x, thresh) {
  exp (-(x^2/thresh^2))
}

distance.decay.exponential.semivariance <- function(x, thresh) {
  exp(-(x/thresh))
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
