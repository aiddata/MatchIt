#Script that returns a PSM weight
#based on the given distance-decay model.

thresh_function <- function(x, thresh) {
  if (x > thresh) {
    NA
  } else {
    1*x
  }
}

sph_function <- function(x, thresh) {
  (3/2) * (x/thresh) - (1/2) * (x/thresh)^3
}

gauss_function <- function(x, thresh) {
  exp (-(x^2/thresh^2))
}

exp_function <- function(x, thresh) {
  exp(-(x/thresh))
}

spatial.effects.distance.decay <- function(thresh, model, dist) {

  if (model == "threshold") {

    if (class(dist) == "matrix") {
      weights <- apply(dist, MARGIN=c(1, 2), FUN=thresh_function,
                       thresh=thresh)
    } else {
      weights <- sapply(dist, FUN=thresh_function, thresh=thresh)
    }

  } else if (model == "spherical") {

    if (class(dist) == "matrix") {
      weights <- apply(dist, MARGIN=c(1, 2), FUN=sph_function, thresh=thresh)
    } else {
      weights <- sapply(dist, FUN=sph_function, thresh=thresh)
    }

  } else if (model == "GausSemiVar") {

    if (class(dist) == "matrix") {
      weights <- apply(dist, MARGIN=c(1, 2), FUN=gauss_function, thresh=thresh)
    } else {
      weights <- sapply(dist, FUN=gauss_function, thresh=thresh)
    }

  } else if (model == "ExpSemiVar") {

    if (class(dist) == "matrix") {
      weights <- apply(dist, MARGIN=c(1, 2), FUN=exp_function, thresh=thresh)
    } else {
      weights <- sapply(dist, FUN=exp_function, thresh=thresh)
    }

  }

  return(weights)
}
