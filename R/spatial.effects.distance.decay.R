#Script that returns a PSM weight
#based on the given distance-decay model.

threshold.function <- function(x, thresh) {
  if (x > thresh) {
    NA
  } else {
    1*x
  }
}

spherical.function <- function(x, thresh) {
  (3/2) * (x/thresh) - (1/2) * (x/thresh)^3
}

gaussian.function <- function(x, thresh) {
  exp (-(x^2/thresh^2))
}

exponential.function <- function(x, thresh) {
  exp(-(x/thresh))
}


distance.decay.threshold <- function(thresh, dist) {
  if (class(dist) == "matrix") {
    weights <- apply(dist, MARGIN=c(1, 2), FUN=threshold.function,
                     thresh=thresh)
  } else {
    weights <- sapply(dist, FUN=threshold.function, thresh=thresh)
  }
  return(weights)
}

distance.decay.spherical <- function(thresh, dist) {
  if (class(dist) == "matrix") {
    weights <- apply(dist, MARGIN=c(1, 2), FUN=spherical.function, thresh=thresh)
  } else {
    weights <- sapply(dist, FUN=spherical.function, thresh=thresh)
  }
  return(weights)
}

distance.decay.gaussian.semivariance <- function(thresh, dist) {
  if (class(dist) == "matrix") {
    weights <- apply(dist, MARGIN=c(1, 2), FUN=gaussian.function, thresh=thresh)
  } else {
    weights <- sapply(dist, FUN=gaussian.function, thresh=thresh)
  }
  return(weights)
}

distance.decay.exponential.semivariance <- function(thresh, dist) {
  if (class(dist) == "matrix") {
    weights <- apply(dist, MARGIN=c(1, 2), FUN=exponential.function, thresh=thresh)
  } else {
    weights <- sapply(dist, FUN=exponential.function, thresh=thresh)
  }
  return(weights)
}

# distance.decay.exponential.morans <- function(thresh, dist) {

# }



# spatial.effects.distance.decay <- function(thresh, model, dist) {

#   if (model == "threshold") {

#     if (class(dist) == "matrix") {
#       weights <- apply(dist, MARGIN=c(1, 2), FUN=threshold.function,
#                        thresh=thresh)
#     } else {
#       weights <- sapply(dist, FUN=threshold.function, thresh=thresh)
#     }

#   } else if (model == "spherical") {

#     if (class(dist) == "matrix") {
#       weights <- apply(dist, MARGIN=c(1, 2), FUN=spherical.function, thresh=thresh)
#     } else {
#       weights <- sapply(dist, FUN=spherical.function, thresh=thresh)
#     }

#   } else if (model == "GausSemiVar") {

#     if (class(dist) == "matrix") {
#       weights <- apply(dist, MARGIN=c(1, 2), FUN=gaussian.function, thresh=thresh)
#     } else {
#       weights <- sapply(dist, FUN=gaussian.function, thresh=thresh)
#     }

#   } else if (model == "ExpSemiVar") {

#     if (class(dist) == "matrix") {
#       weights <- apply(dist, MARGIN=c(1, 2), FUN=exponential.function, thresh=thresh)
#     } else {
#       weights <- sapply(dist, FUN=exponential.function, thresh=thresh)
#     }
#   } # else if (model == "morans") {
#   # }

#   return(weights)
# }
