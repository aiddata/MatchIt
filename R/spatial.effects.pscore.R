# Calculates the full matrix of spatially-weighted PSMs, given
# a dataframe containing a vector of PSM scores (from which deviations are calculated)
# and two vectors with latitude and longitude.
# Returns the matrix-wide, distance-weighted caliper.
# caliper*sqrt(var(distance[in.sample==1]))

spatial.effects.pscore.caliper <- function(spatial.threshold,
                                           spatial.decay.function,
                                           spatial.data,
                                           distance, caliper, treat) {

  # Select the treated being analyzed to calculate distance penalties from
  # treated_unit <- spatial.data[rownames(spatial.data@data)==t.iter.label,]

  spatial.data@data$distance <- distance
  treated <- spatial.data[names(distance[treat == 1]),]
  untreated <- spatial.data[names(distance[treat == 0]),]

  trt.matrix <- matrix(treated@data$distance, nrow=length(treated),
                       ncol=length(untreated), byrow=FALSE)

  untrt.matrix <- matrix(untreated@data$distance, nrow=length(treated),
                         ncol=length(untreated), byrow=TRUE)

  # raw PSM deviation matrix
  # psm.dev.matrix <- abs(trt.matrix - untrt.matrix)

  # # Make an empty copy of the matrix to fill with geographic distances
  # geog.dist.matrix <- psm.dev.matrix
  # geog.dist.matrix[geog.dist.matrix != 0] <- NA

  geog.dist.matrix <- matrix(NA, nrow=length(treated), ncol=length(untreated))

  # Calculate the geographic distances between points
  for (i in 1:length(treated)) {
    geog.dist.matrix[i,] <- spDistsN1(untreated, treated[i,], longlat=TRUE)
  }

  # Permute the geographic distances by the spatial distance-decay function.
  spatial.weights <- run.distance.decay(thresh=spatial.threshold,
                                        dist=geog.dist.matrix,
                                        func=spatial.decay.function)

  # # calculate the weighted P-scores
  # spatial.weighted.pscores <- spatial.weights * psm.dev.matrix

  # return(as.vector(spatial.weighted.pscores))


  return(c(rowMeans((spatial.weights + trt.matrix) / 2, na.rm=TRUE),
           colMeans((spatial.weights + untrt.matrix) / 2, na.rm=TRUE)))

}


# Spatial penalty function
#
# Takes in a spatial dataframe, origin (treatment) case
# Vector of calculated, un-spatially adjusted propensity scores
# And an ID of the treatment case the vector belong to.
# Outputs an adjusted vector of the same type with P-scores adjusted to
# account for spatial autocorrelation by penalizing along a distance-decay
# function. Further, a lower-bounds threshold is applied to mitigate the
# potential for spillovers in treatments.

spatial.effects.pscore.deviation <- function(spatial.threshold,
                                             spatial.decay.function,
                                             spatial.data,
                                             deviation, t.iter.label) {

  # Select the treated being analyzed to calculate distance penalties from
  # treated.unit <- spatial.data[rownames(spatial.data@data) == t.iter.label, ]

  # treated.units <- spatial.data[spatial.data$treament.status == 1, ]


  # Identify candidates to match with
  control.candidates <- spatial.data[names(deviation), ]
  control.candidates$unstd.deviation <- deviation

  # Calculate the geographic distances between points
  # geog.dist.vector <- spDistsN1(control.candidates, treated.unit, longlat=TRUE)

  # geog.dist.vector <- c()
  # for (i in 1:nrow(control.candidates)){
  #   tmp.control <- control.candidates[1, ]
  #   geog.dist.vector[i] <- min(spDistsN1(treated.units, tmp.control, longlat=TRUE))
  # }

  geog.dist.vector <- control.candidates$geog.dist.vector

  # Permute the geographic distances by the spatial distance-decay function.
  spatial.weights <- run.distance.decay(thresh=spatial.threshold,
                                        dist=geog.dist.vector,
                                        func=spatial.decay.function)

  # calculate the weighted P-scores
  spatial.weighted.pscores <- (spatial.weights + deviation) / 2

  return(spatial.weighted.pscores)

}



