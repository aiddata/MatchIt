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
  # treated_unit <- spatial.data[rownames(spatial.data@data)==itert,]

  spatial.data@data$distance <- distance
  treated <- spatial.data[names(distance[treat == 1]),]
  untreated <- spatial.data[names(distance[treat == 0]),]

  trt.matrix <- matrix(treated@data$distance, length(treated),
                       length(untreated), byrow=FALSE)

  untrt.matrix <- matrix(untreated@data$distance, length(treated),
                         length(untreated), byrow=TRUE)

  # raw PSM deviation matrix
  psm.dev.matrix <- abs(trt.matrix - untrt.matrix)

  # Make an empty copy of the matrix to fill with geographic distances
  geog.dist.matrix <- psm.dev.matrix
  geog.dist.matrix[geog.dist.matrix != 0] <- NA

  # Calculate the geographic distances between points
  for (i in 1:length(treated)) {
    geog.dist.matrix[i,] <- spDistsN1(untreated, treated[i,])
  }

  # Permute the geographic distances by the spatial distance-decay function.
  spatial.weights <- do.call(spatial.decay.function,
                             list(thresh=spatial.threshold,
                                  dist=geog.dist.vector))


  # calculate the weighted P-scores
  spatial.weighted.pscores <- spatial.weights * psm.dev.matrix

  return(as.vector(spatial.weighted.pscores))

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

spatial.effects.pscore.itert <- function(spatial.threshold,
                                         spatial.decay.function, spatial.data,
                                         deviation, itert) {

  # Select the treated being analyzed to calculate distance penalties from
  treated.unit <- spatial.data[rownames(spatial.data@data) == itert,]

  # Identify candidates to match with
  pair.candidates <- spatial.data[names(deviation),]
  pair.candidates@data$unstd_deviation <- deviation

  # Calculate the geographic distances between points
  geog.dist.vector <- spDistsN1(pair.candidates, treated.unit)

  # Permute the geographic distances by the spatial distance-decay function.
  spatial.weights <- do.call(spatial.decay.function,
                             list(thresh=spatial.threshold,
                                  dist=geog.dist.vector))

  # calculate the weighted P-scores
  spatial.weighted.pscores <- spatial.weights * deviation
  return(spatial.weighted.pscores)
}
