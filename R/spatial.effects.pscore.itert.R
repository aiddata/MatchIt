#Spatial penalty function
#Takes in a spatial dataframe, origin (treatment) case
#Vector of calculated, un-spatially adjusted propensity scores
#And an ID of the treatment case the vector belong to.
#Outputs an adjusted vector of the same type with P-scores adjusted to account for
#spatial autocorrelation by penalizing along a distance-decay function.
#Further, a lower-bounds threshold is applied to mitigate the potential for spillovers in treatments.

spatial.effects.pscore.itert <- function(spatial.thresholds, spatial.decay.model, spatial.data, deviation, itert)
{
  #Select the treated being analyzed to calculate distance penalties from
  treated_unit <- spatial.data[rownames(spatial.data@data)==itert,]
  
  #Identify candidates to match with
  pair_candidates <- spatial.data[names(deviation),]
  pair_candidates@data$unstd_deviation <- deviation
  
  #Calculate the geographic distances between points
  geog_dist_vector<-spDistsN1(pair_candidates, treated_unit)
  
  
  #Permute the geographic distances by the spatial distance-decay function.
  W <- spatial.effects.distance.decay(thresh=spatial.thresholds, model=spatial.decay.model, dist=geog_dist_vector)
  
  #calculate the weighted P-scores
  W_P <- W * deviation
  return(W_P)
}