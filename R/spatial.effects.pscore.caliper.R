#Calculates the full matrix of spatially-weighted PSMs, given
#a dataframe containing a vector of PSM scores (from which deviations are calculated)
#and two vectors with latitude and longitude.
#Returns the matrix-wide, distance-weighted caliper.
#caliper*sqrt(var(distance[in.sample==1]))
spatial.effects.pscore.caliper <- function(spatial.thresholds, spatial.decay.model, spatial.data, distance, caliper, treat)
{

  #Select the treated being analyzed to calculate distance penalties from
  #treated_unit <- spatial.data[rownames(spatial.data@data)==itert,]

  spatial.data@data$distance <- distance
  treated <- spatial.data[names(distance[treat==1]),]
  untreated <- spatial.data[names(distance[treat==0]),]

  trt_matrix <- matrix(treated@data$distance, length(treated), length(untreated), byrow=FALSE)
  untrt_matrix <- matrix(untreated@data$distance, length(treated), length(untreated), byrow=TRUE) 
  
  #raw PSM deviation matrix
  psm_dev_matrix = abs(trt_matrix - untrt_matrix)
  
  #Make an empty copy of the matrix to fill with geographic distances
  geog_dist_matrix <- psm_dev_matrix 
  geog_dist_matrix[geog_dist_matrix != 0] <- NA
  
  #Calculate the geographic distances between points
  for (i in 1:length(treated))
  {
    geog_dist_matrix[i,]<-spDistsN1(untreated, treated[i,])
  }
  
  #Permute the geographic distances by the spatial distance-decay function.
  W <- spatial.effects.distance.decay(thresh=spatial.thresholds, model=spatial.decay.model, dist=geog_dist_matrix)
  
  #calculate the weighted P-scores
  W_P <- W * psm_dev_matrix
  

  caliper_cut <- caliper * sqrt(var(as.vector(W_P)))
  return(caliper_cut)
}