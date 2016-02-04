#Script that returns a PSM weight
#based on the given distance-decay model.

spatial.effects.distance.decay <- function(thresh, model, dist)
{
  if(model=="threshold")
  {
    thresh_function <- function(x){if(x<thresh){NA}else{1*x}}
    if(class(dist) == "matrix")
    {
    weights <- apply(dist,MARGIN=c(1,2),FUN=thresh_function)
    }
    else
    {
      weights <- sapply(dist, FUN=thresh_function) 
    }
   
  }
  if(model=="spherical")
  {
    sph_function <- function(x){(3/2) * (x/thresh) - (1/2) * (x/thresh)^3}
    if(class(dist) == "matrix")
    {
      weights <- apply(dist,MARGIN=c(1,2),FUN=sph_function)
    }
    else
    {
      weights <- sapply(dist, FUN=sph_function) 
    }
  }
  if(model=="GausSemiVar")
  {
    gauss_function <- function(x){1 - exp (-(x^2/thresh^2))}
    if(class(dist) == "matrix")
    {
      weights <- apply(dist,MARGIN=c(1,2),FUN=gauss_function)
    }
    else
    {
      weights <- sapply(dist, FUN=gauss_function) 
    }
  }
  if(model=="ExpSemiVar")
  {
    exp_function <- function(x){1 - exp(-(x/thresh))}
    if(class(dist) == "matrix")
    {
      weights <- apply(dist,MARGIN=c(1,2),FUN=exp_function)
    }
    else
    {
      weights <- sapply(dist, FUN=exp_function) 
    }
  }
  return(weights)
}