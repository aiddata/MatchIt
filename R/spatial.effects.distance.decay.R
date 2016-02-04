#Script that returns a PSM weight
#based on the given distance-decay model.

spatial.effects.distance.decay <- function(thresh, model, dist)
{
  if(model=="threshold")
  {
    thresh_function <- function(x){if(x<thresh){0*x}else{1*x}}
    if(class(dist) == "matrix")
    {
    weights <- apply(dist,MARGIN=c(1,2),FUN=thresh_function)
    }
    else
    {
      weights <- sapply(dist, FUN=thresh_function) 
    }
    return(weights)
  }
}