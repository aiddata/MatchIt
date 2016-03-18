check.is.spatial <- function(data) {

  is.spatial <- FALSE

  slots <- slotNames(data)

  if (class(data) != "data.frame" && "data" %in% slots &&
      "coords" %in% slots && "bbox" %in% slots && "proj4string" %in% slots) {

    is.spatial <- TRUE

  }

  return(is.spatial)

}
