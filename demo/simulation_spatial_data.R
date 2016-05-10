
# -----------------------------------------------------------------------------
#Generating Dataframe
# -----------------------------------------------------------------------------

# generate random coordinates
random.longitude <- runif(nrandom, minx, maxx)
random.latitude <- runif(nrandom, miny, maxy)

# create dataframe
spdf <- data.frame(id = 1:nrandom,
                   longitude = random.longitude,
                   latitude = random.latitude)

# convert to spatial points dataframe
coordinates(spdf) <- c("longitude", "latitude")
#proj4string(spdf) <- CRS("+proj=longlat +datum=WGS84")
#prj <- proj4string(spdf)
#spdf <- spTransform(spdf, CRS(prj))

# -----------------------------------------------------------------------------
#Generating Ancillary Data
# -----------------------------------------------------------------------------


# using gstat to generate fields with spatial autocorrelation
# source:
#   http://santiago.begueria.es/2010/10/
#     generating-spatially-correlated-random-fields-with-r/

# Define the spatial correlation of the x-variable
var.g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                     model=vgm(psill=xvar_psill, model="Sph", var1.vrange),
                     nmax=20)

# make simulations based on the gstat object
var.sim <- predict(var.g.dummy, newdata=spdf, nsim=1)

spdf@data$trueVar <- var.sim$sim1

# -----------------------------------------------------------------------------
#Generating Ancillary Error Data
# -----------------------------------------------------------------------------

# Define the spatial correlation of the x-variable
var.error.g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                           model=vgm(psill=xvar_error_psill, model="Sph", var1_error.vrange),
                           nmax=20)

# make simulations based on the gstat object
var.error.sim <- predict(var.error.g.dummy, newdata=spdf, nsim=1)

spdf@data$modelVarError <- var.error.sim$sim1* (1-prop_acc)
spdf@data$modelVar <- (var.error.sim$sim1 * (1-prop_acc)) + var.sim$sim1

# -----------------------------------------------------------------------------
#Generating General Model Error Data
# -----------------------------------------------------------------------------

# Define the spatial correlation of the x-variable
mod.error.g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                           model=vgm(psill=mod_error_psill, model="Sph", mod_error.vrange),
                           nmax=20)

# make simulations based on the gstat object
mod.error.sim <- predict(mod.error.g.dummy, newdata=spdf, nsim=1)

spdf@data$modelError <- mod.error.sim$sim1 * mod_error.magnitude

# -----------------------------------------------------------------------------
#Generating the Treatment Variable
# -----------------------------------------------------------------------------

treatment.binary <- ifelse(spdf@data$trueVar > quantile(spdf@data$trueVar, 
                                                        (1-trt_prc)), 1, 0)
spdf$treatment.status = treatment.binary




# -----------------------------------------------------------------------------
#Outcome Model
# -----------------------------------------------------------------------------

# Define the spatial spillover variogram
trt.spillover <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                       model=vgm(psill=trt_spill_sill, model="Sph", spill.vrange),
                       nmax=20)
trt.spillover.sim <- predict(trt.spillover, newdata=spdf, nsim=1)
vgm.spillover <- variogram(sim1~1, trt.spillover.sim, covariogram=TRUE)



vgm.polynomial <- lm(gamma ~ poly(dist, degree=5, raw=TRUE),
                     data=vgm.spillover)

vgm.spillover$fit <- predict(vgm.polynomial, vgm.spillover)




#Calculate the distance to treated units for each unit (excluding itself)
tmp.spillover.weights <- c()
for (i in 1:nrandom) {
  
  tmp.dist <- spDists(spdf[i, ]@coords, spdf@coords)
  
  
  tmp.neighbors <- tmp.dist > 0
  tmp.treated <- spdf$treatment.status
  
  tmp.newdata <- data.frame(
    dist = c(tmp.dist)
  )
  
  #Limit to the min (after the final predicted distance all values are 0)
  tmp.newdata$gamma <- predict(vgm.polynomial, newdata=tmp.newdata)
  tmp.newdata$gamma[tmp.newdata$gamma < 0] <- 0
  
  tmp.spillover.weights[i] <- sum(tmp.newdata$gamma * tmp.neighbors * tmp.treated)
  
}

vgm.spillover$gamma[vgm.spillover$gamma < 0] <- 0

tmp.spillover.weights[is.na(tmp.spillover.weights)] <- 0

#Total treated effect across the study area
tmp.spillover.t1 <- sum(theta * spdf$treatment.status)

#Ratio of spillover to assign to each unit
tmp.spillover.t2 <- c()
sum.spillover.weights <- sum(tmp.spillover.weights, na.rm=TRUE)
for (i in 1:nrandom) {
  tmp.spillover.t2[i] <- tmp.spillover.weights[i] / sum.spillover.weights
  if(is.na(tmp.spillover.t2[i])){
    tmp.spillover.t2[i] = 0
  }
}

#Calculate total spillover
tot_spill <- tmp.spillover.t1 * theta * spill.magnitude

#Calculate per-unit spillover
spdf$trueSpill <-tmp.spillover.t2 * tot_spill

#Calculate per-unit true outcome
spdf$trueOutcome <- spdf$trueSpill + (beta * spdf$trueVar) + spdf$treatment.status

spdf$modelOutcome <- spdf$trueOutcome + spdf$modelError





