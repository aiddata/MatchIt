
library(devtools)

detach("package:MatchIt", unload=TRUE)
load_all("~/git/matchit/R")
# library(devtools)
# install_github("itpir/matchit")
library(MatchIt)

library(sp)
library(ncf)
library(gstat)


# -----------------------------------------------------------------------------
# options

# set dataframe size (number of points)
nrandom <- 1000

# variogram model range (larger = coarser autocorrelation)
vrange <- 1500

# numner of covariates
ncovariates <- 3

# define bounding box
minx <- -45
maxx <- 45
miny <- -22.5
maxy <- 22.5

# correlogram increment size
correlogram.increment <- 500

# -----------------------------------------------------------------------------


# generate random treatment and coordinates
random.treatment <- rbinom(nrandom, 1, 0.5)
random.longitude <- runif(nrandom, minx, maxx)
random.latitude <- runif(nrandom, miny, maxy)

# create dataframe
spdf <- data.frame(treat = random.treatment,
                   longitude = random.longitude,
                   latitude = random.latitude)

# convert to spatial points dataframe
coordinates(spdf) <- c("longitude", "latitude")
proj4string(spdf) <- CRS("+proj=longlat +datum=WGS84")
prj <- proj4string(spdf)
spdf <- spTransform(spdf, CRS(prj))


# using gstat to generate fields with spatial autocorrelation
# source:
# http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/

# define the gstat object (spatial model)
g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, 
                 model=vgm(psill=0.1, model="Sph", range=vrange), 
                 nmax=20)

# make simulations based on the gstat object
# where each simulation will be a variable
sim <- predict(g.dummy, newdata=spdf, nsim=ncovariates)

# add covariates from sim to spatial dataframe
for (i in 1:ncovariates) {
  spdf@data[[paste("var", i, sep="")]] <- sim@data[[paste("sim", i, sep="")]]
}

# plot spatial autocorrelation of covariates
spplot(spdf)



vgm1 <- variogram(var1~1, spdf)
plot(vgm1)


# model.1 <- fit.variogram(vgm1,vgm(1,"Sph",300,1))
# plot(vgm1, model=model.1)
# plot(vgm1, plot.numbers = TRUE, pch = "+")
# vgm2 <- variogram(log(zinc)~1, meuse, alpha=c(0,45,90,135))
# plot(vgm2)
# # the following demonstrates plotting of directional models:
# model.2 <- vgm(.59,"Sph",926,.06,anis=c(0,0.3))
# plot(vgm2, model=model.2)



# -----------------------------------------------------------------------------


# traditional, non-spatial matchit
m1.out <- matchit(treat ~ var1 + var2 + var3, data=spdf@data,
                  method="nearest", distance="logit", caliper=0.25)


# get correlogram for PSM distance, x-intercept and plot
m1.correlogram.data <- correlog(x=spdf@coords[, 1],
                                y=spdf@coords[, 2],
                                z=m1.out$distance, 
                                increment=correlogram.increment, latlon=TRUE,
                                na.rm=TRUE, resamp=50, quiet=TRUE)

m1.correlogram.xintercept <- as.numeric(m1.correlogram.data$x.intercept)

plot.correlog(m1.correlogram.data)
mtext("PSM pscores")

m1.covariate.correlogram.xintercept <- c()
for (i in 1:ncovariates) {

  tmp.m1.correlogram.data <- correlog(x=spdf@coords[, 1],
                                      y=spdf@coords[, 2],
                                      z=spdf@data[[paste("var", i, sep="")]], 
                                      increment=correlogram.increment, latlon=TRUE,
                                      na.rm=TRUE, resamp=50, quiet=TRUE)
  
  m1.covariate.correlogram.xintercept[i] <- as.numeric(tmp.m1.correlogram.data$x.intercept)
  
  plot(tmp.m1.correlogram.data)
  mtext(paste("var", i, sep=""))
  
}

# -------------------------------------

m1.match.distances <- c()
  
for (i in 1:length(m1.out$match.matrix)) {
  control.id <- as.numeric(labels(m1.out$match.matrix)[[1]][i])
  treated.id <- m1.out$match.matrix[i]
  
  if (!is.na(treated.id)) {
    treated.id <- as.numeric(treated.id)
    
    control.coords <- spdf[control.id, ]@coords
    treated.coords <- spdf[treated.id, ]@coords
    
    m1.match.distances[i] <- spDists(control.coords, treated.coords, 
                                     longlat=FALSE, segments=FALSE, 
                                     diagonal=FALSE)
  }
}

m1.autocorrelation.match.count <- length(m1.match.distances[!is.na(m1.match.distances) & m1.match.distances < (m1.correlogram.xintercept/110)])


# -----------------------------------------------------------------------------


# spatial matchit
spatial.opts <- list(decay.model = "gaussian.semivariance",
                     threshold = 0.001)

# spatial.opts <- list(decay.model = "threshold",
#                      threshold = (m1.correlogram.xintercept/110))

m2.out <- matchit(treat ~ var1 + var2 + var3, data=spdf,
                  method = "nearest", distance = "logit", caliper=0.25,
                  spatial.options=spatial.opts)


# # get correlogram for PSM distance, x-intercept and plot
# m2.correlogram.data <- correlog(x=spdf@coords[, 1],
#                                 y=spdf@coords[, 2],
#                                 z=m2.out$distance, 
#                                 increment=correlogram.increment, latlon=TRUE, 
#                                 na.rm=TRUE, resamp=50, quiet=TRUE)
# 
# m2.correlogram.xintercept <- as.numeric(m2.correlogram.data$x.intercept)
# 
# plot.correlog(m2.correlogram.data)

# -------------------------------------

m2.match.distances <- c()

for (i in 1:length(m2.out$match.matrix)) {
  control.id <- as.numeric(labels(m2.out$match.matrix)[[1]][i])
  treated.id <- m2.out$match.matrix[i]
  
  if (!is.na(treated.id)) {
    treated.id <- as.numeric(treated.id)
    
    control.coords <- spdf[control.id, ]@coords
    treated.coords <- spdf[treated.id, ]@coords
    
    m2.match.distances[i] <- spDists(control.coords, treated.coords, 
                                     longlat=FALSE, segments=FALSE, 
                                     diagonal=FALSE)
  }
}

m2.autocorrelation.match.count <- length(m2.match.distances[!is.na(m2.match.distances) & m2.match.distances < (m1.correlogram.xintercept/110)])





