
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


# generate random coordinates
random.longitude <- runif(nrandom, minx, maxx)
random.latitude <- runif(nrandom, miny, maxy)

# create dataframe
spdf <- data.frame(longitude = random.longitude,
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

# initialize vars for generating treatment var
tmp.var <- rep(0, nrandom)
tmp.random <- runif(nrandom, min=0, max=1)

# add covariates from sim to spatial dataframe
for (i in 1:ncovariates) {
  spdf@data[[paste("var", i, sep="")]] <- sim@data[[paste("sim", i, sep="")]]

  tmp.cov <- spdf@data[[paste("var", i, sep="")]]
  tmp.var <- tmp.var + tmp.cov / max(tmp.cov) 
}

# incorporate covariate distribution into treatment
tmp.treat.val <- tmp.random + tmp.var / ncovariates
treatment.binary <- ifelse(tmp.treat.val > quantile(tmp.treat.val, 0.50), 1, 0)
spdf$treat = treatment.binary



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


# ====================================

# get correlogram for PSM distance, x-intercept and plot
correlogram.data <- correlog(x=spdf@coords[, 1],
                                y=spdf@coords[, 2],
                                z=m1.out$distance, 
                                increment=correlogram.increment, latlon=TRUE,
                                na.rm=TRUE, resamp=50, quiet=TRUE)

correlogram.xintercept <- as.numeric(correlogram.data$x.intercept)

plot.correlog(correlogram.data)
mtext("PSM pscores")

# covariate.correlogram.xintercept <- c()
# for (i in 1:ncovariates) {
# 
#   tmp.correlogram.data <- correlog(x=spdf@coords[, 1],
#                                    y=spdf@coords[, 2],
#                                    z=spdf@data[[paste("var", i, sep="")]], 
#                                    increment=correlogram.increment, latlon=TRUE,
#                                    na.rm=TRUE, resamp=50, quiet=TRUE)
#   
#   covariate.correlogram.xintercept[i] <- as.numeric(tmp.correlogram.data$x.intercept)
#   
#   plot(tmp.correlogram.data)
#   mtext(paste("var", i, sep=""))
#   
# }

model.correlogram.data <- data.frame(
  y = correlogram.data$correlation,
  x = correlogram.data$mean.of.class
)

correlogram.polynomial <- lm(y ~ poly(x, 10, raw=TRUE), data=model.correlogram.data)


neighbor.threshold <- 500

tmp.weights <- c()
for (i in 1:nrandom) {
  
  tmp.dist <- spDists(spdf[i, ]@coords, spdf@coords, 
                      longlat=TRUE, segments=FALSE, 
                      diagonal=FALSE)

    
  tmp.neighbors <- tmp.dist < neighbor.threshold
  tmp.treated <- spdf$treat
  
  tmp.newdata <- data.frame(
    x = tmp.dist
  )
  tmp.newdata$y <- predict(correlogram.polynomial, tmp.newdata)
  
  tmp.sum <- sum(abs(tmp.newdata$y) * tmp.neighbors * tmp.treated)
  
  tmp.weights[i] <- tmp.sum
}

spdf$weights <- tmp.weights

spplot(spdf)

# ====================================


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





