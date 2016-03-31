
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
nrandom <- 2000

theta <- 1

# variogram model range (larger = coarser autocorrelation) and psill
var.vrange <- 2000
var.psill <- 1

# numner of covariates
ncovariates <- 1

# define bounding box
minx <- -45
maxx <- 45
miny <- -22.5
maxy <- 22.5

# correlogram increment size
correlogram.increment <- 100

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
proj4string(spdf) <- CRS("+proj=longlat +datum=WGS84")
prj <- proj4string(spdf)
spdf <- spTransform(spdf, CRS(prj))


# using gstat to generate fields with spatial autocorrelation
# source:
#   http://santiago.begueria.es/2010/10/
#     generating-spatially-correlated-random-fields-with-r/

# define the gstat object (spatial model)
var.g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                 model=vgm(psill=var.psill, model="Sph", range=var.vrange),
                 nmax=20)

# make simulations based on the gstat object
# where each simulation will be a variable
var.sim <- predict(var.g.dummy, newdata=spdf, nsim=ncovariates)

# scale 0:1
for (i in names(var.sim)) {
  var.sim@data[[i]] <- var.sim@data[[i]] + abs(min(var.sim@data[[i]]))
  var.sim@data[[i]] <- var.sim@data[[i]] / max(var.sim@data[[i]])
}


# initialize vars for generating treatment var
tmp.var <- rep(0, nrandom)
tmp.random <- runif(nrandom, min=0, max=1)

# add covariates from sim to spatial dataframe
for (i in 1:ncovariates) {
  spdf@data[[paste("var", i, sep="")]] <- 
    var.sim@data[[paste("sim", i, sep="")]]

  tmp.cov <- spdf@data[[paste("var", i, sep="")]]
  tmp.var <- tmp.var + tmp.cov / max(tmp.cov)
}

# incorporate covariate distribution into treatment
tmp.var.avg <- tmp.var / ncovariates
tmp.treat.val <- tmp.random/(1+ncovariates) + tmp.var.avg
treatment.binary <- ifelse(tmp.treat.val > quantile(tmp.treat.val, 0.50), 1, 0)
spdf$treatment.status = treatment.binary


# ---------------------------
# random error

error.vrange <- 1
error.psill <- 0.5
error.ratio <- 1
error.scale <- 0

error.g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                 model=vgm(psill=error.psill, model="Sph", range=error.vrange),
                 nmax=20)

error.sim <- predict(error.g.dummy, newdata=spdf, nsim=1)

error.sim$sim1 <- error.sim$sim1 + abs(min(error.sim$sim1))
error.sim$sim1 <- error.sim$sim1 / max(error.sim$sim1)
error.var.base <- error.sim$sim1


error.tmp <- rep(0, nrandom)
error.random <- runif(nrandom, min=0, max=1)

error.var.modified <- error.random + 
                      error.var.base * error.ratio + 
                      tmp.var.avg * (1 - error.ratio)

error.var <- (error.var.modified / 3) * error.scale * ncovariates

spdf$error <-  error.var

# ---------------------------


# plot spatial autocorrelation of covariates
spplot(spdf, zcol=names(spdf)[names(spdf) != "id"])



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
m1.out <- matchit(treatment.status ~ var1, data=spdf@data,
                  method="nearest", distance="logit", 
                  caliper=1, calclosest=FALSE, calrandom=FALSE)


# ====================================


# get correlogram for PSM distance, x-intercept and plot
correlogram.data <- correlog(x=spdf@coords[, 1],
                             y=spdf@coords[, 2],
                             z=m1.out$distance,
                             increment=correlogram.increment,
                             latlon=TRUE, na.rm=TRUE, resamp=10,
                             quiet=TRUE)

correlogram.xintercept <- as.numeric(correlogram.data$x.intercept)

plot.correlog(correlogram.data)
mtext("PSM pscores")

# covariate.correlogram.xintercept <- c()
# for (i in 1:ncovariates) {
#
#   tmp.correlogram.data <- correlog(x=spdf@coords[, 1],
#                                    y=spdf@coords[, 2],
#                                    z=spdf@data[[paste("var", i, sep="")]],
#                                    increment=correlogram.increment,
#                                    latlon=TRUE, na.rm=TRUE, resamp=50,
#                                    quiet=TRUE)
#
#   covariate.correlogram.xintercept[i] <-
#     as.numeric(tmp.correlogram.data$x.intercept)
#
#   plot(tmp.correlogram.data)
#   mtext(paste("var", i, sep=""))
#
# }

model.correlogram.data <- data.frame(
  y = correlogram.data$correlation,
  x = correlogram.data$mean.of.class
)

correlogram.polynomial <- lm(y ~ poly(x, 10, raw=TRUE),
                             data=model.correlogram.data)



neighbor.threshold <- correlogram.xintercept
if (neighbor.threshold == 0 ) {
  neighbor.threshold = max(correlogram.data$mean.of.class)
}

tmp.spillover.weights <- c()
for (i in 1:nrandom) {

  tmp.dist <- spDists(spdf[i, ]@coords, spdf@coords, longlat=TRUE)


  tmp.neighbors <- tmp.dist < neighbor.threshold & tmp.dist > 0
  tmp.treated <- spdf$treatment.status

  tmp.newdata <- data.frame(
    x = c(tmp.dist)
  )
  tmp.newdata$y <- predict(correlogram.polynomial, tmp.newdata)

  tmp.sum <- sum(abs(tmp.newdata$y) * tmp.neighbors * tmp.treated)

  tmp.spillover.weights[i] <- tmp.sum / sum(tmp.neighbors * tmp.treated)

}

tmp.spillover.t1 <- sum(theta * spdf$treatment.status)

tmp.spillover.t2 <- c()
sum.spillover.weights <- sum(tmp.spillover.weights)
for (i in 1:nrandom) {
  tmp.spillover.t2[i] <- tmp.spillover.weights[i] / sum.spillover.weights
}

# z <- 5000
#
# spdf$treatment.effect <- theta * spdf$treatment.status
# spdf$spillover <- z * tmp.spillover.t1 * tmp.spillover.t2
# spdf$ancillary <- tmp.var.avg
# spdf$error <- 0
# spdf$intercept <- 0
#
# spdf$outcome <- spdf$treatment.effect + spdf$spillover +
#                 spdf$ancillary + spdf$error + spdf$intercept

# spplot(spdf, zcol=names(spdf)[names(spdf) != "id"])

# spplot(spdf, zcol=c('treatment.effect'), main='treatment.effect')
# spplot(spdf, zcol=c('spillover'), main='spillover')
# spplot(spdf, zcol=c('ancillary'), main='ancillary')
#
# spplot(spdf, zcol=c('treatment.effect', 'spillover', 'ancillary'),
#        main='treatment.effect, spillover, ancillary')
#
# spplot(spdf, zcol=c('outcome'), main='outcome')

# spplot(spdf, zcol=c('error'))
# spplot(spdf, zcol=c('intercept'))

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
                                     longlat=TRUE)
  }
}

m1.autocorrelation.match.count <- length(
  m1.match.distances[!is.na(m1.match.distances) &
                     m1.match.distances < correlogram.xintercept])


# -----------------------------------------------------------------------------


# spatial matchit
# spatial.opts <- list(decay.model = "gaussian.semivariance",
#                      threshold = 0.05)

spatial.opts <- list(decay.model = "threshold",
                     threshold = 4000)#correlogram.xintercept)


detach("package:MatchIt", unload=TRUE)
load_all("~/git/matchit/R")
library(MatchIt)

m2.out <- matchit(treatment.status ~ var1, data=spdf,
                  method = "nearest", distance = "logit", 
                  caliper=1, calclosest=FALSE, calrandom=FALSE,
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
                                     longlat=TRUE)
  }
}

m2.autocorrelation.match.count <- length(
  m2.match.distances[!is.na(m2.match.distances) &
                     m2.match.distances < correlogram.xintercept])



m1.out
m2.out
mean(m1.match.distances, na.rm=T)
mean(m2.match.distances, na.rm=T)
hist(m1.match.distances)
hist(m2.match.distances)

# -----------------------------------------------------------------------------

###

# matchit.traditional.data <- match.data(m1.out)
# tmp.traditional.pairs <- rep(0, length(matchit.traditional.data))
# 
# for (i in 1:length(m1.out$match.matrix)) {
#   control.id <- as.numeric(labels(m1.out$match.matrix)[[1]][i])
#   treated.id <- m1.out$match.matrix[i]
# 
#   if (!is.na(treated.id)) {
#     treated.id <- as.numeric(treated.id)
# 
#     pair.id <- factor(paste(as.character(control.id),
#                             as.character(treated.id), sep="_"))
# 
#     tmp.traditional.pairs[control.id] <- pair.id
#     tmp.traditional.pairs[treated.id] <- pair.id
# 
#   }
# }
# 
# matchit.traditional.data$pairs <- tmp.traditional.pairs
# 
# 
# matchit.spatial.data <- match.data(m2.out)
# tmp.spatial.pairs <- rep(0, length(matchit.spatial.data))
# 
# for (i in 1:length(m2.out$match.matrix)) {
#   control.id <- as.numeric(labels(m2.out$match.matrix)[[1]][i])
#   treated.id <- m1.out$match.matrix[i]
# 
#   if (!is.na(treated.id)) {
#     treated.id <- as.numeric(treated.id)
# 
#     pair.id <- factor(paste(as.character(control.id),
#                             as.character(treated.id), sep="_"))
# 
#     tmp.spatial.pairs[control.id] <- pair.id
#     tmp.spatial.pairs[treated.id] <- pair.id
# 
#   }
# }
# 
# matchit.spatial.data$pairs <- tmp.spatial.pairs

###


# traditional.diff.a <- c()
# spatial.diff.a <- c()

traditional.diff.b <- c()
spatial.diff.b <- c()


true.treatment <- c()

z.vals <- c(seq(0, 3000, 1))

# spdf$spillover.var <- spdf$spillover.weights
for (i in 1:length(z.vals)) {

  z <- z.vals[i]
  
  spdf$treatment.effect <- 1 * spdf$treatment.status
  
  spdf$z <- z
  spdf$spillover.t1 <- tmp.spillover.t1

  spdf$spillover <- z * tmp.spillover.t1 * tmp.spillover.t2

  
  spdf$ancillary <- tmp.var
  # spdf$error <- 0
  spdf$intercept <- 0
  
  spdf$outcome <- spdf$treatment.effect + spdf$spillover +
    spdf$ancillary + spdf$error + spdf$intercept

  true.treatment[i] <- theta + (z * tmp.spillover.t1 / nrandom)

  # spplot(spdf, zcol=c('treatment.effect'), main='treatment.effect')
  # spplot(spdf, zcol=c('spillover'), main='spillover')
  # spplot(spdf, zcol=c('ancillary'), main='ancillary')

  # spplot(spdf, zcol=c('treatment.effect', 'spillover', 'ancillary'),
  #        main='treatment.effect, spillover, ancillary')

  # spplot(spdf, zcol=c('outcome'), main='outcome')

  
  
#   test.model.traditional.a <- lm(outcome ~ treatment.status + var1 + var2 + var3,
#                                  data=match.data(m1.out))
#   test.model.spatial.a <- lm(outcome ~ treatment.status + var1 + var2 + var3,
#                              data=match.data(m2.out))
# 
#   test.predict.traditional.a <- predict(test.model.traditional.a, spdf@data)
#   test.predict.spatial.a <- predict(test.model.spatial.a, spdf@data)
# 
#   traditional.diff.a[i] <- mean(abs(test.predict.traditional.a - spdf$outcome))
#   spatial.diff.a[i] <- mean(abs(test.predict.spatial.a - spdf$outcome))
#   
#   if (i == 1) {
#     traditional.coef.a <- summary(test.model.traditional.a)$coef[, 1]
#     spatial.coef.a <- summary(test.model.spatial.a)$coef[, 1]
#     
#   } else {
#     traditional.coef.a <- rbind(traditional.coef.a, 
#                                 summary(test.model.traditional.a)$coef[, 1])
#     spatial.coef.a <- rbind(spatial.coef.a, 
#                             summary(test.model.spatial.a)$coef[, 1])
#   }
  
  
  test.model.traditional.b <- lm(outcome ~ 0 + treatment.status + var1 + var2 + var3,# + spillover.var,
                                 data=match.data(m1.out))
  test.model.spatial.b <- lm(outcome ~ 0 + treatment.status + var1 + var2 + var3,# + spillover.var,
                             data=match.data(m2.out))

  test.predict.traditional.b <- predict(test.model.traditional.b, spdf@data)
  test.predict.spatial.b <- predict(test.model.spatial.b, spdf@data)

  traditional.diff.b[i] <- mean(abs(test.predict.traditional.b - spdf$outcome))
  spatial.diff.b[i] <- mean(abs(test.predict.spatial.b - spdf$outcome))
  
  
  if (i == 1) {
    traditional.coef.b <- summary(test.model.traditional.b)$coef[, 1]
    spatial.coef.b <- summary(test.model.spatial.b)$coef[, 1]
    
  } else {
    traditional.coef.b <- rbind(traditional.coef.b, 
                                summary(test.model.traditional.b)$coef[, 1])
    spatial.coef.b <- rbind(spatial.coef.b, 
                            summary(test.model.spatial.b)$coef[, 1])
  }

  
}


plot(z.vals, true.treatment, main='true')



# plot(z.vals, spatial.diff.a, col="green", type="l", 
#      ylim=c(-10, 10), main='diff a')
# lines(z.vals, traditional.diff.a, col="blue")
# 
# legend(0, 10, 
#        c('traditional', 'spatial'), 
#        lty=c(1, 1), 
#        lwd=c(2.5, 2.5),
#        col=c('blue', 'green')) 



plot(z.vals, spatial.diff.b, col="green", type="l", 
     ylim=c(min(spatial.diff.b), max(traditional.diff.b)), main='diff b')
lines(z.vals, traditional.diff.b, col="blue")

legend(0, 20, 
       c('traditional', 'spatial'), 
       lty=c(1, 1), 
       lwd=c(2.5, 2.5),
       col=c('blue', 'green')) 



# plot(x=z.vals, type="n", ylim=c(-10, 1000), main='Coef a', ylab="", xlab="z")
# 
# lines(z.vals, traditional.coef.a[, '(Intercept)'], col="black", lty=2, lwd=1)
# lines(z.vals, spatial.coef.a[, '(Intercept)'], col="black", lty=1, lwd=1)
# 
# lines(z.vals, traditional.coef.a[, 'treatment.status'], col="red", lty=2, lwd=1)
# lines(z.vals, spatial.coef.a[, 'treatment.status'], col="red", lty=1, lwd=1)
# 
# lines(z.vals, traditional.coef.a[, 'var1'], col="green", lty=2, lwd=1)
# lines(z.vals, spatial.coef.a[, 'var1'], col="green", lty=1, lwd=1)
# 
# lines(z.vals, traditional.coef.a[, 'var2'], col="blue", lty=2, lwd=1)
# lines(z.vals, spatial.coef.a[, 'var2'], col="blue", lty=1, lwd=1)
# 
# lines(z.vals, traditional.coef.a[, 'var3'], col="purple", lty=2, lwd=1)
# lines(z.vals, spatial.coef.a[, 'var3'], col="purple", lty=1, lwd=1)
# 
# lines(z.vals, traditional.coef.a[, 'spillover.var'], col="brown", lty=2, lwd=1)
# lines(z.vals, spatial.coef.a[, 'spillover.var'], col="brown", lty=1, lwd=1)

# legend(0, 25, 
#        c('traditional', 'spatial'), 
#        lty=c(2, 1), 
#        lwd=c(2.5)) 
# 
# legend(0, 20, 
#        c('intercept', 'treatment.status', 'var1', 'var2', 'var3', 'spillover.var'), 
#        lty=c(1), 
#        lwd=c(2.5), 
#        col=c("black", "red", "green", "blue", "purple", "brown")) 



plot(x=z.vals, type="n", ylim=c(-100, 1000), xlim=c(min(z.vals), max(z.vals)), main='Coef b', ylab="", xlab="z")

lines(z.vals, traditional.coef.b[, 'treatment.status'], col="red", lty=2, lwd=1)
lines(z.vals, spatial.coef.b[, 'treatment.status'], col="red", lty=1, lwd=1)

lines(z.vals, traditional.coef.b[, 'var1'], col="green", lty=2, lwd=1)
lines(z.vals, spatial.coef.b[, 'var1'], col="green", lty=1, lwd=1)

lines(z.vals, traditional.coef.b[, 'var2'], col="blue", lty=2, lwd=1)
lines(z.vals, spatial.coef.b[, 'var2'], col="blue", lty=1, lwd=1)

lines(z.vals, traditional.coef.b[, 'var3'], col="purple", lty=2, lwd=1)
lines(z.vals, spatial.coef.b[, 'var3'], col="purple", lty=1, lwd=1)

lines(z.vals, traditional.coef.b[, 'spillover.var'], col="brown", lty=2, lwd=1)
lines(z.vals, spatial.coef.b[, 'spillover.var'], col="brown", lty=1, lwd=1)

legend(0, 500, 
       c('traditional', 'spatial'), 
       lty=c(2, 1), 
       lwd=c(2.5)) 

legend(0, 350, 
       c('treatment.status', 'var1', 'var2', 'var3', 'spillover.var'), 
       lty=c(1), 
       lwd=c(2.5), 
       col=c("red", "green", "blue", "purple", "brown")) 





