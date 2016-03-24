
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

theta <- 1

# variogram model range (larger = coarser autocorrelation)
vrange <- 3500

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
g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1,
                 model=vgm(psill=0.5, model="Sph", range=vrange),
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
tmp.var.avg <- tmp.var / ncovariates
tmp.treat.val <- tmp.random/(1+ncovariates) + tmp.var.avg
treatment.binary <- ifelse(tmp.treat.val > quantile(tmp.treat.val, 0.50), 1, 0)
spdf$treatment.status = treatment.binary



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
m1.out <- matchit(treatment.status ~ var1 + var2 + var3, data=spdf@data,
                  method="nearest", distance="logit", caliper=0.25)


# ====================================


# get correlogram for PSM distance, x-intercept and plot
correlogram.data <- correlog(x=spdf@coords[, 1],
                                y=spdf@coords[, 2],
                                z=m1.out$distance,
                                increment=correlogram.increment,
                                latlon=TRUE, na.rm=TRUE, resamp=50,
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

tmp.spillover.weights <- c()
for (i in 1:nrandom) {

  tmp.dist <- spDists(spdf[i, ]@coords, spdf@coords,
                      longlat=TRUE, segments=FALSE,
                      diagonal=FALSE)


  tmp.neighbors <- tmp.dist < neighbor.threshold & tmp.dist > 0
  tmp.treated <- spdf$treatment.status

  tmp.newdata <- data.frame(
    x = c(tmp.dist)
  )
  tmp.newdata$y <- predict(correlogram.polynomial, tmp.newdata)

  tmp.sum <- sum(abs(tmp.newdata$y) * tmp.neighbors * tmp.treated)

  tmp.spillover.weights[i] <- tmp.sum / sum(tmp.neighbors * tmp.treated)

}
spdf$spillover.weights <- tmp.spillover.weights

tmp.spillover.t1 <- sum(theta * spdf$treatment.status)

tmp.spillover.t2 <- c()
for (i in 1:nrandom) {
  tmp.spillover.t2[i] <-
    spdf$spillover.weights[i] / sum(spdf$spillover.weights)
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
                                     longlat=FALSE, segments=FALSE,
                                     diagonal=FALSE)
  }
}

m1.autocorrelation.match.count <- length(
  m1.match.distances[!is.na(m1.match.distances) &
                     m1.match.distances < (correlogram.xintercept/110)])


# -----------------------------------------------------------------------------


# spatial matchit
# spatial.opts <- list(decay.model = "gaussian.semivariance",
#                      threshold = 0.05)

spatial.opts <- list(decay.model = "spherical",
                     threshold = (correlogram.xintercept/110))

m2.out <- matchit(treatment.status ~ var1 + var2 + var3, data=spdf,
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

m2.autocorrelation.match.count <- length(
  m2.match.distances[!is.na(m2.match.distances) &
                     m2.match.distances < (correlogram.xintercept/110)])


# -----------------------------------------------------------------------------


#
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



traditional.diff <- c()
spatial.diff <- c()

traditional.coef <- c()
spatial.coef <- c()

true.treatment <- c()

z.vals <- c(seq(0,500,1))

for (i in 1:length(z.vals)) {

  z <- z.vals[i]

  spdf$treatment.effect <- 1 * spdf$treatment.status
  spdf$spillover <- z * tmp.spillover.t1 * tmp.spillover.t2
  spdf$ancillary <- tmp.var.avg
  spdf$error <- 0
  spdf$intercept <- 0

  spdf$outcome <- spdf$treatment.effect + spdf$spillover +
    spdf$ancillary + spdf$error + spdf$intercept

  true.treatment[i] <- theta + mean(spdf$spillover)

#   spplot(spdf, zcol=c('treatment.effect'), main='treatment.effect')
#   spplot(spdf, zcol=c('spillover'), main='spillover')
#   spplot(spdf, zcol=c('ancillary'), main='ancillary')
#
#   spplot(spdf, zcol=c('treatment.effect', 'spillover', 'ancillary'),
#          main='treatment.effect, spillover, ancillary')
#
#   spplot(spdf, zcol=c('outcome'), main='outcome')

  test.model.traditional <- lm(outcome ~ treatment.status + var1 + var2 + var3,
                               data=match.data(m1.out))
  test.model.spatial <- lm(outcome ~ treatment.status + var1 + var2 + var3,
                           data=match.data(m2.out))

  test.predict.traditional <- predict(test.model.traditional, spdf@data)
  test.predict.spatial <- predict(test.model.spatial, spdf@data)

  traditional.diff[i] <- mean(test.predict.traditional - spdf$outcome)
  spatial.diff[i] <- mean(test.predict.spatial - spdf$outcome)

  traditional.coef[i] <- summary(test.model.traditional)$coef[2,1]
  spatial.coef[i] <- summary(test.model.spatial)$coef[2,1]

}

plot(z.vals, spatial.diff, col="green", ylim=c(-10, 10))
lines(z.vals, traditional.diff, col="blue")

plot(z.vals, spatial.coef, col="green", ylim=c(-10, 10))
lines(z.vals, traditional.coef, col="blue")

plot(z.vals, true.treatment)











