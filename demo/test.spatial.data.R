
library(devtools)

detach("package:MatchIt", unload=TRUE)
load_all("~/Desktop/Github/MatchIt/R")
# library(devtools)
# install_github("itpir/matchit")
library(MatchIt)

library(sp)
library(ncf)
library(gstat)

rm(list = ls())
# -----------------------------------------------------------------------------
# options

# set dataframe size (number of points)
nrandom <- 3000
control.ratio <- 0.90

theta <- 1

# variogram model range (larger = coarser autocorrelation) and psill
var.vrange <- 2000
var.psill <- 1.0

# numner of covariates
ncovariates <- 1

# define bounding box
minx <- -45
maxx <- 45
miny <- -22.5
maxy <- 22.5

# correlogram increment size
correlogram.increment <- 100

# set = 1 for balanced ratio based on number of covariates
# set < 1 for more random
# set > 1 for greater spatial autocorrelation based on covariates
treatment.autocorrelation <- 0.5

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
tmp.treat.val <- tmp.random/((1+ncovariates) * treatment.autocorrelation) + tmp.var.avg
treatment.binary <- ifelse(tmp.treat.val > quantile(tmp.treat.val, control.ratio), 1, 0)
spdf$treatment.status = treatment.binary


# ---------------------------
# random error

error.vrange <- 1
error.psill <- 0.5
error.ratio <- 1
error.scale <- 5.0

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
spplot(spdf, zcol=names(spdf)[names(spdf) %in% c("treatment.status", "var1")])


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
                  caliper=0.25, calclosest=FALSE, calrandom=FALSE)


# ====================================


# get correlogram for PSM distance, x-intercept and plot
correlogram.data <- correlog(x=spdf@coords[, 1],
                             y=spdf@coords[, 2],
                             z=m1.out$distance,
                             increment=correlogram.increment,
                             latlon=TRUE, na.rm=TRUE, resamp=10,
                             quiet=FALSE)

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
tmp.spillover.weights[is.na(tmp.spillover.weights)] <- 0

tmp.spillover.t1 <- sum(theta * spdf$treatment.status)

tmp.spillover.t2 <- c()
sum.spillover.weights <- sum(tmp.spillover.weights, na.rm=TRUE)
for (i in 1:nrandom) {
  tmp.spillover.t2[i] <- tmp.spillover.weights[i] / sum.spillover.weights
  if(is.na(tmp.spillover.t2[i])){
    tmp.spillover.t2[i] = 0
    }
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

spatial.opts <- list(decay.model = "spherical",
                     threshold = neighbor.threshold)


m2.out <- matchit(treatment.status ~ var1, data=spdf,
                  method = "nearest", distance = "logit", 
                  caliper=0.25, calclosest=FALSE, calrandom=FALSE,
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


# -------------------------------------


# traditional.diff.b <- c()
# spatial.diff.b <- c()


true.treatment <- c()

z.vals <- c(seq(0, 25, 0.1))

spdf$spillover.var <- tmp.spillover.weights
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

  # true.treatment[i] <- theta + (z * tmp.spillover.t1 / nrandom)
  true.treatment[i] <- sum(spdf$treatment.effect + spdf$spillover) / nrandom
  

  test.model.traditional.b.noSpill <- lm(outcome ~ 0 + treatment.status + var1,
                                       data=match.data(m1.out)) 
  test.model.spatial.b.noSpill <- lm(outcome ~ 0 + treatment.status + var1,
                                       data=match.data(m2.out)) 
  
  
  test.model.traditional.b.spill <- lm(outcome ~ 0 + treatment.status + var1 + spillover.var,
                                 data=match.data(m1.out))
  test.model.spatial.b.spill <- lm(outcome ~ 0 + treatment.status + var1 + spillover.var,
                             data=match.data(m2.out))

  test.predict.traditional.b.spill <- predict(test.model.traditional.b, spdf@data)
  test.predict.spatial.b.spill <- predict(test.model.spatial.b, spdf@data)

  test.predict.traditional.b.noSpill <- predict(test.model.traditional.b.noSpill, spdf@data)
  test.predict.spatial.b.noSpill <- predict(test.model.spatial.b.noSpill, spdf@data)
  
  # traditional.diff.b[i] <- mean(abs(test.predict.traditional.b - spdf$outcome))
  # spatial.diff.b[i] <- mean(abs(test.predict.spatial.b - spdf$outcome))

  
  if (i == 1) {
    traditional.coef.b.spill <- summary(test.model.traditional.b.spill)$coef[, 1]
    spatial.coef.b.spill <- summary(test.model.spatial.b.spill)$coef[, 1]
    
    traditional.coef.b.noSpill <- summary(test.model.traditional.b.noSpill)$coef[, 1]
    spatial.coef.b.noSpill <- summary(test.model.spatial.b.noSpill)$coef[, 1]
    
  } else {
    traditional.coef.b.spill <- rbind(traditional.coef.b.spill, 
                                summary(test.model.traditional.b.spill)$coef[, 1])
    spatial.coef.b.spill <- rbind(spatial.coef.b.spill, 
                            summary(test.model.spatial.b.spill)$coef[, 1])
    
    traditional.coef.b.noSpill <- rbind(traditional.coef.b.noSpill, 
                                      summary(test.model.traditional.b.noSpill)$coef[, 1])
    spatial.coef.b.noSpill <- rbind(spatial.coef.b.noSpill, 
                                  summary(test.model.spatial.b.noSpill)$coef[, 1])
  }

  
}


plot(z.vals, true.treatment, main='True Average Treatment Effect\n(Including Spillovers)')



max.y.val <- max(traditional.coef.b.spill, spatial.coef.b.spill) * 2
min.y.val <- min(traditional.coef.b.spill, spatial.coef.b.spill)

plot(x=z.vals, type="n", main='Coeffients \n(Treatment Lag Model)', ylab="", xlab="z", 
     ylim=c(min.y.val, max.y.val), 
     xlim=c(min(z.vals), max(z.vals)))

lines(z.vals, traditional.coef.b.spill[, 'treatment.status'] + traditional.coef.b.spill[, 'spillover.var'], 
      col="red", lty=2, lwd=1)
lines(z.vals, spatial.coef.b.spill[, 'treatment.status'] + spatial.coef.b.spill[, 'spillover.var'], 
      col="red", lty=1, lwd=1)

lines(z.vals, abs(traditional.coef.b.spill[, 'var1']), 
      col="green", lty=2, lwd=1)
lines(z.vals, abs(spatial.coef.b.spill[, 'var1']), 
      col="green", lty=1, lwd=1)

lines(z.vals, traditional.coef.b.spill[, 'spillover.var'], 
      col="blue", lty=2, lwd=1)
 lines(z.vals, spatial.coef.b.spill[, 'spillover.var'], 
       col="blue", lty=1, lwd=1)
 
 lines(z.vals, traditional.coef.b.spill[, 'treatment.status'], 
       col="orange", lty=2, lwd=1)
 lines(z.vals, spatial.coef.b.spill[, 'treatment.status'], 
       col="orange", lty=1, lwd=1)


legend(min(z.vals), max.y.val - 0.35*max.y.val, 
       c('traditional', 'spatial'), 
       lty=c(2, 1), lwd=c(2.5)) 

legend(min(z.vals), max.y.val, 
       c('trt + spill', 'var1', 'spill', 'trt'), 
       lty=c(1), lwd=c(2.5), 
       col=c("red", "green", "blue", "orange")) 

trad.treatment.spill <- vector()
spatial.treatment.spill <- vector()

for(i in 1:length(traditional.coef.b.spill[,'treatment.status'])){
  trad.treatment.spill[i] <- sum((spdf@data$treatment.status * traditional.coef.b.spill[, 'treatment.status'][i]) + 
    (spdf@data$spillover.var * traditional.coef.b.spill[, 'spillover.var'][i])) / nrandom
  spatial.treatment.spill[i] <- sum((spdf@data$treatment.status * spatial.coef.b.spill[, 'treatment.status'][i]) + 
    (spdf@data$spillover.var * spatial.coef.b.spill[, 'spillover.var'][i])) / nrandom
}


tmp.traditional.diff <- (true.treatment - trad.treatment.spill)
tmp.spatial.diff <- (true.treatment - spatial.treatment.spill)

  
plot(x=z.vals, type="n", main='Estimated vs. True \n Averge Treatment Effect \n (Model Including Spillover)', ylab="", xlab="z", 
     ylim=c(min(tmp.traditional.diff, tmp.spatial.diff), 
            max(tmp.traditional.diff, tmp.spatial.diff)),
     xlim=c(min(z.vals), max(z.vals)))

lines(z.vals, tmp.traditional.diff, col="green")
lines(z.vals, tmp.spatial.diff, col="blue")

legend(min(z.vals), max(tmp.traditional.diff, tmp.spatial.diff), 
       c('spatial', 'traditional'), 
       lty=c(1), 
       lwd=c(2.5),
       col=c('blue', 'green')) 




#----------------
#No Spillover

max.y.val <- max(traditional.coef.b.noSpill, spatial.coef.b.noSpill)
min.y.val <- min(traditional.coef.b.noSpill, spatial.coef.b.noSpill)

plot(x=z.vals, type="n", main='Coeffients \n(No Treatment Lag)', ylab="", xlab="z", 
     ylim=c(min.y.val, max.y.val), 
     xlim=c(min(z.vals), max(z.vals)))

lines(z.vals, traditional.coef.b.noSpill[, 'treatment.status'], 
      col="red", lty=2, lwd=1)
lines(z.vals, spatial.coef.b.noSpill[, 'treatment.status'], 
      col="red", lty=1, lwd=1)

lines(z.vals, abs(traditional.coef.b.noSpill[, 'var1']), 
      col="green", lty=2, lwd=1)
lines(z.vals, abs(spatial.coef.b.noSpill[, 'var1']), 
      col="green", lty=1, lwd=1)


legend(min(z.vals), max.y.val - 0.35*max.y.val, 
       c('traditional', 'spatial'), 
       lty=c(2, 1), lwd=c(2.5)) 

legend(min(z.vals), max.y.val, 
       c('treatment.status', 'var1'), 
       lty=c(1), lwd=c(2.5), 
       col=c("red", "green")) 

trad.treatment.noSpill <- vector()
spatial.treatment.noSpill <- vector()

for(i in 1:length(traditional.coef.b[,'treatment.status'])){
  trad.treatment.noSpill[i] <- sum(spdf@data$treatment.status * 
                                     traditional.coef.b.noSpill[, 'treatment.status'][i]) / nrandom
  spatial.treatment.noSpill[i] <- sum(spdf@data$treatment.status * 
                                         spatial.coef.b.noSpill[, 'treatment.status'][i]) / nrandom
}


tmp.traditional.diff.noSpill <- (true.treatment - trad.treatment.noSpill)
tmp.spatial.diff.noSpill <- (true.treatment - spatial.treatment.noSpill)


plot(x=z.vals, type="n", main='Estimated vs. True \n Averge Treatment Effect \n (Model Without Spillover)', ylab="", xlab="z", 
     ylim=c(min(tmp.traditional.diff.noSpill, tmp.spatial.diff.noSpill), 
            max(tmp.traditional.diff.noSpill, tmp.spatial.diff.noSpill)),
     xlim=c(min(z.vals), max(z.vals)))

lines(z.vals, tmp.traditional.diff.noSpill, col="green")
lines(z.vals, tmp.spatial.diff.noSpill, col="blue")

legend(min(z.vals), max(tmp.traditional.diff.noSpill, tmp.spatial.diff.noSpill), 
       c('spatial', 'traditional'), 
       lty=c(1), 
       lwd=c(2.5),
       col=c('blue', 'green')) 







