

# -----------------------------------------------------------------------------
#Spillover Correlograms
# -----------------------------------------------------------------------------


m1.out <- matchit(treatment.status ~ modelVar, data=spdf@data,
                  method="nearest", distance="logit", 
                  caliper=caliper, calclosest=FALSE, calrandom=FALSE)

#Calculate the variogram of the PSM scores for matched pairs



spdf@data$m1.pscore <- m1.out$distance

spdf_pairs <- spdf[spdf@data$id %in% match.data(m1.out)[[1]],]

spdf_pairs <- spdf


#Propensity Correlogram
correlog.pscore.spillover <- correlog(x=spdf@coords[, 1],
                                      y=spdf@coords[, 2],
                                      z=spdf$m1.pscore,
                                      increment=500,
                                      latlon=TRUE, na.rm=TRUE, resamp=2,
                                      quiet=FALSE)

pscore.spillover.model.dta <- data.frame(
  mean.of.class = c(correlog.pscore.spillover$mean.of.class),
  correlation = c(correlog.pscore.spillover$correlation)
)

correlog.m1.polynomial <- lm(correlation ~ poly(mean.of.class, degree=5, raw=TRUE),
                             data=pscore.spillover.model.dta)
correlog.pscore.spillover$est_gamma <- predict(correlog.m1.polynomial, 
                                               pscore.spillover.model.dta)


#True Spillover Correlogram
correlog.spillover.true <- correlog(x=spdf@coords[, 1],
                                    y=spdf@coords[, 2],
                                    z=spdf$trueSpill,
                                    increment=500,
                                    latlon=TRUE, na.rm=TRUE, resamp=2,
                                    quiet=FALSE)

correlog.polynomial.true.model.dta <- data.frame(
  mean.of.class = c(correlog.spillover.true$mean.of.class),
  correlation = c(correlog.spillover.true$correlation)
)

correlog.polynomial.true <- lm(correlation ~ poly(mean.of.class, degree=5, raw=TRUE),
                               data=correlog.polynomial.true.model.dta)

correlog.spillover.true$est_gamma <- predict(correlog.polynomial.true, 
                                             correlog.polynomial.true.model.dta)

#Parameterized Correlogram

correlog.parameterized <- correlog(x=trt.spillover.sim@coords[, 1],
                                   y=trt.spillover.sim@coords[, 2],
                                   z=trt.spillover.sim$sim1,
                                   increment=500,
                                   latlon=TRUE, na.rm=TRUE, resamp=2,
                                   quiet=FALSE)

# -----------------------------------------------------------------------------
#Figures
# -----------------------------------------------------------------------------
# plot spatial autocorrelation of covariates (true and with error)
#spplot(spdf, zcol=names(spdf_pairs)[names(spdf_pairs) != "id"])
spplot(spdf, zcol=names(spdf_pairs)[names(spdf_pairs) == "trueSpill"])
spplot(spdf, zcol=names(spdf_pairs)[names(spdf_pairs) == "treatment.status"])

#vgm1 <- variogram(trueVar~1, spdf, covariogram=TRUE)
#vgm2 <- variogram(modelVar~1, spdf, covariogram=TRUE)

#plot(vgm1, main="True Ancillary Spatial Decay")
#plot(vgm2, main="Model Ancillary Spatial Decay")



plot(ylim=c(-1.0,1.0), correlog.parameterized$mean.of.class, correlog.parameterized$correlation, col="red", main="Spillover Correlograms")
#lines(vgm.spillover$distkm, vgm.spillover$fit, col="red")

lines(correlog.spillover.true$mean.of.class, 
      correlog.spillover.true$est_gamma, col="green")
points(correlog.spillover.true$mean.of.class, 
       correlog.spillover.true$correlation, col="green")

lines(correlog.pscore.spillover$mean.of.class, 
      correlog.pscore.spillover$est_gamma, col="blue")
points(correlog.pscore.spillover$mean.of.class, 
       correlog.pscore.spillover$correlation, col="blue")

legend("topright", legend=c("Parameterized Spillover","True Spillover", "PSM-approximated"), pch=c(pch = 1, pch=1, pch=1),
       col=c(col="red", col="green", col="blue"), title = "Legend")
