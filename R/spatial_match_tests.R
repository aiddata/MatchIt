# 
# 
# 
# 
# # traditional, non-spatial matchit
# 
# 
# 
# # ====================================
# 
# 
# # get correlogram for PSM distance, x-intercept and plot
# correlogram.data <- correlog(x=spdf@coords[, 1],
#                              y=spdf@coords[, 2],
#                              z=m1.out$distance,
#                              increment=correlogram.increment,
#                              latlon=TRUE, na.rm=TRUE, resamp=10,
#                              quiet=FALSE)
# 
# correlogram.xintercept <- as.numeric(correlogram.data$x.intercept)
# 
# plot.correlog(correlogram.data)
# mtext("PSM pscores (with error)")
# 
# 
# model.correlogram.data <- data.frame(
#   y = correlogram.data$correlation,
#   x = correlogram.data$mean.of.class
# )
# 
# correlogram.polynomial <- lm(y ~ poly(x, 10, raw=TRUE),
#                              data=model.correlogram.data)
# 
# 
# #Second model using the model without error for data generation.
# 
# m1.out.true <- matchit(treatment.status ~ var1.true, data=spdf@data,
#                        method="nearest", distance="logit", 
#                        caliper=0.25, calclosest=FALSE, calrandom=FALSE)
# 
# correlogram.data.true <- correlog(x=spdf@coords[, 1],
#                                   y=spdf@coords[, 2],
#                                   z=m1.out.true$distance,
#                                   increment=correlogram.increment,
#                                   latlon=TRUE, na.rm=TRUE, resamp=10,
#                                   quiet=FALSE)
# 
# correlogram.xintercept.true  <- as.numeric(correlogram.data.true $x.intercept)
# 
# plot.correlog(correlogram.data.true )
# mtext("PSM pscores (true)")
# 
# 
# model.correlogram.data.true  <- data.frame(
#   y = correlogram.data.true $correlation,
#   x = correlogram.data.true $mean.of.class
# )
# 
# correlogram.polynomial.true <- lm(y ~ poly(x, 10, raw=TRUE),
#                                   data=model.correlogram.data.true )
# 
# 
# 
# 
# 
# 
# m1.match.distances <- c()
# 
# for (i in 1:length(m1.out$match.matrix)) {
#   control.id <- as.numeric(labels(m1.out$match.matrix)[[1]][i])
#   treated.id <- m1.out$match.matrix[i]
#   
#   if (!is.na(treated.id)) {
#     treated.id <- as.numeric(treated.id)
#     
#     control.coords <- spdf[control.id, ]@coords
#     treated.coords <- spdf[treated.id, ]@coords
#     
#     m1.match.distances[i] <- spDists(control.coords, treated.coords,
#                                      longlat=TRUE)
#   }
# }
# 
# m1.autocorrelation.match.count <- length(
#   m1.match.distances[!is.na(m1.match.distances) &
#                        m1.match.distances < correlogram.xintercept])
# 
# 
# # -----------------------------------------------------------------------------
# 
# 
# # spatial matchit
# # spatial.opts <- list(decay.model = "gaussian.semivariance",
# #                      threshold = 0.05)
# 
# spatial.opts <- list(decay.model = "threshold",
#                      threshold = neighbor.threshold)
# 
# 
# m2.out <- matchit(treatment.status ~ var1, data=spdf,
#                   method = "nearest", distance = "logit", 
#                   caliper=0, calclosest=FALSE, calrandom=FALSE,
#                   spatial.options=spatial.opts)
# 
# # -------------------------------------
# 
# m2.match.distances <- c()
# 
# for (i in 1:length(m2.out$match.matrix)) {
#   control.id <- as.numeric(labels(m2.out$match.matrix)[[1]][i])
#   treated.id <- m2.out$match.matrix[i]
#   
#   if (!is.na(treated.id)) {
#     treated.id <- as.numeric(treated.id)
#     
#     control.coords <- spdf[control.id, ]@coords
#     treated.coords <- spdf[treated.id, ]@coords
#     
#     m2.match.distances[i] <- spDists(control.coords, treated.coords,
#                                      longlat=TRUE)
#   }
# }
# 
# m2.autocorrelation.match.count <- length(
#   m2.match.distances[!is.na(m2.match.distances) &
#                        m2.match.distances < correlogram.xintercept])
# 
# 
# 
# m1.out
# m2.out
# mean(m1.match.distances, na.rm=T)
# mean(m2.match.distances, na.rm=T)
# hist(m1.match.distances)
# hist(m2.match.distances)
# 
# 
# 
# 
# true.treatment <- c()
# 
# z.vals <- c(seq(0, 50, 0.1))
# 
# spdf$spillover.var <- tmp.spillover.weights
# for (i in 1:length(z.vals)) {
#   
#   z <- z.vals[i]
#   
#   spdf$treatment.effect <- 1 * spdf$treatment.status
#   
#   spdf$z <- z
#   spdf$spillover.t1 <- tmp.spillover.t1
#   
#   spdf$spillover <- z * tmp.spillover.t1 * tmp.spillover.t2
#   
#   
#   spdf$ancillary <- tmp.var
#   
#   spdf$intercept <- 0
#   
#   spdf$outcome <- spdf$treatment.effect + spdf$spillover +
#     spdf$ancillary + spdf$intercept + spdf$error.outcome
#   
#   # true.treatment[i] <- theta + (z * tmp.spillover.t1 / nrandom)
#   true.treatment[i] <- sum(spdf$treatment.effect + spdf$spillover) / nrandom
#   
#   
#   test.model.traditional.b.noSpill <- lm(outcome ~ 0 + treatment.status + var1,
#                                          data=match.data(m1.out)) 
#   test.model.spatial.b.noSpill <- lm(outcome ~ 0 + treatment.status + var1,
#                                      data=match.data(m2.out)) 
#   
#   
#   test.model.traditional.b.spill <- lm(outcome ~ 0 + treatment.status + var1 + spillover.var,
#                                        data=match.data(m1.out))
#   test.model.spatial.b.spill <- lm(outcome ~ 0 + treatment.status + var1 + spillover.var,
#                                    data=match.data(m2.out))
#   
#   test.predict.traditional.b.spill <- predict(test.model.traditional.b.spill, spdf@data)
#   test.predict.spatial.b.spill <- predict(test.model.spatial.b.spill, spdf@data)
#   
#   test.predict.traditional.b.noSpill <- predict(test.model.traditional.b.noSpill, spdf@data)
#   test.predict.spatial.b.noSpill <- predict(test.model.spatial.b.noSpill, spdf@data)
#   
#   # traditional.diff.b[i] <- mean(abs(test.predict.traditional.b - spdf$outcome))
#   # spatial.diff.b[i] <- mean(abs(test.predict.spatial.b - spdf$outcome))
#   
#   
#   if (i == 1) {
#     traditional.coef.b.spill <- summary(test.model.traditional.b.spill)$coef[, 1]
#     spatial.coef.b.spill <- summary(test.model.spatial.b.spill)$coef[, 1]
#     
#     traditional.coef.b.noSpill <- summary(test.model.traditional.b.noSpill)$coef[, 1]
#     spatial.coef.b.noSpill <- summary(test.model.spatial.b.noSpill)$coef[, 1]
#     
#   } else {
#     traditional.coef.b.spill <- rbind(traditional.coef.b.spill, 
#                                       summary(test.model.traditional.b.spill)$coef[, 1])
#     spatial.coef.b.spill <- rbind(spatial.coef.b.spill, 
#                                   summary(test.model.spatial.b.spill)$coef[, 1])
#     
#     traditional.coef.b.noSpill <- rbind(traditional.coef.b.noSpill, 
#                                         summary(test.model.traditional.b.noSpill)$coef[, 1])
#     spatial.coef.b.noSpill <- rbind(spatial.coef.b.noSpill, 
#                                     summary(test.model.spatial.b.noSpill)$coef[, 1])
#   }
#   
#   
# }
# 
# 
# plot(z.vals, true.treatment, main='True Average Treatment Effect\n(Including Spillovers)')
# 
# 
# 
# max.y.val <- max(traditional.coef.b.spill, spatial.coef.b.spill) * 2
# min.y.val <- min(traditional.coef.b.spill, spatial.coef.b.spill)
# 
# plot(x=z.vals, type="n", main='Coeffients \n(Treatment Lag Model)', ylab="", xlab="z", 
#      ylim=c(min.y.val, max.y.val), 
#      xlim=c(min(z.vals), max(z.vals)))
# 
# lines(z.vals, traditional.coef.b.spill[, 'treatment.status'] + traditional.coef.b.spill[, 'spillover.var'], 
#       col="red", lty=2, lwd=1)
# lines(z.vals, spatial.coef.b.spill[, 'treatment.status'] + spatial.coef.b.spill[, 'spillover.var'], 
#       col="red", lty=1, lwd=1)
# 
# lines(z.vals, abs(traditional.coef.b.spill[, 'var1']), 
#       col="green", lty=2, lwd=1)
# lines(z.vals, abs(spatial.coef.b.spill[, 'var1']), 
#       col="green", lty=1, lwd=1)
# 
# lines(z.vals, traditional.coef.b.spill[, 'spillover.var'], 
#       col="blue", lty=2, lwd=1)
# lines(z.vals, spatial.coef.b.spill[, 'spillover.var'], 
#       col="blue", lty=1, lwd=1)
# 
# lines(z.vals, traditional.coef.b.spill[, 'treatment.status'], 
#       col="orange", lty=2, lwd=1)
# lines(z.vals, spatial.coef.b.spill[, 'treatment.status'], 
#       col="orange", lty=1, lwd=1)
# 
# 
# legend(min(z.vals), max.y.val - 0.35*max.y.val, 
#        c('traditional', 'spatial'), 
#        lty=c(2, 1), lwd=c(2.5)) 
# 
# legend(min(z.vals), max.y.val, 
#        c('trt + spill', 'var1', 'spill', 'trt'), 
#        lty=c(1), lwd=c(2.5), 
#        col=c("red", "green", "blue", "orange")) 
# 
# trad.treatment.spill <- vector()
# spatial.treatment.spill <- vector()
# 
# for(j in 1:length(traditional.coef.b.spill[,'treatment.status'])){
#   trad.treatment.spill[j] <- sum((spdf@data$treatment.status * traditional.coef.b.spill[, 'treatment.status'][j]) + 
#                                    (spdf@data$spillover.var * traditional.coef.b.spill[, 'spillover.var'][j])) / nrandom
#   spatial.treatment.spill[j] <- sum((spdf@data$treatment.status * spatial.coef.b.spill[, 'treatment.status'][j]) + 
#                                       (spdf@data$spillover.var * spatial.coef.b.spill[, 'spillover.var'][j])) / nrandom
# }
# 
# 
# tmp.traditional.diff <- (true.treatment - trad.treatment.spill)
# tmp.spatial.diff <- (true.treatment - spatial.treatment.spill)
# 
# 
# plot(x=z.vals, type="n", main='Estimated vs. True \n Averge Treatment Effect \n (Model Including Spillover)', ylab="", xlab="z", 
#      ylim=c(min(tmp.traditional.diff, tmp.spatial.diff), 
#             max(tmp.traditional.diff, tmp.spatial.diff)),
#      xlim=c(min(z.vals), max(z.vals)))
# 
# lines(z.vals, tmp.traditional.diff, col="green")
# lines(z.vals, tmp.spatial.diff, col="blue")
# 
# legend(min(z.vals), max(tmp.traditional.diff, tmp.spatial.diff), 
#        c('spatial', 'traditional'), 
#        lty=c(1), 
#        lwd=c(2.5),
#        col=c('blue', 'green')) 
# 
# 
# 
# 
# #----------------
# #No Spillover
# 
# max.y.val <- max(traditional.coef.b.noSpill, spatial.coef.b.noSpill)
# min.y.val <- min(traditional.coef.b.noSpill, spatial.coef.b.noSpill)
# 
# plot(x=z.vals, type="n", main='Coeffients \n(No Treatment Lag)', ylab="", xlab="z", 
#      ylim=c(min.y.val, max.y.val), 
#      xlim=c(min(z.vals), max(z.vals)))
# 
# lines(z.vals, traditional.coef.b.noSpill[, 'treatment.status'], 
#       col="red", lty=2, lwd=1)
# lines(z.vals, spatial.coef.b.noSpill[, 'treatment.status'], 
#       col="red", lty=1, lwd=1)
# 
# lines(z.vals, abs(traditional.coef.b.noSpill[, 'var1']), 
#       col="green", lty=2, lwd=1)
# lines(z.vals, abs(spatial.coef.b.noSpill[, 'var1']), 
#       col="green", lty=1, lwd=1)
# 
# 
# legend(min(z.vals), max.y.val - 0.35*max.y.val, 
#        c('traditional', 'spatial'), 
#        lty=c(2, 1), lwd=c(2.5)) 
# 
# legend(min(z.vals), max.y.val, 
#        c('treatment.status', 'var1'), 
#        lty=c(1), lwd=c(2.5), 
#        col=c("red", "green")) 
# 
# trad.treatment.noSpill <- vector()
# spatial.treatment.noSpill <- vector()
# 
# for(j in 1:length(traditional.coef.b.noSpill[,'treatment.status'])){
#   trad.treatment.noSpill[j] <- sum(spdf@data$treatment.status * 
#                                      traditional.coef.b.noSpill[, 'treatment.status'][j]) / nrandom
#   spatial.treatment.noSpill[j] <- sum(spdf@data$treatment.status * 
#                                         spatial.coef.b.noSpill[, 'treatment.status'][j]) / nrandom
# }
# 
# 
# tmp.traditional.diff.noSpill <- (true.treatment - trad.treatment.noSpill)
# tmp.spatial.diff.noSpill <- (true.treatment - spatial.treatment.noSpill)
# 
# 
# plot(x=z.vals, type="n", main='Estimated vs. True \n Averge Treatment Effect \n (Model Without Spillover)', ylab="", xlab="z", 
#      ylim=c(min(tmp.traditional.diff.noSpill, tmp.spatial.diff.noSpill), 
#             max(tmp.traditional.diff.noSpill, tmp.spatial.diff.noSpill)),
#      xlim=c(min(z.vals), max(z.vals)))
# 
# lines(z.vals, tmp.traditional.diff.noSpill, col="green")
# lines(z.vals, tmp.spatial.diff.noSpill, col="blue")
# 
# legend(min(z.vals), max(tmp.traditional.diff.noSpill, tmp.spatial.diff.noSpill), 
#        c('spatial', 'traditional'), 
#        lty=c(1), 
#        lwd=c(2.5),
#        col=c('blue', 'green')) 
# 
# 
# 
