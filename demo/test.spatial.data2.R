library(sp)
library(ncf)
library(gstat)

rm(list = ls())

iterations <- 5

results <- data.frame(
  id=c(1:iterations)
)

results_out <- data.frame(
  id=c(1:iterations)
)

for(p in 1:iterations)
{
  results["id"] <- p
  results_out["id"] <- p
  # -----------------------------------------------------------------------------
  # Simulation Settings
  # -----------------------------------------------------------------------------
  
  
  #General Options
  # set dataframe size (number of points)
  nrandom <- 2000
  
  #General psill
  #Note: setting this to 0 will result in no data randomization.
  #Larger values indicate more autocorrelation.
  psill <- 1.0 + runif(1, -.95, 10)
  
  # define bounding box
  minx <- -45
  maxx <- 45
  miny <- -22.5
  maxy <- 22.5
  
  #Covariate Spatial Correlation
  var1.vrange <- 1.0 + runif(1, -.95, 5)
  
  
  #Degree to which the covariate
  #explains the propensity to receive
  #treatment (1.0 = perfect correlation, or no error)
  prop_acc = runif(1, 0.5, .95)
  
  #Spatial pattern in the PSM error, if any 
  #(vrange = 1 approximates random noise)
  var1_error.vrange <- runif(1, 0.1, .95)
  
  #Define the spatial pattern of any model error.
  #Magnitue of error is 0-1 (0 = no error)
  mod_error.vrange <-  1.0 
  mod_error.magnitude <- 0.0 
  
  #Percent of locations which are to be defined as
  #treated
  trt_prc = .2
  
  #Beta coefficient for ancillary data
  #in defining the treatment
  beta <- 1.0 
  
  #Theta coefficient for treatment effect
  #used for defining the outcome
  theta <- 1.0 
  
  #Spillover Range
  spill.vrange <- 1.0 + runif(1, -.95, 10)
  
  #Spillover Magnitude (relative to theta)
  spill.magnitude <- 1 * runif(1, 0, 50)
  
  #Caliper for Matching
  caliper = 0.5
  
  # -----------------------------------------------------------------------------
  # Data Simulation
  # -----------------------------------------------------------------------------
  
  #Given your settings, outputs:
  #spdf@data$modelVar - an ancillary variable with error (defined by var1_error)
  #spdf@data$treatment.status - Binary 1/0 treated
  #spdf@data$trueOutcome - the measured outcome
  #spdf@data$modelOutcome - measured outcome including measurement error (model.error)
  source("/home/aiddata/Desktop/Github/MatchIt/demo/simulation_spatial_data.R")
  
  
  
  # -----------------------------------------------------------------------------
  # Data Generation Visualizations (Per-iteration; currently not saved)
  # Disable for a large speed boost.
  # Note visualizations use MatchIt for PSM calculations (nearest/logit/caliper=0.25)
  # -----------------------------------------------------------------------------
  #Creates a figure describing your iteration (PSM vs. observed and parameterized spillovers
  #and maps).
  #source("/home/aiddata/Desktop/Github/MatchIt/demo/test.spatial.figures.R")
  
  
  # -----------------------------------------------------------------------------
  # Model Tests
  # Save all Model Predictions into a SPDF for comparison.
  # -----------------------------------------------------------------------------
  outcome.predictions <- spdf
  outcome.predictions@data <- outcome.predictions@data[c(1,8)]
  
  treatment.predictions <- spdf
  treatment.predictions@data <- treatment.predictions@data[c(1,6,7)]
  treatment.predictions@data$trueTreatment <- (treatment.predictions@data$treatment.status *
                                               theta) + treatment.predictions$trueSpill
  
  
  #No Matching
  baseline <- lm(modelOutcome ~ treatment.status +  modelVar, data=spdf@data)
  outcome.predictions@data$baseline <- predict(baseline, newdata=spdf@data)
  treatment.predictions@data$baseline <- summary(baseline)$coefficients[2]
  
  #Baseline for Comparison
  baseline.matchit <- matchit(treatment.status ~ modelVar, data=spdf@data,
                    method="nearest", distance="logit", 
                    caliper=caliper, calclosest=FALSE, calrandom=FALSE)
  
  baseline.model <- lm(modelOutcome ~ treatment.status +  modelVar, 
                       data=match.data(baseline.matchit))
  
  outcome.predictions@data$baseline.matchit <- predict(baseline.model, newdata=spdf@data)
  treatment.predictions@data$baseline.matchit <- summary(baseline.model)$coefficients[2]
  
  
  #Cheating Spatial PSM - we give the accurate vrange, and use it as a threshold.
  spatial.opts <- list(decay.model = "threshold",
                       threshold = spill.vrange)

  spatial.trueThreshold <- matchit(treatment.status ~ modelVar, data=spdf,
                    method = "nearest", distance = "logit", 
                    caliper=caliper, calclosest=FALSE, calrandom=FALSE,
                    spatial.options=spatial.opts)
  
  spatial.trueThreshold.model <- lm(modelOutcome ~ treatment.status +  modelVar, 
                              data=match.data(spatial.trueThreshold))
  outcome.predictions@data$spatial.trueThreshold <- predict(spatial.trueThreshold.model,
                                                            newdata=spdf@data)
  treatment.predictions@data$spatial.trueThreshold <- 
    summary(spatial.trueThreshold.model)$coefficients[2]
  
  
  
  
  #PSM-approximating Traditional and Spatial PSMs

  
  #Propensity Correlogram
  p_cor_spdf <- spdf
  p_cor_spdf$m1.pscore <- baseline.matchit$distance
  correlog.pscore.spillover <- correlog(x=p_cor_spdf@coords[, 1],
                                        y=p_cor_spdf@coords[, 2],
                                        z=p_cor_spdf$m1.pscore,
                                        increment=500,
                                        latlon=TRUE, na.rm=TRUE, resamp=2,
                                        quiet=FALSE)
  
  pscore.spillover.model.dta <- data.frame(
    mean.of.class = c(correlog.pscore.spillover$mean.of.class),
    correlation = c(correlog.pscore.spillover$correlation)
  )
  
  correlog.m1.polynomial <- lm(correlation ~ poly(mean.of.class, degree=5, raw=TRUE),
                               data=pscore.spillover.model.dta)
  
  estimated.spillover.weights <- c()
  for (k in 1:nrandom) {
    
    correlog.dist <- spDists(spdf[k, ]@coords, spdf@coords, longlat=TRUE)
    
    
    correlog.neighbors <- correlog.dist > 0
    correlog.treated <- spdf$treatment.status
    
    correlog.newdata <- data.frame(
      mean.of.class = c(correlog.dist)
    )
    
    #Limit to the min (after the final predicted distance all values are 0)
    correlog.newdata$weight_dist <- predict(correlog.m1.polynomial, newdata=correlog.newdata)
    correlog.newdata$weight_dist[correlog.newdata$weight_dist < 0] <- 0
    
    estimated.spillover.weights[k] <- sum(correlog.newdata$weight_dist * correlog.neighbors * 
                                           correlog.treated)
    
  }
  
  p_cor_spdf$spillover.est <- estimated.spillover.weights
  
  #Baseline with Spillover Term
  baseline.spill.model <- lm(modelOutcome ~ treatment.status +  modelVar + spillover.est, 
                       data=p_cor_spdf)
  
  outcome.predictions@data$baseline.spill.model <- predict(baseline.spill.model, 
                                                           newdata=p_cor_spdf@data)
  treatment.predictions@data$baseline.spill.model <- 
    summary(baseline.spill.model)$coefficients[2]
  
  
  #Baseline Matchit with Spillover
  baseline.matchit.spill <- matchit(treatment.status ~ modelVar, data=p_cor_spdf@data,
                              method="nearest", distance="logit", 
                              caliper=caliper, calclosest=FALSE, calrandom=FALSE)
  
  baseline.matchit.spill.model <- lm(modelOutcome ~ treatment.status +  modelVar + spillover.est, 
                             data=match.data(baseline.matchit.spill))
  
  outcome.predictions@data$matchit.spill.model <- predict(baseline.matchit.spill.model, 
                                                           newdata=p_cor_spdf@data)
  treatment.predictions@data$matchit.spill.model <- 
    summary(baseline.matchit.spill.model)$coefficients[2]
  
  #Spatial Matchit with Spillover
  spatial.opts <- list(decay.model = "spherical",
                       threshold = (correlog.pscore.spillover$x.intercept+100))
  
  spatial.matchit.spill <- matchit(treatment.status ~ modelVar, data=p_cor_spdf@data,
                                    method="nearest", distance="logit", 
                                    caliper=caliper, calclosest=FALSE, calrandom=FALSE)
  
  spatial.matchit.spill.model <- lm(modelOutcome ~ treatment.status +  modelVar + spillover.est, 
                                     data=match.data(spatial.matchit.spill))
  
  outcome.predictions@data$spatial.matchit.spill <- predict(spatial.matchit.spill.model, 
                                                          newdata=p_cor_spdf@data)
  treatment.predictions@data$spatial.matchit.spill <- 
    summary(spatial.matchit.spill.model)$coefficients[2]
  
  
  #Save Summary Results
  for(i in 4:length(treatment.predictions@data))
  {
    if(p == 1)
    {
      results[names(treatment.predictions@data)[i]] <- NA
      results["trueTreatment"] <- NA
    }
    results[names(treatment.predictions@data)[i]][p,] <- mean(treatment.predictions@data[,i])
    results["trueTreatment"][p,] <- mean(treatment.predictions@data$trueTreatment)
  }

  for(i in 3:length(outcome.predictions@data))
  {
    if(p == 1)
    {
      results_out[names(outcome.predictions@data)[i]] <- NA
      results_out["trueTreatment"] <- NA
    }
    results_out[names(outcome.predictions@data)[i]][p,] <- mean(outcome.predictions@data[,i])
    results_out["trueTreatment"][p,] <- mean(treatment.predictions@data$trueTreatment)
    }

  #Save relevant parameters
  if(p == 1)
  {
    results["spill.magnitude"] <- NA
    results["psill"] <- NA
    results["var1.vrange"] <- NA
    results["prop_acc"] <- NA
    results["mod_error.vrange"] <- NA
    results["mod_error.magnitude"] <- NA
    results["spill.vrange"] <- NA
    results["beta"] <- NA
    results["var1_error.vrange"] <- NA
    results["theta"] <- NA
    
    results_out["spill.magnitude"] <- NA
    results_out["psill"] <- NA
    results_out["var1.vrange"] <- NA
    results_out["prop_acc"] <- NA
    results_out["mod_error.vrange"] <- NA
    results_out["mod_error.magnitude"] <- NA
    results_out["spill.vrange"] <- NA
    results_out["beta"] <- NA
    results_out["var1_error.vrange"] <- NA
    results_out["theta"] <- NA
  }
  results["spill.magnitude"][p,] <- spill.magnitude
  results["psill"][p,] <- psill
  results["var1.vrange"][p,] <- var1.vrange
  results["prop_acc"][p,] <- prop_acc
  results["mod_error.vrange"][p,] <- mod_error.vrange
  results["mod_error.magnitude"][p,] <- mod_error.magnitude
  results["spill.vrange"][p,] <- spill.vrange
  results["beta"][p,] <- beta
  results["theta"][p,] <- theta
  results["var1_error.vrange"][p,] <- var1_error.vrange
  
  results_out["spill.magnitude"][p,] <- spill.magnitude
  results_out["psill"][p,] <- psill
  results_out["var1.vrange"][p,] <- var1.vrange
  results_out["prop_acc"][p,] <- prop_acc
  results_out["mod_error.vrange"][p,] <- mod_error.vrange
  results_out["mod_error.magnitude"][p,] <- mod_error.magnitude
  results_out["spill.vrange"][p,] <- spill.vrange
  results_out["beta"][p,] <- beta
  results_out["theta"][p,] <- theta
  results_out["var1_error.vrange"][p,] <- var1_error.vrange
  
  

  
#Compare Maps
#spplot(outcome.predictions, zcol=names(outcome.predictions)
#       [names(outcome.predictions) != "id"])
#spplot(treatment.predictions, zcol=names(treatment.predictions)
#       [names(treatment.predictions) != "id"])

}


#visualizaiton function
viz.sims <- function(results, varH, mtitle)
{
  results.plot <- results
  eval(parse(text=paste("results.plot$v1 <- results.plot$",varH,sep="")))
  results.plot <- results.plot[order(results.plot$v1),]
  ylower <- min(results.plot$baseline)-min(results.plot$baseline)
  yupper <- max(results.plot$baseline)*2
  plot(ylim=c(ylower,yupper), 
       results.plot$v1, 
       results.plot$baseline, 
       col="red", pch=3, 
       main=mtitle,
       ylab="Estimate",
       xlab=varH)
  lines(results.plot$v1,results.plot$baseline, col="red", pch=3)
  
  lines(results.plot$v1, 
        results.plot$baseline.matchit, col="blue", pch=4)
  points(results.plot$v1, 
         results.plot$baseline.matchit , col="blue", pch=4)
  
  lines(results.plot$v1, 
        results.plot$trueTreatment, col="green")
  points(results.plot$v1, 
         results.plot$trueTreatment, col="green")
  
  lines(results.plot$v1, 
        results.plot$spatial.trueThreshold, col="orange", pch=2)
  points(results.plot$v1, 
         results.plot$spatial.trueThreshold , col="orange", pch=2)
  
  lines(results.plot$v1, 
        results.plot$baseline.spill.model, col="black", pch=3)
  points(results.plot$v1, 
         results.plot$baseline.spill.model , col="black", pch=3)
  
  lines(results.plot$v1, 
        results.plot$matchit.spill.model, col=109, pch=4)
  points(results.plot$v1, 
         results.plot$matchit.spill.model , col=109, pch=4)
  
  lines(results.plot$v1, 
        results.plot$spatial.matchit.spill, col=144, pch=2)
  points(results.plot$v1, 
         results.plot$spatial.matchit.spill , col=144, pch=2)
  
  legend("topleft",
         cex = 0.65,
         legend=c("Baseline LM","True ATE", "Baseline MatchIt", 
                  "Spatial True Thresh", "Baseline LM Spill", 
                  "Matchit Spill", "Spatial Spill"), 
         pch=c(pch = 3, pch=1, pch=4, pch=2, pch=3, pch=4, pch=2),
         col=c(col="red", col="green", col="blue", col="orange", 
               col="black", col=109, col=144), title = "Legend")
}

#Variable Spillover Magnitude
viz.sims(results, "spill.magnitude", "ATE by Model")
viz.sims(results, "var1.vrange", "ATE by Model")
viz.sims(results, "psill", "ATE by Model")
viz.sims(results, "prop_acc", "ATE by Model")
viz.sims(results, "spill.vrange", "ATE by Model")


viz.sims(results_out, "spill.magnitude", "Predicted Avg. Outcome")
viz.sims(results_out, "var1.vrange", "Predicted Avg. Outcome")
viz.sims(results_out, "psill", "Predicted Avg. Outcome")
viz.sims(results_out, "prop_acc", "Predicted Avg. Outcome")
viz.sims(results_out, "spill.vrange", "Predicted Avg. Outcome")







