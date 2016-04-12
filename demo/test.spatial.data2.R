library(sp)
library(ncf)
library(gstat)

rm(list = ls())

iterations <- 25

results <- data.frame(
  id=c(1:iterations)
)

for(p in 1:iterations)
{
  results["id"] <- p
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
  prop_acc = runif(1, 0, .95)
  
  #Spatial pattern in the PSM error, if any 
  #(vrange = 1 approximates random noise)
  var1_error.vrange <- 0.1
  
  #Define the spatial pattern of any model error.
  #Magnitue of error is 0-1 (0 = no error)
  mod_error.vrange <-  1.0 + runif(1, -.95, 5)
  mod_error.magnitude <- 0.1
  
  #Percent of locations which are to be defined as
  #treated
  trt_prc = .2
  
  #Beta coefficient for ancillary data
  #in defining the treatment
  beta <- 1.0 * runif(1, -5, 5)
  
  #Theta coefficient for treatment effect
  #used for defining the outcome
  theta <- 1.0
  
  #Spillover Range
  spill.vrange <- 1.0 + runif(1, -.95, 5)
  
  #Spillover Magnitude (relative to theta)
  spill.magnitude <- 1 * runif(1, -10, 10)
  
  #Caliper for Matching
  caliper = 1.0
  
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
  
  
  
  
  #Save Summary Results
  for(i in 4:length(treatment.predictions@data))
  {
    if(p == 1)
    {
      results[names(treatment.predictions@data)[i]] <- NA
    }
    results[names(treatment.predictions@data)[i]][p,] <- mean(treatment.predictions@data[,i])
  }

  for(i in 4:length(outcome.predictions@data))
  {
    if(p == 1)
    {
      results[names(outcome.predictions@data)[i]] <- NA
    }
    results[names(outcome.predictions@data)[i]][p,] <- mean(outcome.predictions@data[,i])
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
  }
  results["spill.magnitude"][p,] <- spill.magnitude
  results["psill"][p,] <- psill
  results["var1.vrange"][p,] <- var1.vrange
  results["prop_acc"][p,] <- prop_acc
  results["mod_error.vrange"][p,] <- mod_error.vrange
  results["mod_error.magnitude"][p,] <- mod_error.magnitude
  results["spill.vrange"][p,] <- spill.vrange
  results["beta"][p,] <- beta
  

  
#Compare Maps
#spplot(outcome.predictions, zcol=names(outcome.predictions)
#       [names(outcome.predictions) != "id"])
#spplot(treatment.predictions, zcol=names(treatment.predictions)
#       [names(treatment.predictions) != "id"])

}




#Variable Spillover Magnitude
results.plot <- results[order(results$spill.magnitude),]
plot(ylim=c(-20,20), results.plot$spill.magnitude, results.plot$baseline, col="red", main="ATE by Spillover")
lines(results.plot$spill.magnitude,results.plot$baseline, col="red")
lines(results.plot$spill.magnitude, 
      results.plot$trueTreatment, col="green")
points(results.plot$spill.magnitude, 
       results.plot$trueTreatment, col="green")

lines(results.plot$spill.magnitude, 
      results.plot$baseline.matchit, col="blue")
points(results.plot$spill.magnitude, 
       results.plot$baseline.matchit , col="blue")

legend("topright", legend=c("Baseline LM","True ATE", "Baseline MatchIt"), pch=c(pch = 1, pch=1, pch=1),
       col=c(col="red", col="green", col="blue"), title = "Legend")


#Variable V1 Range
results.plot <- results[order(results$var1.vrange),]
plot(ylim=c(-20,20), results.plot$var1.vrange, results.plot$baseline, col="red", main="ATE by Beta Autocor.")
lines(results.plot$var1.vrange,results.plot$baseline, col="red")
lines(results.plot$var1.vrange, 
      results.plot$trueTreatment, col="green")
points(results.plot$var1.vrange, 
       results.plot$trueTreatment, col="green")

lines(results.plot$var1.vrange, 
      results.plot$baseline.matchit, col="blue")
points(results.plot$var1.vrange, 
       results.plot$baseline.matchit , col="blue")

legend("topright", legend=c("Baseline LM","True ATE", "Baseline MatchIt"), pch=c(pch = 1, pch=1, pch=1),
       col=c(col="red", col="green", col="blue"), title = "Legend")

#Sill
results.plot <- results[order(results$psill),]
plot(ylim=c(-20,20), results.plot$psill, results.plot$baseline, col="red", main="ATE by Sill")
lines(results.plot$psill,results.plot$baseline, col="red")
lines(results.plot$psill, 
      results.plot$trueTreatment, col="green")
points(results.plot$psill, 
       results.plot$trueTreatment, col="green")

lines(results.plot$psill, 
      results.plot$baseline.matchit, col="blue")
points(results.plot$psill, 
       results.plot$baseline.matchit , col="blue")

legend("topright", legend=c("Baseline LM","True ATE", "Baseline MatchIt"), pch=c(pch = 1, pch=1, pch=1),
       col=c(col="red", col="green", col="blue"), title = "Legend")

#prop_acc
results.plot <- results[order(results$prop_acc),]
plot(ylim=c(-20,20), results.plot$prop_acc, results.plot$baseline, col="red", main="ATE by Treatment/Beta Corr")
lines(results.plot$prop_acc,results.plot$baseline, col="red")
lines(results.plot$prop_acc, 
      results.plot$trueTreatment, col="green")
points(results.plot$prop_acc, 
       results.plot$trueTreatment, col="green")

lines(results.plot$prop_acc, 
      results.plot$baseline.matchit, col="blue")
points(results.plot$prop_acc, 
       results.plot$baseline.matchit , col="blue")

legend("topright", legend=c("Baseline LM","True ATE", "Baseline MatchIt"), pch=c(pch = 1, pch=1, pch=1),
       col=c(col="red", col="green", col="blue"), title = "Legend")

