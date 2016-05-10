library(sp)
library(ncf)
library(gstat)
library(devtools)
library(rpart)
library(rpart.plot)

#detach("package:MatchIt", unload=TRUE)
#load_all("~/Desktop/Github/MatchIt/R")
library(devtools)
install_github("itpir/matchit")
library(MatchIt)


rm(list = ls())

iterations <- 2

results <- data.frame(
  id=c(1:iterations)
)

results_out <- data.frame(
  id=c(1:iterations)
)


results_nospill <- data.frame(
  id=c(1:iterations)
)


ptm <- proc.time()

  for(p in 1:iterations)
{

  results["id"] <- p
  results_out["id"] <- p
  # -----------------------------------------------------------------------------
  # Simulation Settings
  # -----------------------------------------------------------------------------
  
  
  #General Options
  # set dataframe size (number of points)
  nrandom <- 1000
  
  #General psill
  #Note: setting this to 0 will result in no data randomization.
  #Larger values indicate more autocorrelation.
  psill <- 1 + runif(1, -.95, 20)
  
  # define bounding box
  minx <- -45
  maxx <- 45
  miny <- -22.5
  maxy <- 22.5
  
  #Covariate Spatial Correlation
  var1.vrange <- 1 + runif(1, -.95, 20)
  
  
  #Degree to which the covariate
  #explains the propensity to receive
  #treatment (1.0 = perfect correlation, or no error)
  prop_acc = runif(1, 0.05, .95)
  
  #Spatial pattern in the PSM error, if any 
  #(vrange = 1 approximates random noise)
  var1_error.vrange <- runif(1, 0.1, 20)
  
  #Define the spatial pattern of any model error.
  #Magnitue of error is 0-1 (0 = no error)
  mod_error.vrange <-  1.0 
  mod_error.magnitude <- runif(1, 0.1, 5)
  
  #Percent of locations which are to be defined as
  #treated
  trt_prc = runif(1, 0.2, 0.5)
  
  
  #Theta coefficient for treatment effect
  #used for defining the outcome
  theta <- runif(1, -.95, 10)
  
  #Beta coefficient for ancillary data
  #in defining the treatment
  beta <- runif(1,1,5) * theta
  
  #Spillover Range
  spill.vrange <- 1.0 + runif(1, -.95, 2)
  
  #Spillover Magnitude (relative to theta)
  spill.magnitude <- 1 * runif(1, 0, 3)
  
  #Caliper for Matching
  cal = runif(1, 0.25, 1.0)

  sample_size = runif(1,0.05,1.0)
  
  tree_split_lim = sample(2:50,1) 
  
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
 # source("/home/aiddata/Desktop/Github/MatchIt/demo/test.spatial.figures.R")
  
  
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
  
  nospill.t.pred <- spdf
  nospill.t.pred@data <- nospill.t.pred@data[1]
  nospill.t.pred@data$trueTreatment <- theta
  
  model_dta <- spdf[sample(nrow(spdf), (nrandom * sample_size)), ]
  
  #No Matching
  baseline <- lm(modelOutcome ~ treatment.status +  modelVar, data=model_dta@data)
  outcome.predictions@data$baseline <- predict(baseline, newdata=spdf@data)
  treatment.predictions@data$baseline <- summary(baseline)$coefficients[2]
  nospill.t.pred@data$baseline <- summary(baseline)$coefficients[2]
    
  #Baseline for Comparison
  baseline.matchit <- matchit(treatment.status ~ modelVar, data= model_dta@data,
                    method="nearest", distance="logit", 
                    caliper=cal, calclosest=FALSE, calrandom=FALSE)
  
  baseline.model <- lm(modelOutcome ~ treatment.status +  modelVar, 
                       data=match.data(baseline.matchit))
  
  outcome.predictions@data$baseline.matchit <- predict(baseline.model, newdata=spdf@data)
  treatment.predictions@data$baseline.matchit <- summary(baseline.model)$coefficients[2]
  nospill.t.pred@data$baseline.matchit <- summary(baseline.model)$coefficients[2]
  
  
  #Cheating Spatial PSM - we give the accurate vrange, and use it as a threshold.
  spatial.opts <- list(decay.model = "threshold",
                       threshold = spill.vrange)

  spatial.trueThreshold <- matchit(treatment.status ~ modelVar, data= model_dta,
                    method = "nearest", distance = "logit", 
                    caliper=cal, calclosest=FALSE, calrandom=FALSE,
                    spatial.options=spatial.opts)
  
  spatial.trueThreshold.model <- lm(modelOutcome ~ treatment.status +  modelVar, 
                              data=match.data(spatial.trueThreshold))
  outcome.predictions@data$spatial.trueThreshold <- predict(spatial.trueThreshold.model,
                                                            newdata=spdf@data)
  treatment.predictions@data$spatial.trueThreshold <- 
    summary(spatial.trueThreshold.model)$coefficients[2]
  
  nospill.t.pred@data$spatial.trueThreshold <- summary(spatial.trueThreshold.model)$coefficients[2]
  
  
  #PSM-approximating Traditional and Spatial PSMs

  
  #Propensity Correlogram
  p_cor_spdf <- model_dta
  p_cor_spdf$m1.pscore <- baseline.matchit$distance
  correlog.pscore.spillover <- correlog(x=p_cor_spdf@coords[, 1],
                                        y=p_cor_spdf@coords[, 2],
                                        z=p_cor_spdf$treatment.status,
                                        increment=500,
                                        latlon=TRUE, na.rm=TRUE, resamp=5,
                                        quiet=FALSE)
  
  pscore.spillover.model.dta <- data.frame(
    mean.of.class = c(correlog.pscore.spillover$mean.of.class),
    correlation = c(correlog.pscore.spillover$correlation)
  )
  
  correlog.m1.polynomial <- lm(correlation ~ poly(mean.of.class, degree=5, raw=TRUE),
                               data=pscore.spillover.model.dta)
  
  estimated.spillover.weights <- c()
  for (k in 1:length(p_cor_spdf)) {
    
    correlog.dist <- spDists(p_cor_spdf[k, ]@coords, p_cor_spdf@coords, longlat=TRUE)
    
    
    correlog.neighbors <- correlog.dist > 0
    correlog.treated <- p_cor_spdf$treatment.status
    
    
    
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
  
  #Spatial Matchit
  spatial.opts <- list(decay.model = "threshold",
                       threshold = correlog.pscore.spillover$x.intercept+1)#(correlog.pscore.spillover$x.intercept+1))
  print(correlog.pscore.spillover$x.intercept)
  spatial.matchit.spill <- matchit(treatment.status ~ modelVar, data=p_cor_spdf@data,
                                    method="nearest", distance="logit", 
                                    caliper=cal, calclosest=FALSE, calrandom=FALSE)
  
  spatial.matchit.spill.model <- lm(modelOutcome ~ treatment.status +  modelVar, 
                                     data=match.data(spatial.matchit.spill))
  
  outcome.predictions@data$spatial.matchit.spill <- predict(spatial.matchit.spill.model, 
                                                          newdata=spdf@data)
  treatment.predictions@data$spatial.matchit.spill <- 
    ((summary(spatial.matchit.spill.model)$coefficients[2] * nrandom)) / 
    nrandom
  
  nospill.t.pred@data$spatial.matchit.spill <- summary(spatial.matchit.spill.model)$coefficients[2]
  
  p_cor_spdf@data["coord1"] <- coordinates(p_cor_spdf)[,1]
  p_cor_spdf@data["coord2"] <- coordinates(p_cor_spdf)[,2]
  
  #TOT - Non Random Forest
  trans_dta <- p_cor_spdf
  trans_dta <- trans_dta[(trans_dta@data$m1.pscore != 0 & 
                            trans_dta@data$m1.pscore != 1),]
  
  transOutcome <- list(rep(0,nrow(trans_dta)))
  
  for(i in 1:nrow(p_cor_spdf))
  {
    if(trans_dta$treatment.status[i] == 1)
    {
      transOutcome[i] = trans_dta@data$modelOutcome[i] / trans_dta@data$m1.pscore[i]
    }
    else
    {
      transOutcome[i] = trans_dta@data$modelOutcome[i] / (1-trans_dta@data$m1.pscore[i])
    }
  }
  trans_dta@data$transOutcome <- unlist(transOutcome)
  
  tot.fit.spill <- rpart(transOutcome ~ modelVar + coord1 + coord2,
                         data = trans_dta@data,
                         control=rpart.control(cp=0, minsplit=tree_split_lim),
                         method="anova")
  cpbest <- tot.fit.spill$cptable[which.min(tot.fit.spill$cptable[,"xerror"]),"CP"]
  
  tot.fit.spillB <- rpart(transOutcome ~ modelVar + coord1 + coord2,
                          data = trans_dta@data,
                          control=rpart.control(cp=cpbest, minsplit=tree_split_lim),
                          method="anova")
  
  spdf@data["coord1"] <- coordinates(spdf)[,1]
  spdf@data["coord2"] <- coordinates(spdf)[,2]
  
  res <- predict(tot.fit.spillB, newdata=spdf@data)
  
  # Add by Jianing 
 
  
  ## end
  
  
  #Total Outcome
  outcome.predictions@data$tot.spill <- NA

  treatment.predictions@data$tot.spill <-  res * treatment.predictions@data$treatment.status 

  nospill.t.pred@data$tot.spill <- NA

  
  #Save Summary Results
  
  for(i in 3:length(nospill.t.pred@data))
  {

    if(p == 1)
    {
      results_nospill[names(nospill.t.pred@data)[i]] <- NA
      results_nospill["trueTreatment"] <- NA
    }
    results_nospill[names(nospill.t.pred@data)[i]][p,] <- mean(nospill.t.pred@data[,i])
    results_nospill["trueTreatment"][p,] <- theta
  }
  
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
    results["caliper"] <- NA
    results["sample_size"] <- NA
    results["tree_split_lim"] <- NA
    
    
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
    results_out["caliper"] <- NA
    results_out["sample_size"] <- NA
    results_out["tree_split_lim"] <- NA
    
    results_nospill["spill.magnitude"] <- NA
    results_nospill["psill"] <- NA
    results_nospill["var1.vrange"] <- NA
    results_nospill["prop_acc"] <- NA
    results_nospill["mod_error.vrange"] <- NA
    results_nospill["mod_error.magnitude"] <- NA
    results_nospill["spill.vrange"] <- NA
    results_nospill["beta"] <- NA
    results_nospill["var1_error.vrange"] <- NA
    results_nospill["theta"] <- NA
    results_nospill["caliper"] <- NA
    results_nospill["sample_size"] <- NA
    results_nospill["tree_split_lim"] <- NA
    
    
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
  results["caliper"][p,] <- cal
  results["sample_size"][p,] <- sample_size
  results["tree_split_lim"][p,] <- tree_split_lim
  
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
  results_out["caliper"][p,] <- cal
  results_out["sample_size"][p,] <- sample_size
  results_out["tree_split_lim"][p,] <- tree_split_lim
  
  results_nospill["spill.magnitude"][p,] <- spill.magnitude
  results_nospill["psill"][p,] <- psill
  results_nospill["var1.vrange"][p,] <- var1.vrange
  results_nospill["prop_acc"][p,] <- prop_acc
  results_nospill["mod_error.vrange"][p,] <- mod_error.vrange
  results_nospill["mod_error.magnitude"][p,] <- mod_error.magnitude
  results_nospill["spill.vrange"][p,] <- spill.vrange
  results_nospill["beta"][p,] <- beta
  results_nospill["theta"][p,] <- theta
  results_nospill["var1_error.vrange"][p,] <- var1_error.vrange
  results_nospill["caliper"][p,] <- cal
  results_nospill["sample_size"][p,] <- sample_size
  results_nospill["tree_split_lim"][p,] <- tree_split_lim

  
#Compare Maps
#spplot(outcome.predictions, zcol=names(outcome.predictions)
#       [names(outcome.predictions) != "id"])
#spplot(treatment.predictions, zcol=names(treatment.predictions)
#       [names(treatment.predictions) != "id"])
  print("Iteration Complete")
  print(p)
  print("----------------")
}

print(proc.time() - ptm)



