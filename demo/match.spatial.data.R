library(devtools)
library(sp)

detach("package:MatchIt", unload=TRUE)
load_all("/home/aid_data/Desktop/GitRepo/MatchIt/R")
#library(devtools)
#install_github("itpir/matchit")
library(MatchIt)

###
### An Example Script for Obtaining Matched Data when you have
### Spatial information
###
data(lalonde)

##Simulate Latitude and Longtiude information for each point
set.seed(424)
coords = cbind(runif(614,37.1708,37.3708), runif(614,76.6069,76.8069))

##Create a spatial points data frame
spdf_LL <- SpatialPointsDataFrame(coords, lalonde)

##Traditional, non-spatially weighted matching
m.out1 <- matchit(treat ~ re74 + re75 + age + educ, data = lalonde,
                  method = "nearest", distance = "logit", caliper=.25)

##Matching accounting for spatial spillover and autocorrelation
spatial_opts <- list(spatial.decay.model="GausSemiVar",
                     spatial.thresholds=c(.05))
m.out2 <- matchit(treat ~ re74 + re75 + age + educ, data = spdf_LL,
                  method = "nearest", distance = "logit",
                  spatial.options=spatial_opts, caliper=0.25)

#Next things to work on:
#In spatial case, return a spatial data frame rather than a standard data.frame
#Include a mechanism to automatically model the spatial thresholds
#Expose the funcitonality to generate a correlogram chart, with fitted curves
#Integrate spatial functionality into methods other than nearest.
