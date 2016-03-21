library(devtools)
library(sp)

detach("package:MatchIt", unload=TRUE)
load_all("~/git/matchit/R")
#library(devtools)
#install_github("itpir/matchit")
library(MatchIt)

library(ncf)
# # for genetic
# library(rgenoud)
library(Matching)
# # for cem
library(cem)
library(Zelig)

# # for demo/analysis.R
# # install Rgraphviz dependency of MCMCpack which is dependency of Zelig
# source("https://bioconductor.org/biocLite.R")
# biocLite("Rgraphviz")
# library(MCMCpack)
# library(Zelig)

###
### An Example Script for Obtaining Matched Data when you have
### Spatial information
###
data(lalonde)

##Simulate Latitude and Longtiude information for each point
set.seed(424)
coords = cbind(runif(nrow(lalonde),37.1708,37.3708), runif(nrow(lalonde),76.6069,76.8069))

##Create a spatial points data frame
spdf_LL <- SpatialPointsDataFrame(coords, lalonde)

##Traditional, non-spatially weighted matching
m.out1 <- matchit(treat ~ re74 + re75 + age + educ, data=lalonde,
                  method="nearest", distance="logit", caliper=0.25)

##Matching accounting for spatial spillover and autocorrelation
spatial_opts <- list(decay.model = "gaussian.semivariance",
                     thresholds = c(.05),
                     caliper = 0.25)

m.out2 <- matchit(treat ~ re74 + re75 + age + educ, data = spdf_LL,
                  method = "nearest", distance = "logit",
                  spatial.options=spatial_opts)

#Next things to work on:
#In spatial case, return a spatial data frame rather than a standard data.frame
#Include a mechanism to automatically model the spatial thresholds
#Expose the funcitonality to generate a correlogram chart, with fitted curves
#Integrate spatial functionality into methods other than nearest.


