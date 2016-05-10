library(sp)
library(ncf)
library(gstat)
library(devtools)
library(rpart)
library(rpart.plot)
library(Rcpp)
library(classInt)
library(RColorBrewer)
library(methods)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp("/sciclone/home00/geogdan/MatchIt/demo/splitc.cpp")
CT_src <- "/sciclone/home00/geogdan/MatchIt/demo/CT_functions.R"
sim_src <- "/sciclone/home00/geogdan/MatchIt/demo/simulation_spatial_data.R"

#detach("package:MatchIt", unload=TRUE)
load_all("/sciclone/home00/geogdan/MatchIt/R")

#1 3800.41856568 0.90376081592 -45.0 45.0 -22.5 22.5 3.21749654825 0.250852506018 0.448021052911 4.27592030555 0.0684864449219 0.29100048171 1 0.330411927736 3.83573033709 1.88067542642 0.698254286741 0.437623061042 10 2.58494466138 /sciclone/home00/geogdan/AlphaSims/test_0.csv 0.954552979835 0.539550663469 0.164665770447
#Args <- c("1", "3800.41856568", "0.90376081592", "-45.0", "45.0", "-22.5", "22.5", "3.21749654825", "0.250852506018", "0.448021052911", "4.27592030555", "0.0684864449219", "0.29100048171", "1", "0.330411927736", "3.83573033709", "1.88067542642", "0.698254286741", "0.437623061042", "10", "2.58494466138", "/sciclone/home00/geogdan/AlphaSims/test_0.csv", "0.954552979835", "0.539550663469", "0.164665770447")
Args <- commandArgs(trailingOnly = TRUE)
print(Args)
out_path=Args[22]
nums = as.numeric(Args)



# -----------------------------------------------------------------------------
# Simulation Settings
# -----------------------------------------------------------------------------

version = Args[1]
nrandom = as.numeric(Args[2])
xvar_error_psill = as.numeric(Args[3])
minx = as.numeric(Args[4])
maxx = as.numeric(Args[5])
miny = as.numeric(Args[6])
maxy = as.numeric(Args[7])

var1.vrange = as.numeric(Args[8])
var1.error = as.numeric(Args[9])
prop_acc = as.numeric(Args[10])
var1_error.vrange = as.numeric(Args[11])
mod_error.magnitude = as.numeric(Args[12])
trt_prc = as.numeric(Args[13])
theta = as.numeric(Args[14])
beta =as.numeric( Args[15])
spill.vrange = as.numeric(Args[16])
spill.magnitude= as.numeric(Args[17])
cal= as.numeric(Args[18])
sample_size = as.numeric(Args[19])
tree_split_lim= as.numeric(Args[20])
mod_error.vrange= as.numeric(Args[21])
xvar_psill=as.numeric(Args[23])
mod_error_psill=as.numeric(Args[24])
trt_spill_sill=as.numeric(Args[25])

p <- 1

iterations <- 1
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

# -----------------------------------------------------------------------------
# Data Simulation
# -----------------------------------------------------------------------------

#Given your settings, outputs:
#spdf@data$modelVar - an ancillary variable with error (defined by var1_error)
#spdf@data$treatment.status - Binary 1/0 treated
#spdf@data$trueOutcome - the measured outcome
#spdf@data$modelOutcome - measured outcome including measurement error (model.error)
source(sim_src)



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
treatment.predictions@data$trueTreatment <- (treatment.predictions@data$treatment.status)*
  (theta + treatment.predictions$trueSpill)


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
#print(correlog.pscore.spillover$x.intercept)
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
db = spdf@data
tree2 = tot.fit.spillB
tree2$frame$yval = as.numeric(rownames(tree2$frame))
res2 = predict(tree2,newdata=spdf@data)

leaf = unique(tot.fit.spillB$where)
res3 = 0
for(i in 1:length(leaf)){
  
  treatcount = 0
  untreatcount = 0
  count = 0
  temp  = 0
  
  for(j in 1:nrow(db)){
    if(res2[j] == leaf[i]){
      temp = c(temp,j)
      count = count + 1
      if(db$treatment.status[j] == 1){
        treatcount = treatcount + 1
      }
      if(db$treatment.status[j] == 0){
        untreatcount = untreatcount + 1
      }
    }
    
  }
  
  if(count == treatcount | count == untreatcount ){
    temp = temp[-1]
    res3 = c(res3,temp)
  }
}

res_leaf2prune = res3[-1]
#res[res_leaf2prune] = NA

## end


#Total Outcome


treatment.predictions@data$tot.spill <-  res * treatment.predictions@data$treatment.status





#CT
source(CT_src)
alist <- list(eval=ctev, split=ctsplit, init=ctinit)

#trans_dta
dbb = trans_dta@data
k = 10
n = dim(dbb)[1]
#sample_size = floor(n)
#ridx = sample(1:n,sample_size,replace=FALSE)
#crxvdata= dbb[ridx,]
crxvdata = dbb
crxvdata$id <- sample(1:k, nrow(crxvdata), replace = TRUE)
list = 1:k
fit1 = rpart(cbind(modelOutcome,treatment.status,m1.pscore,transOutcome) ~ modelVar + coord1 + coord2,
             crxvdata,
             control = rpart.control(cp = 0,minsplit = 5),
             method=alist)
fit = data.matrix(fit1$frame)
index = as.numeric(rownames(fit1$frame))
tsize = dim(fit1$frame[which(fit1$frame$var=="<leaf>"),])[1]

alpha = 0
alphalist = 0
alphalist = cross_validate(fit, index,alphalist)

res = rep(0,length(alphalist)-1)
for(j in 2:(length(alphalist)-1)){
  res[j] = sqrt(alphalist[j]*alphalist[j+1])
}

alphacandidate = res
alphaset = rep(0,length(alphacandidate))
errset = rep(0,length(alphacandidate))
tsize = 0
for(l in 1:length(alphacandidate)){
  alpha = alphacandidate[l]
  error = 0
  treesize = 0
  for (i in 1:k){
    trainingset <- subset(crxvdata, id %in% list[-i])
    testset <- subset(crxvdata, id %in% c(i))
    fit1 = rpart (cbind(modelOutcome,treatment.status,m1.pscore,transOutcome)  ~ modelVar + coord1 + coord2,
                  trainingset,
                  control = rpart.control(cp = alpha,minsplit = 10),
                  method=alist)
    
    if(dim(fit1$frame)[1] == 1){
      error = 0
      break
    }
    
    else{
      treesize = treesize + dim(fit1$frame[which(fit1$frame$var=="<leaf>"),])[1]
      pt = predict(fit1,testset,type = "matrix")
      y = data.frame(pt)
      val = data.matrix(y)
      idx = as.numeric(rownames(y))
      dbidx = as.numeric(rownames(dbb))
      
      for(pid in 1:(dim(y)[1])){
        id = match(idx[pid],dbidx)
        error = error + (crxvdata$transOutcome[id] - val[pid])^2
        #print(error)
      }
    }
    
  }
  
  tsize = c(tsize,treesize/k)
  if(error == 0){
    errset[l] = 1000000
  }
  else{
    errset[l] = error/k
  }
  msg = paste(l,": ",errset[l]*k,sep="")
  #print(msg)
}

tsize = tsize[-1]
alpha_res = alphacandidate[which.min(errset)]
fit_ctpred <- rpart(cbind(modelOutcome,treatment.status,m1.pscore,transOutcome) ~ modelVar + coord1 + coord2,
                    crxvdata, control=rpart.control(minsplit=10,cp=alpha_res),
                    method=alist)
#prp(fit_ctpred)
#res = rep(0,length(alphalist)-1)
#for(j in 2:(length(alphalist)-1)){
#  res[j] = sqrt(alphalist[j]*alphalist[j+1])
#}

#Total Outcome - CT

treatment.predictions@data$ct.spill <-  predict(fit_ctpred,newdata=spdf@data) * treatment.predictions@data$treatment.status
#print("CT Nodes:")
#print(length(unique(fit_ctpred$where)))



#Save Summary Results



for(i in 4:length(treatment.predictions@data))
{
  if(p == 1)
  {
    results[names(treatment.predictions@data)[i]] <- NA
    results["trueTreatment"] <- NA
  }
  results[names(treatment.predictions@data)[i]][p,] <- sum(treatment.predictions@data[,i],na.rm=T) / (nrandom*trt_prc)
  results["trueTreatment"][p,] <- sum(treatment.predictions@data$trueTreatment) / (nrandom*trt_prc)
}






results["spill.magnitude"] <- NA
results["xvar_error_psill"] <- NA
results["mod_error_psill"] <- NA
results["trt_spill_sill"] <- NA
results["xvar_psill"] <- NA
results["var1.vrange"] <- NA
results["prop_acc"] <- NA
results["mod_error.vrange"] <- NA
results["mod_error.magnitude"] <- NA
results["spill.vrange"]<- NA
results["beta"] <- NA
results["theta"] <- NA
results["var1_error.vrange"] <- NA
results["caliper"] <- NA
results["sample_size"]<- NA
results["tree_split_lim"] <- NA

results["spill.magnitude"][p,] <- spill.magnitude
results["xvar_error_psill"][p,]  <- xvar_error_psill
results["mod_error_psill"][p,]  <- mod_error_psill
results["trt_spill_sill"][p,]  <- trt_spill_sill
results["xvar_psill"][p,]  <- xvar_psill
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



#Save Trees
#png()
#prp(fit_ctpred)
#prp(tot.fit.spillB)


#Compare Maps
# map_trt <- treatment.predictions[names(treatment.predictions) != "trueSpill"]
# map_trt <- map_trt[names(map_trt) != "spatial.trueThreshold"]
# map_trt <- map_trt[names(map_trt) != "spatial.matchit.spill"]
# map_trt <- map_trt[names(map_trt) != "baseline"]
# map_trt <- map_trt[names(map_trt) != "baseline.matchit"]
# map_trt <- map_trt[names(map_trt) != "treatment.status"]
# 
# 
# pal = brewer.pal(10,"Greens")
# brks = c(0.0,1,1.25,1.5,1.75,2.0,2.5,3.0,3.5,100)
# spplot(map_trt, zcol=names(map_trt)[names(map_trt) != "id"],cuts=brks,col.regions=pal, col="transparent",
#        main = list(label="Treatment Estimates"), legendEntries=c("0-.25",".25-.5",".5-.75",".75-1","1-1.25","1.25-1.5","1.5-1.75","1.75-2.0","2.0+"), cex=0.5)


#print("Iteration Complete")


#print(proc.time() - ptm)

print(out_path)

write.csv(results,file=out_path)


