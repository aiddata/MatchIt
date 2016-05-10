#Folder that contains CSVs to be visualized
res_folder <- "/mnt/sc/AlphaSims"
files <- list.files(path=res_folder, full.names=T, recursive=FALSE)

print(length(files))

for(path in 1:length(files))
{
  if(path == 1)
  {
    full_df <- read.csv(files[path])
  }
  else
  {
    full_df <-  rbind(full_df,read.csv(files[path]))
  }
}

results <- full_df

#Visualization of Simulation results
library(plotly)
viz.sims <- function(results, varH, mtitle, pre="")
{
  results.plot <- results

  
  eval(parse(text=paste("results.plot$v1 <- results.plot$",varH,sep="")))
  results.plot <- results.plot[order(results.plot$v1),]
  ylower <- min(results.plot[paste(pre,"trueTreatment",sep="")]) - (2*abs(min(results.plot[paste(pre,"trueTreatment",sep="")])))
  yupper <- max(results.plot[paste(pre,"baseline",sep="")]) 
  if(mtitle == "ATE by Model")
  {
    yupper <- max(results.plot[paste(pre,"baseline",sep="")]) * 1.25
  }
  plot(ylim=c(ylower,yupper), 
       results.plot$v1, 
       results.plot[paste(pre,"baseline",sep="")][[1]], 
       col=rgb(1,0,0,alpha=0.1), pch=3, cex=0.5,
       main=mtitle,
       ylab="Estimate",
       xlab=varH)
  lines(lowess(results.plot$v1,results.plot[paste(pre,"baseline",sep="")][[1]]), col=rgb(1,0,0), pch=3)
  
  lines(lowess(results.plot$v1, 
        results.plot[paste(pre,"baseline.matchit",sep="")][[1]]), col=rgb(0,0,1), pch=4)
  points(results.plot$v1, 
         results.plot[paste(pre,"baseline.matchit",sep="")][[1]], col=rgb(0,0,1,alpha=0.1), pch=4, cex=0.5)
  
  lines(lowess(results.plot$v1, 
               results.plot[paste(pre,"trueTreatment",sep="")][[1]]), col=rgb(0,1,0), pch=4)
  points(results.plot$v1, 
         results.plot[paste(pre,"trueTreatment",sep="")][[1]], col=rgb(0,1,0,alpha=0.1), pch=4, cex=0.5)
  
  lines(lowess(results.plot$v1, 
               results.plot[paste(pre,"ct.spill",sep="")][[1]]), col=rgb(1,0.5,0), pch=2)
  points(results.plot$v1, 
         results.plot[paste(pre,"ct.spill",sep="")][[1]], col=rgb(1,0.5,0,alpha=0.1), pch=2, cex=0.5)

  #lines(lowess(results.plot$v1, 
  #             results.plot[paste(pre,"matchit.spill.model",sep="")][[1]]), col=rgb(0,1,1), pch=4)
  #points(results.plot$v1, 
  #       results.plot[paste(pre,"matchit.spill.model",sep="")][[1]], col=rgb(0,1,1,alpha=0.1), pch=4, cex=0.5)
  
  lines(lowess(results.plot$v1, 
               results.plot[paste(pre,"spatial.matchit.spill",sep="")][[1]]), col=rgb(0.5,0.5,0.5), pch=4)
  points(results.plot$v1, 
         results.plot[paste(pre,"spatial.matchit.spill",sep="")][[1]], col=rgb(0.5,0.5,0.5,alpha=0.1), pch=4, cex=0.5)
  
  if(mtitle == "ATE by Model")
  {
    lines(lowess(results.plot$v1, 
               results.plot[paste(pre,"tot.spill",sep="")][[1]]), col=rgb(0,0,0), pch=3)
  points(results.plot$v1, 
         results.plot[paste(pre,"tot.spill",sep="")][[1]], col=rgb(0,0,0,alpha=0.1), pch=3, cex=0.5)
  }
  
  legend("topleft",
         cex = 0.65,
         legend=c("Baseline LM","True ATE", "Baseline MatchIt", 
                  "CT", "TOT", 
                  "NA", "Spatial Thresh"), 
         pch=c(pch = 3, pch=1, pch=4, pch=2, pch=3, pch=4, pch=2),
         col=c(col="red", col="green", col="blue", col="orange", 
               col="black", col=109, col=144), title = "Legend")
}

#Variable Spillover Magnitude
viz.sims(results, "spill.magnitude", "ATE by Model")
viz.sims(results, "var1.vrange", "ATE by Model")
viz.sims(results, "prop_acc", "ATE by Model")
viz.sims(results, "spill.vrange", "ATE by Model")
viz.sims(results, "caliper", "ATE by Model")
viz.sims(results, "sample_size", "ATE by Model")
viz.sims(results, "mod_error.magnitude", "ATE by Model")
viz.sims(results, "tree_split_lim", "ATE by Model")

#viz.sims(results_out, "spill.magnitude", "Predicted Avg. Outcome")
#viz.sims(results_out, "var1.vrange", "Predicted Avg. Outcome")
#viz.sims(results_out, "psill", "Predicted Avg. Outcome")
#viz.sims(results_out, "prop_acc", "Predicted Avg. Outcome")
#viz.sims(results_out, "spill.vrange", "Predicted Avg. Outcome")
#viz.sims(results_out, "caliper", "Predicted Avg. Outcome")

#viz.sims(results_nospill, "spill.magnitude", "Predicted Theta")
#viz.sims(results_nospill, "var1.vrange", "Predicted Theta") #5
#viz.sims(results_nospill, "psill", "Predicted Theta") #200
#viz.sims(results_nospill, "prop_acc", "Predicted Theta") #0.95
#viz.sims(results_nospill, "spill.vrange", "Predicted Theta") #1.0
#viz.sims(results_nospill, "caliper", "Predicted Theta")
#viz.sims(results_nospill, "theta", "Predicted Theta") #1.0


#Compare to Truth
dif_func <- function(results, comparison)
{
  results.dif <- results
  abs.comp <- paste("results.dif$dif.abs.", comparison, 
                    "<- abs(results.dif$trueTreatment - results.dif$",
                    comparison,")", sep="")
  rel.comp <- paste("results.dif$dif.", comparison, 
                    "<- (results.dif$trueTreatment - results.dif$",
                    comparison,")", sep="")
  eval(parse(text=abs.comp))
  eval(parse(text=rel.comp))
  return(results.dif)
  
}

results <- dif_func(results, "baseline")
results <- dif_func(results, "baseline.matchit")
results <- dif_func(results, "spatial.matchit.spill")
results <- dif_func(results, "spatial.trueThreshold")
results <- dif_func(results, "trueTreatment")
results <- dif_func(results, "tot.spill")
results <- dif_func(results, "ct.spill")



viz.sims(results, "spill.magnitude", "Abs Dif ATE by Model", "dif.abs.")
viz.sims(results, "var1.vrange", "Abs Dif ATE by Model", "dif.abs.")
viz.sims(results, "psill", "Abs Dif ATE by Model", "dif.abs.")
viz.sims(results, "prop_acc", "Abs Dif ATE by Model", "dif.abs.")
viz.sims(results, "spill.vrange", "Abs Dif ATE by Model", "dif.abs.")
viz.sims(results, "caliper", "Abs Dif ATE by Model", "dif.abs.")
viz.sims(results, "sample_size", "Abs Dif ATE by Model", "dif.abs.")

#viz.sims(results, "spill.magnitude", "ATE by Model", "dif.")
#viz.sims(results, "var1.vrange", "ATE by Model", "dif.")
#viz.sims(results, "psill", "ATE by Model", "dif.")
#viz.sims(results, "prop_acc", "ATE by Model", "dif.")
#viz.sims(results, "spill.vrange", "ATE by Model", "dif.")
#viz.sims(results, "caliper", "ATE by Model", "dif.")



# plot_ly(results, x=spill.magnitude, y=spill.vrange, 
#         z=baseline, type="scatter3d", mode="markers", opacity=0.5, name="Baseline")
# 
# add_trace(results, x=spill.magnitude, y=spill.vrange, 
#           z=spatial.trueThreshold, 
#           type="scatter3d",
#           mode="markers",
#           opacity=0.5, 
#           name="Spatial Threshold")
# 
# add_trace(results, x=spill.magnitude, y=spill.vrange, 
#           z=trueTreatment, 
#           type="scatter3d",
#           mode="markers",
#           opacity=0.5, 
#           name="True (0 Surface)")
# 
# add_trace(results, x=spill.magnitude, y=spill.vrange, 
#           z=spatial.matchit.spill, 
#           type="scatter3d",
#           mode="markers",
#           opacity=0.5, 
#           name="Spatial Spill")
