library(sp)

#qtest branch
#TemporaryWorkingFile to add random latitude and longitude data
#to the lalonde dataset to illustrate mechanisms for spatial matching.

LL <- read.table('/home/aid_data/Desktop/GitRepo/MatchIt/data/lalonde.tab.gz')

coords = cbind(runif(614,37.1708,37.3708), runif(614,76.6069,76.8069))

spdf_LL <- SpatialPointsDataFrame(coords, LL)

write.table(spdf_LL,'/home/aid_data/Desktop/GitRepo/MatchIt/data/lalonde_spatial.tab.gz' )
