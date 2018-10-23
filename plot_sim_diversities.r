library(raster)

setwd("C:/users/u3579238/Dropbox/Simulation/test_runs_July2018/real225_nicherate0.02_dispersal1.5_200gens_v2")

steps <- 1

rich.ras <- raster(paste("richness_at_", steps, ".asc", sep=''))
abund.ras <- raster(paste("abundance_at_", steps, ".asc", sep=''))
sp_end.ras <- raster(paste("speciesEnd_at_", steps, ".asc", sep=''))
abund_end.ras <- raster(paste("abundanceEnd_at_", steps, ".asc", sep=''))

windows()
par(mfcol=c(2,2))

plot(rich.ras, main=paste("Richness ", steps), cex.main=0.7)
plot(abund.ras, main=paste("Abundance ", steps), cex.main=0.7)
plot(sp_end.ras, main=paste("Species endemism (WE) ", steps), cex.main=0.7)
plot(rich.ras, main=paste("Abundance endemism ", steps), cex.main=0.7)

