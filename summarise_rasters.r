# this file combines and summarises output rasters from the simulations
library(raster)

raster_dir     <- "/short/ka2/dfr805/simulation/test_runs/real225_nicherate0.02_dispersal1.5_200gens/raster/"

file_template  <- "*time100.asc"

doPlot <- FALSE

files <- list.files(path=raster_dir, pattern=file_template, include.dirs = TRUE, full.names = TRUE)

abund.stack <- stack(files)
layer.count <- nlayers(abund.stack)

cat("\nThere are", layer.count, "simulation results to summarise for:\n", raster_dir, "\n", file_template, "\n")

abund.ras    <- sum(abund.stack)
presence.stack <- abund.stack
presence.stack[which(presence.stack[]>0)] <- 1
richness.ras <- sum(presence.stack)

speciesEnd.stack <- presence.stack
for (lyr in 1:layer.count) {
  layer_sum <- sum(speciesEnd.stack[[lyr]][],na.rm=T)
  speciesEnd.stack[[lyr]] <- speciesEnd.stack[[lyr]] / layer_sum
}
speciesEnd.ras <- sum(speciesEnd.stack, na.rm=T)

abundEnd.stack <- abund.stack
for (lyr in 1:layer.count) {
  layer_sum <- sum(abundEnd.stack[[lyr]][],na.rm=T)
  abundEnd.stack[[lyr]] <- abundEnd.stack[[lyr]] / layer_sum
}
abundEnd.ras <- sum(abundEnd.stack, na.rm=T)

writeRaster(abund.ras, paste(raster_dir, "abundance_at_100.asc", sep=""), overwrite=TRUE)
writeRaster(richness.ras, paste(raster_dir, "richness_at_100.asc", sep=""), overwrite=TRUE)
writeRaster(speciesEnd.ras, paste(raster_dir, "speciesEnd_at_100.asc", sep=""), overwrite=TRUE)
writeRaster(abundEnd.ras, paste(raster_dir, "abundanceEnd_at_100.asc", sep=""), overwrite=TRUE)


if (doPlot) {

  # lines for windows desktop
  abund.sum   <- raster("C:/Work/Simulation/test_runs_Sept/real225_nicherate0.02_dispersal1.5_200gens/raster/abundance_at_100.asc")
  plot(abund.sum)

  richness.sum   <- raster("C:/Work/Simulation/test_runs_Sept/real225_nicherate0.02_dispersal1.5_200gens/raster/richness_at_100.asc")
  plot(richness.sum)

  speciesEndemism.sum   <- raster("C:/Work/Simulation/test_runs_Sept/real225_nicherate0.02_dispersal1.5_200gens/raster/speciesEnd_at_100.asc")
  plot(speciesEndemism.sum)

  abundanceEndemism.sum   <- raster("C:/Work/Simulation/test_runs_Sept/real225_nicherate0.02_dispersal1.5_200gens/raster/abundanceEnd_at_100.asc")
  plot(log(abundanceEndemism.sum))
}
