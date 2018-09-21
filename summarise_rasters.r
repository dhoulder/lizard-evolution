# this file combines and summarises output rasters from the simulations
library(raster)

raster_dir    <- "/short/ka2/dfr805/simulation/test_runs/real225_nicherate0.02_dispersal1.5_200gens/raster/"

summary.times <- c(1, 10, 50, 100, 150, 200)
#summary.time  <- 50

for (summary.time in summary.times) {

  file_template <- paste("*time", summary.time, ".asc", sep="")

  doPlot <- FALSE

  files <- list.files(path=raster_dir, pattern=file_template, include.dirs = TRUE, full.names = TRUE)

  abund.stack <- stack(files)
  layer.count <- nlayers(abund.stack)

  cat("\nThere are", layer.count, "simulation results to summarise for:\n", raster_dir, "\n", file_template, "\n\n")

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

  writeRaster(abund.ras, paste(raster_dir, "abundance_at_", summary.time ,".asc", sep=""), overwrite=TRUE)
  writeRaster(richness.ras, paste(raster_dir, "richness_at_", summary.time ,".asc", sep=""), overwrite=TRUE)
  writeRaster(speciesEnd.ras, paste(raster_dir, "speciesEnd_at_", summary.time ,".asc", sep=""), overwrite=TRUE)
  writeRaster(abundEnd.ras, paste(raster_dir, "abundanceEnd_at_", summary.time ,".asc", sep=""), overwrite=TRUE)


  if (doPlot) {

    windows(10,10)
    par(mfcol=c(2,2))

    # lines for windows desktop
    abund.sum   <- raster(paste("C:/Work/Simulation/test_runs_Sept/real225_nicherate0.02_dispersal1.5_200gens/raster/abundance_at_", summary.time, ".asc", sep=""))
    main_header <- paste("Abundance after", summary.time, "steps from", layer.count, "runs")
    plot(abund.sum, main=main_header)

    richness.sum   <- raster(paste("C:/Work/Simulation/test_runs_Sept/real225_nicherate0.02_dispersal1.5_200gens/raster/richness_at_", summary.time, ".asc", sep=""))
    main_header <- paste("Species richness after", summary.time, "steps")
    plot(richness.sum, main=main_header)

    speciesEndemism.sum   <- raster(paste("C:/Work/Simulation/test_runs_Sept/real225_nicherate0.02_dispersal1.5_200gens/raster/speciesEnd_at_", summary.time, ".asc", sep=""))
    main_header <- paste("Species endemism after", summary.time, "steps")
    plot(speciesEndemism.sum, main=main_header)

    abundanceEndemism.sum   <- raster(paste("C:/Work/Simulation/test_runs_Sept/real225_nicherate0.02_dispersal1.5_200gens/raster/abundanceEnd_at_", summary.time, ".asc", sep=""))
    main_header <- paste("log(Abundance endemism) after", summary.time, "steps")
    plot(log(abundanceEndemism.sum), main=main_header)
  }

}