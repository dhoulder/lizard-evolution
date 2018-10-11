# this file combines and summarises output rasters from the simulations
library(raster)

model_name[1]    <- "real225_nicherate0.02_dispersal1.5_200gens_v1"
model_name[2]    <- "real225_nicherate0.02_dispersal1.5_200gens_v2"

output_image

summary.time  <- 50
layer.count   <- 100

output_image     <- paste("C:/Work/Simulation/test_runs_Sept/real225_nicherate0.02_dispersal1.5_", summary.time, ".pdf", sep="")

pdf(file=output_image, width=8.27, height=10, onefile = T)

for (name in model_name) {

  raster_dir      <- paste("C:/Work/Simulation/test_runs_Sept/", name, "/raster/", sep="")
  abund.ras       <- raster(paste(raster_dir, "abundance_at_", summary.time ,".asc", sep=""))
  richness.ras    <- raster(paste(raster_dir, "richness_at_", summary.time ,".asc", sep=""))
  speciesEnd.ras  <- raster(paste(raster_dir, "speciesEnd_at_", summary.time ,".asc", sep=""))
  abundEnd.ras    <- raster(paste(raster_dir, "abundanceEnd_at_", summary.time ,".asc", sep=""))

  #windows(10,10)
  par(mfcol=c(2,2))

  main_header <- paste("Species richness after", summary.time, "steps")
  plot(richness.ras, main=main_header)

  main_header <- paste("Abundance after", summary.time, "steps from", layer.count, "runs")
  plot(abund.ras, main=main_header)

  main_header <- paste("Species endemism after", summary.time, "steps")
  plot(speciesEnd.ras, main=main_header)

  main_header <- paste("log(Abundance endemism) after", summary.time, "steps")
  plot(log(abundEnd.ras), main=main_header)

  mtext(name, side = 3, line = -1, outer = TRUE, cex=0.6)

}

# and now look at the relationship between separate sets of runs
raster_dir    <- paste("C:/Work/Simulation/test_runs_Sept/", model_name[1], "/raster/", sep="")
abund1.ras       <- raster(paste(raster_dir, "abundance_at_", summary.time ,".asc", sep=""))
richness1.ras    <- raster(paste(raster_dir, "richness_at_", summary.time ,".asc", sep=""))
speciesEnd1.ras  <- raster(paste(raster_dir, "speciesEnd_at_", summary.time ,".asc", sep=""))
abundEnd1.ras    <- raster(paste(raster_dir, "abundanceEnd_at_", summary.time ,".asc", sep=""))

raster_dir    <- paste("C:/Work/Simulation/test_runs_Sept/", model_name[2], "/raster/", sep="")
abund2.ras       <- raster(paste(raster_dir, "abundance_at_", summary.time ,".asc", sep=""))
richness2.ras    <- raster(paste(raster_dir, "richness_at_", summary.time ,".asc", sep=""))
speciesEnd2.ras  <- raster(paste(raster_dir, "speciesEnd_at_", summary.time ,".asc", sep=""))
abundEnd2.ras    <- raster(paste(raster_dir, "abundanceEnd_at_", summary.time ,".asc", sep=""))

#windows(10,10)
par(mfcol=c(2,2))

main_header <- paste("Sp. richness @", summary.time, "steps Diff between runs")
plot((richness1.ras - richness2.ras), main=main_header)

main_header <- paste("Sp. richness @", summary.time, "steps cells run1 v run2")
plot(richness1.ras[], richness2.ras[], main=main_header, xlab="richness run 1", ylab="richness run 2", col="blue")

# text for r squared
rich.lm <- lm(richness1.ras[] ~ richness2.ras[])
r2 <- round(summary(rich.lm)$r.squared, 3)
text(50, 8, paste("r-squared:", r2), pos=4, font=2)

main_header <- paste("Abundance @", summary.time, "steps Diff between runs")
plot((abund1.ras - abund2.ras), main=main_header)

main_header <- paste("Abundance @", summary.time, "steps cells run1 v run2")
plot(abund1.ras[], abund2.ras[], main=main_header, xlab="abundance run 1", ylab="abundance run 2", col="blue")

# text for r squared
abund.lm <- lm(abund1.ras[] ~ abund2.ras[])
r2 <- round(summary(abund.lm)$r.squared, 3)
text(18, 3, paste("r-squared:", r2), pos=4, font=2)

# now compare endemism results
#windows(10,10)
par(mfcol=c(2,2))

main_header <- paste("Sp. endemism @", summary.time, "steps Diff between runs")
plot((speciesEnd1.ras - speciesEnd2.ras), main=main_header)

main_header <- paste("Sp. endemism @", summary.time, "steps   run1 v run2")
plot(speciesEnd1.ras[], speciesEnd2.ras[], main=main_header, xlab="Species endemism run 1", ylab="Species endemism run 2", col="blue")

# text for r squared
sp_end.lm <- lm(speciesEnd1.ras[] ~ speciesEnd2.ras[])
r2 <- round(summary(sp_end.lm)$r.squared, 3)
text(0.0005, 0.0003, paste("r-squared:", r2), pos=4, font=2)

main_header <- paste("Abund. endemism @", summary.time, "steps Diff between runs")
plot((log(abundEnd1.ras) - log(abundEnd2.ras)), main=main_header)

main_header <- paste("Abund. endemism @", summary.time, "steps   run1 v run2")
plot(log(abundEnd1.ras[]), log(abundEnd2.ras[]), main=main_header, xlab="log(abund. endemism run 1)", ylab="log(abund. endemism run 2)", col="blue")

# text for r squared
abundEnd.lm <- lm(abundEnd1.ras[] ~ abundEnd2.ras[])
r2 <- round(summary(abundEnd.lm)$r.squared, 3)
text(-12, -14, paste("r-squared:", r2), pos=4, font=2)

dev.off()
