rm(list=ls())
gc()

scripts <- c("disperse_evolve_compete_ds.r","seedSpecies.R","environmentalChange.R", "DREaD_ds.R",
             "generateSummaryStatistics.R", "helperFunctions.R", "dataStructures.r", "dynamicDisplay.r")
lapply(scripts, source)

#sample parameters
total.time            <- 100
dispersal             <- 2                    # dispersal distance
niche.evolution.rate  <- 0.1
env.amp               <- 0 #runif(1, 0.25, 2)
env.freq              <- runif(1, 10, 25)
NEb                   <- runif(1, 0.0025, 1)
niche.blocksize       <- 0.05
suitability.mode      <- "sine"
speciation.gene.distance <- 10  # this parameter will need to be set with the drift rate
do.display            <- TRUE
do.display.diff       <- TRUE
do.text.output        <- TRUE
environment.source    <- "~/Work/Software/dan-github/DREaD_extras/circular.asc"  # 'internal to generate in the code
                          # or a raster file to load
#environment.source    <- "internal"

if (do.display) {
  if (do.display.diff) {
    my.colours.double <- display.initialise.double()
    my.colours      <- my.colours.double[[1]]
    my.coloursdiff  <- my.colours.double[[2]]
    rm(my.colours.double)
  } else {
    my.colours <- display.initialise()
  }
}

# Run model
simulation.1 <- DREaD_ds(total.time=total.time, dispersal=dispersal, amp=env.amp, freq=env.freq,
                  niche.evolution.rate = niche.evolution.rate, breadth.ev.rate=NEb, enviro.hetero=T,
                  enviro.mode="sine", suitability.mode=suitability.mode,
                  speciation.gene.distance=speciation.gene.distance, environment.source=environment.source)
