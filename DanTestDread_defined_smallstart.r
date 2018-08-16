rm(list=ls())
gc()

scripts <- c("disperse_evolve_compete_ds.r","seedSpecies.R","environmentalChange.R", "DREaD_ds.R",
             "generateSummaryStatistics.R", "helperFunctions.R", "dataStructures.r", "dynamicDisplay.r")
lapply(scripts, source)

#sample parameters
total.time            <- 150
dispersal             <- 1.5              # dispersal distance
niche.evolution.rate  <- 0.015
env.amp               <- 0 #runif(1, 0.25, 2)
env.freq              <- runif(1, 10, 25)
NEb                   <- runif(1, 0.0025, 1)
niche.blocksize       <- 0.05
suitability.mode      <- "sine"
speciation.gene.distance <- 25  # this parameter will need to be set with the drift rate
do.display            <- TRUE
do.display.diff       <- TRUE
do.display.genome     <- TRUE
do.text.output        <- TRUE
environment.source    <- "~/Work/Software/dan-github/DREaD_extras/circular.asc"  # 'internal to generate in the code
                          # or a raster file to load
#environment.source    <- "internal"
initial.breadth       <- 4
initial.cell          <- 6121
initial.extent        <- c(8.18, 37.98, 22.09, 55)
initial.species.defined = list(initial.breadth = initial.breadth,
                               initial.cell = initial.cell,
                               initial.extent = initial.extent)

if (do.display) {
  if (do.display.diff) {
    my.colours.double <- display.initialise.2by2()
    my.colours      <- my.colours.double[[1]]
    my.coloursdiff  <- my.colours.double[[2]]
    rm(my.colours.double)
  } else {
    my.colours <- display.initialise()
  }
}

# Run model
simulation.1 <- DREaD_ds(total.time = total.time, dispersal = dispersal, amp = env.amp, freq = env.freq,
                  niche.evolution.rate = niche.evolution.rate, breadth.ev.rate = NEb,  enviro.hetero = T,
                  enviro.mode = "sine", suitability.mode = suitability.mode,
                  speciation.gene.distance = speciation.gene.distance, environment.source = environment.source,
                  initial.species.defined = initial.species.defined)
