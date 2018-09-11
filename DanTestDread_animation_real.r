rm(list=ls())
gc()

scripts <- c("disperse_evolve_compete_ds.r","seedSpecies.R","environmentalChange.R", "DREaD_ds.R",
             "generateSummaryStatistics.R", "helperFunctions.R", "dataStructures.r", "dynamicDisplay.r")
lapply(scripts, source)

library(animation)

run.name              <- "anim_real225_nicherate0.02_dispersal1.5_250gens_v1"

#sample parameters
total.time            <- 250
dispersal             <- 1.5              # dispersal distance
niche.evolution.rate  <- 0.02
env.amp               <- 0 #runif(1, 0.25, 2)
env.freq              <- runif(1, 10, 25)
NEb                   <- runif(1, 0.0025, 1)
niche.blocksize       <- 0.05
suitability.mode      <- "sine"
speciation.gene.distance <- 50  # this parameter will need to be set with the drift rate
environment.source    <- "~/Work/Software/dan-github/DREaD_extras/realAlps225.asc"  # 'internal to generate in the code
# or a raster file to load
#environment.source    <- "internal"
environment.dimension <- 225

initial.breadth       <- 3
#initial.cell          <- 6121
initial.cell          <- -1
#initial.extent        <- c(8.18, 37.98, 22.09, 55)
initial.extent        <- NA
initial.species.defined = list(initial.breadth = initial.breadth,
                               initial.cell = initial.cell,
                               initial.extent = initial.extent)

# display settings
do.display            <- TRUE
do.display.diff       <- TRUE
do.display.genome     <- TRUE
do.text.output        <- TRUE
do.animate            <- TRUE
output_to_file        <- FALSE

if (do.display) {
  if (do.display.diff) {
    my.colours.double <- display.initialise.2by2()
    my.colours      <- my.colours.double[[1]]
    my.coloursdiff  <- my.colours.double[[2]]
    rm(my.colours.double)
  } else {
    my.colours <- display.initialise()
  }

  if (do.animate) {
    ani.record(reset = TRUE) # clear animation history before recording
  }

}

# Run model
simulation.1 <- DREaD_ds(total.time = total.time, dispersal = dispersal, amp = env.amp, freq = env.freq,
                  niche.evolution.rate = niche.evolution.rate, breadth.ev.rate = NEb,  enviro.hetero = T,
                  enviro.mode = "sine", suitability.mode = suitability.mode,
                  speciation.gene.distance = speciation.gene.distance, environment.source = environment.source,
                  initial.species.defined = initial.species.defined, environment.dimension=environment.dimension)

gc()

# save the animation
#ani.replay()

ani.options('ani.width'=900, 'ani.height'=900, 'interval'=0.4)

save.dir <- paste("D:/Simulation/test_runs/", run.name, sep="")

if (dir.exists(save.dir) == FALSE) {
  dir.create(save.dir)
}

setwd(save.dir)

saveHTML(ani.replay(), img.name="animation", htmlfile = "index.html", navigator = FALSE)