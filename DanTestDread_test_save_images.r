rm(list=ls())
gc()

scripts <- c("disperse_evolve_compete_ds.r","seedSpecies.R","environmentalChange.R", "DREaD_ds.R",
             "generateSummaryStatistics.R", "helperFunctions.R", "dataStructures.r", "dynamicDisplay.r")
lapply(scripts, source)

library(animation)

#sample parameters
total.time            <- 3
dispersal             <- 2              # dispersal distance
niche.evolution.rate  <- 0.04
env.amp               <- 0 #runif(1, 0.25, 2)
env.freq              <- runif(1, 10, 25)
NEb                   <- runif(1, 0.0025, 1)
niche.blocksize       <- 0.05
suitability.mode      <- "sine"
speciation.gene.distance <- 50  # this parameter will need to be set with the drift rate
#environment.source    <- "~/Work/Software/dan-github/DREaD_extras/circular.asc"  # 'internal to generate in the code
# or a raster file to load
environment.source    <- "internal"

initial.breadth       <- -1
initial.cell          <- -1
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
output_to_file        <- TRUE
image_dir             <- "~/Dropbox/Simulation/test_runs_July2018/anim_10gens_nicherate0.02_dispersal1/images/"

if (do.display) {
  if (do.display.diff) {
    my.display <- display.initialise.2by2(output_to_file)
    my.colours      <- my.display[[1]]
    my.coloursdiff  <- my.display[[2]]
    my.par          <- my.display[[3]]
  } else {
    my.display <- display.initialise()
    my.colours <- my.display[[1]]
  }

  if (do.animate) {
    ani.record(reset = TRUE) # clear animation history before recording
    #display.to.file.start(image_dir, 0, image_filename = "animation")
  }

}

# Run model
simulation.1 <- DREaD_ds(total.time = total.time, dispersal = dispersal, amp = env.amp, freq = env.freq,
                  niche.evolution.rate = niche.evolution.rate, breadth.ev.rate = NEb,  enviro.hetero = T,
                  enviro.mode = "sine", suitability.mode = suitability.mode,
                  speciation.gene.distance = speciation.gene.distance, environment.source = environment.source,
                  initial.species.defined = initial.species.defined)

ani.replay()

ani.options('ani.width'=960, 'ani.height'=960, 'interval'=0.5, 'ani.dev'='png')
setwd("~/Dropbox/Simulation/test_runs_July2018/anim_10gens_nicherate0.02_dispersal1")

saveHTML(ani.replay(), img.name="animation", htmlfile = "anim_10gens_nicherate0.02_dispersal1.html", navigator = FALSE)

#display.to.file.stop()
