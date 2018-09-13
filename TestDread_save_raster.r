rm(list=ls())
gc()

scripts <- c("disperse_evolve_compete_ds.r","seedSpecies.R","environmentalChange.R", "DREaD_ds.R",
             "generateSummaryStatistics.R", "helperFunctions.R", "dataStructures.r", "dynamicDisplay.r", "output.r")
lapply(scripts, source)

run.name              <- "real225_nicherate0.02_dispersal1.5_250gens"

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
initial.cell          <- -1
initial.extent        <- NA
initial.species.defined = list(initial.breadth = initial.breadth,
                               initial.cell = initial.cell,
                               initial.extent = initial.extent)

# display and output settings
do.display            <- FALSE
do.display.diff       <- FALSE
do.display.genome     <- FALSE
do.text.output        <- TRUE
do.animate            <- FALSE
image_to_file         <- FALSE
raster_to_file        <- TRUE
generations_to_save   <- c(50, 100)

# create the raster output directory
output_dir            <- paste("C:/Work/Simulation/test_runs_Sept/", run.name, sep="")
raster_dir             <- paste(output_dir, "/raster/", sep="")

if (dir.exists(raster_dir) == FALSE) {

  if (dir.exists(output_dir) == FALSE) {
    dir.create(output_dir)
    cat("\nCreated directory", output_dir, "\n")
  }

  dir.create(raster_dir)
  cat("\nCreated directory", raster_dir, "\n")
}

if (do.display) {
  if (do.display.diff) {
    my.display <- display.initialise.2by2(image_to_file)
    my.colours      <- my.display[[1]]
    my.coloursdiff  <- my.display[[2]]
    my.par          <- my.display[[3]]
  } else {
    my.display <- display.initialise()
    my.colours <- my.display[[1]]
  }

  # if (do.animate) {
  #   ani.record(reset = TRUE) # clear animation history before recording
  #   #display.to.file.start(image_dir, 0, image_filename = "animation")
  # }

}

# Run model
simulation.1 <- DREaD_ds(total.time = total.time, dispersal = dispersal, amp = env.amp, freq = env.freq,
                  niche.evolution.rate = niche.evolution.rate, breadth.ev.rate = NEb,  enviro.hetero = T,
                  enviro.mode = "sine", suitability.mode = suitability.mode,
                  speciation.gene.distance = speciation.gene.distance, environment.source = environment.source,
                  initial.species.defined = initial.species.defined)
