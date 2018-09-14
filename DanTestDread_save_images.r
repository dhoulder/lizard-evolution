rm(list=ls())
gc()

scripts <- c("disperse_evolve_compete_ds.r","seedSpecies.R","environmentalChange.R", "DREaD_ds.R",
             "generateSummaryStatistics.R", "helperFunctions.R", "dataStructures.r", "dynamicDisplay.r")
lapply(scripts, source)

run.name              <- "anim_real225_nicherate0.02_dispersal3_300gens_v5"

#sample parameters
total.time            <- 300
dispersal             <- 3              # dispersal distance
niche.evolution.rate  <- 0.02
env.amp               <- 0 #runif(1, 0.25, 2)
env.freq              <- runif(1, 10, 25)
NEb                   <- runif(1, 0.0025, 1)
niche.blocksize       <- 0.05
suitability.mode      <- "sine"
speciation.gene.distance <- 50  # this parameter will need to be set with the drift rate
environment.source    <- "~/code/DREaD_extras/realAlps225.asc"  # 'internal to generate in the code
# or a raster file to load
#environment.source    <- "internal"
environment.dimension <- 225

initial.breadth       <- 3
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
do.animate            <- FALSE
image_to_file        <- TRUE

# create the image output directory
output_dir            <- paste("/short/ka2/dfr805/simulation/test_runs/", run.name, sep="")
image_dir             <- paste(output_dir, "/images/", sep="")

if (dir.exists(image_dir) == FALSE) {

  if (dir.exists(output_dir) == FALSE) {
    dir.create(output_dir)
    cat("\nCreated directory", output_dir, "\n")
  }

  dir.create(image_dir)
  cat("\nCreated directory", image_dir, "\n")
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

# ADD CODE HERE TO TURN THE IMAGES INTO AN ANIMATION