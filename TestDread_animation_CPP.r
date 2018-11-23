rm(list=ls())
gc()

scripts <- c("disperse_evolve_compete_ds.r", "environmentalChange.R", "DREaD_ds.R", "helperFunctions.R", "dataStructures.r", "dynamicDisplay.r")
lapply(scripts, source)

run.name              <- "anim_real225_nicherate0.02_dispersal3_300gens_v5"

#sample parameters
total.time            <- 500
dispersal             <- 4              # dispersal distance
niche.evolution.rate  <- 0.02
env.amp               <- 1   #runif(1, 0.25, 2)
env.freq              <- 50  #runif(1, 10, 25)
breadth.evolution.rate  <- 0
suitability.mode      <- "sine"
speciation.gene.distance <- 50  # this parameter will need to be set with the drift rate
environment.source    <- "E:/Work/Software/dan-github/DREaD_extras/realAlps225.asc"  # 'internal to generate in the code
# or a raster file to load
#environment.source    <- "internal"
environment.dimension <- 225

initial.breadth       <- 4
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
image_to_file         <- TRUE
raster_to_file        <- FALSE
generations_to_save     <- c(0)

# input and output directories
input.dir             <- "E:/Work/Simulation/test_runs_from_CPP/test_output9Nov_disp4v1"
output_dir            <- input.dir  #paste("/short/ka2/dfr805/simulation/test_runs/", run.name, sep="")
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
simulation.1 <- DREaD_read_plot(total.time = total.time, dispersal = dispersal, amp = env.amp, freq = env.freq,
                  niche.evolution.rate = niche.evolution.rate, breadth.ev.rate = breadth.evolution.rate,
                  enviro.mode = "sine", suitability.mode = suitability.mode,
                  speciation.gene.distance = speciation.gene.distance, environment.source = environment.source,
                  initial.species.defined = initial.species.defined,
                  input.dir = input.dir, input.prefix = '')

# ADD CODE HERE TO TURN THE IMAGES INTO AN ANIMATION