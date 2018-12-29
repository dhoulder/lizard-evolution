rm(list=ls())
gc()

scripts <- c("disperse_evolve_compete_ds.r", "environmentalChange.r", "DREaD_ds.r", "helperFunctions.R", "dataStructures.r", "dynamicDisplay.r")
lapply(scripts, source)

# turn the input arguments into local variables with relevant names
# turn the input arguments into local variables with relevant names
input.args	          <- commandArgs(trailingOnly = TRUE)
dispersal_dist        <- as.numeric(input.args[1])
timesteps        			<- as.numeric(input.args[2])
output_dir						<- input.args[3]

input_prexix					<- "out"

#run.name              <- "anim_real225_nicherate0.02_dispersal3_300gens_v5"

#sample parameters
total.time            <- timesteps
dispersal             <- dispersal_dist              # dispersal distance in cells

# MOST OF THE FOLLOWING ARGUMENTS SHOULD BE LOADED FROM THE c++ OUTPUT FILES

niche.evolution.rate  <- 0.02
env.amp               <- 0   #runif(1, 0.25, 2)
env.freq              <- 50  #runif(1, 10, 25)
breadth.evolution.rate  <- 0
suitability.mode      <- "sine"
speciation.gene.distance <- 50  # this parameter will need to be set with the drift rate
environment.source    <- "/short/ka2/simulation/input_data/realAlps250_rescaled_100_bio01.asc"  # 'internal to generate in the code
# or a raster file to load
#environment.source    <- "internal"
environment.dimension <- 250

initial.breadth       <- 20
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
input.dir             <- output_dir
output.dir            <- input.dir  #paste("/short/ka2/dfr805/simulation/test_runs/", run.name, sep="")
image_dir             <- paste(output.dir, "/images/", sep="")

if (dir.exists(image_dir) == FALSE) {

  if (dir.exists(output.dir) == FALSE) {
    dir.create(output.dir)
    cat("\nCreated directory", output.dir, "\n")
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
                  input.dir = input.dir, input.prefix = input_prexix)

# ADD CODE HERE TO TURN THE IMAGES INTO AN ANIMATION