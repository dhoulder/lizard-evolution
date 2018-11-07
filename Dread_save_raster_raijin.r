rm(list=ls())
gc()

setwd("~/code/DREaD_ds")

scripts <- c("disperse_evolve_compete_ds.r","seedSpecies.R","environmentalChange.R", "DREaD_ds.R",
             "generateSummaryStatistics.R", "helperFunctions.R", "dataStructures.r", "dynamicDisplay.r", "output.r")
lapply(scripts, source)

# turn the input arguments into local variables with relevant names
input.args			<- commandArgs(trailingOnly = TRUE)
dispersal_dist		<- as.numeric(input.args[1])
timesteps			<- as.numeric(input.args[2])
run_number			<- input.args[3]  # not passed to the function - used as a global variable - should fix this
output_dir			<- input.args[4]


#run.name              <- "real225_nicherate0.02_disp2_amp1_freq40_200gens_1r"

#parameters
total.time            <- timesteps
dispersal             <- dispersal_dist			# dispersal distance
niche.evolution.rate  <- 0.02
env.amp               <- 0  #runif(1, 0.25, 2)
env.freq              <- 50 #runif(1, 10, 25)
NEb                   <- runif(1, 0.0025, 1)
niche.blocksize       <- 0.05
suitability.mode      <- "sine"
speciation.gene.distance <- 50  # this parameter will need to be set with the drift rate

#environment.source    <- "internal"
# or a raster file to load
environment.source    <- "~/code/DREaD_extras/realAlps225.asc"  # 'internal to generate in the code

environment.dimension <- 225

initial.breadth       <- 3
initial.cell          <- -1
initial.extent        <- NA
initial.species.defined = list(initial.breadth = initial.breadth,
                               initial.cell = initial.cell,
                               initial.extent = initial.extent)

# display and output settings
do.display              <- FALSE
do.display.diff         <- FALSE
do.display.genome       <- FALSE
do.text.output          <- TRUE
do.animate              <- FALSE
image_to_file           <- FALSE
raster_to_file          <- TRUE
generations_to_save     <- c(1, 10, 50, 100, 150, 200)
run_number_in_filename  <- TRUE	# not passed to the function - used as a global variable - should fix this

# create the raster output directory
#output_dir            <- paste("/short/ka2/dfr805/simulation/test_runs/", run.name, sep="")
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
}

# Run model
simulation.1 <- DREaD_ds(total.time = total.time, dispersal = dispersal, amp = env.amp, freq = env.freq,
                         niche.evolution.rate = niche.evolution.rate, breadth.ev.rate = NEb, enviro.hetero = FALSE,
                         enviro.mode = "sine", suitability.mode = suitability.mode,
                         speciation.gene.distance = speciation.gene.distance, environment.source = environment.source,
                         initial.species.defined = initial.species.defined)
