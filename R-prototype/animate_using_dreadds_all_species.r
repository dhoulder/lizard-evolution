rm(list=ls())

library(dreadds)
library(raster)
library(data.table)

source("~/code/DREaD_ds/R-prototype/dynamicDisplay.r")

# turn the input arguments into local variables with relevant names
input.args	<- commandArgs(trailingOnly = TRUE)
output.dir	<- input.args[1]
config.file <- input.args[2]
#dispersal	 <- as.numeric(input.args[2])
cat("\nExternal arguments received in R\n\toutput.dir\t", output.dir, "\n\tconfig.file\t", config.file, "\n")

rasterFromEnv <- function(dreadModel, envTemplate, envNumber=1) {
  env         <- dreadModel$getEnv()[[envNumber]]
  env.extent  <- extent(envTemplate)
  env.ras     <- raster(env.extent, nrows=envTemplate@nrows, ncols=envTemplate@ncols)
  env.ras[]   <- env[]
  return(env.ras)
}

getAllDemes   <- function(sp.df, model) {
  # first add a species name to each species' deme table
  demesDataFrames <- model$getDemes()
  for (species_name in rownames(sp.df)) {
    species_name_column <- rep(species_name, nrow(demesDataFrames[[species_name]]))
    demesDataFrames[[species_name]] <- cbind(species_name=species_name_column, demesDataFrames[[species_name]])
  }
  allDemes <-rbindlist(demesDataFrames)
  return(allDemes)
}

randomTime   <- function(max) {
  # this function aims to introduce a random number of steps to ensure that each run is different
  x <- round(runif(1, 0, max))
  for (i in 0:x) {
	y <- i * 2
  }
  cat("\n", x, "random steps completed\n")
}

envOrig <- raster("~/simulation/input_data/AMT2.5min-grids/T5000/Dan AMT2.5min_bio PaleoClimDat - MinTemp.txt.tif")
environment.rows <- nrow(envOrig)

niche.evolution.rate    <- 0.02  # this should be extracted from the model arguments, but not yet possible
dispersal               <- 2     # this should be extracted from the model arguments, but not yet possible
gene.flow.max.distance  <- 20    # this should be extracted from the model arguments, but not yet possible
genome.column.name.format   <- "genetic_position_"
paleoTime.start         <- 5000
paleoTime.step          <- 2.5

# set up the plotting environment
# display settings
do.display            <- TRUE
do.display.diff       <- TRUE
do.display.genome     <- TRUE
do.text.output        <- TRUE
image_to_file         <- TRUE
raster_to_file        <- FALSE
image_counter		  <- 1  # give images consecutive integer names to work with ffmpeg movies
image_frequency   <- 1

run_interval      <- 1  # this is the number of steps to run before returning to R
nSteps 				  <- 2000

# define and create the directory structure
base.dir              <- "~/simulation/multisize_test/"
#output.dir            <- paste(base.dir, "output1/", sep='')
image.dir             <- paste(output.dir, "images/", sep='')

if (dir.exists(image.dir) == FALSE) {
  if (dir.exists(output.dir) == FALSE) {
    dir.create(output.dir)
    cat("\nCreated directory", output.dir, "\n")
  }
  
  dir.create(image.dir)
  cat("\nCreated directory", image.dir, "\n")
}

if (do.display) {
  if (do.display.diff) {
    my.display <- display.initialise.colours()
    my.colours      <- my.display[[1]]
    my.coloursdiff  <- my.display[[2]]
  } else {
    my.display <- display.initialise.colours()
    my.colours <- my.display[[1]]
  }
}

starttime_global <- Sys.time()

m <- createDreadDS(
        config.file = config.file,
        output.dir = output.dir,
        iterations = nSteps
      )

# run a process of variable length here to (hopefully) introduce randomness
  set.seed(Sys.time())
  randomTime(100000)
	  
# step through the model 
currentStep <- 0
while (currentStep <= nSteps) {
  
  starttime_timestep <- Sys.time()
  
  m$runSteps(run_interval)
  currentStep <- currentStep + run_interval

  ##### species summary and plots ######
  
  # get the species data frame and number of species
  sp.df <- m$getSpecies()
  sp_all_count <- nrow(sp.df)  # all species ever (ie count of branches)

  # subset sp.df to only extant species
  sp.df <- sp.df[which(sp.df$extinction == -1 & sp.df$split == -1), ]
  sp_count <- nrow(sp.df)
  
  # get all demes
  all.demes <- getAllDemes(sp.df, model=m)
 
  # on first step, identify the genome position columns
  if (exists("genome.columns")==FALSE) {
    genome.columns <- grep(genome.column.name.format, names(all.demes), value=FALSE)
  }  
  
  # record the time before present for map and text outputs
  paleoTime <- paleoTime.start - (paleoTime.step * currentStep)

  for (sp_idx in 1:sp_count) {
    sp_current <- sp.df[sp_idx,]
      
    # calculate and print a species summary - probably drop this once development is done
    if (do.text.output) {
      sp.summary <- data.frame(speciesID = sp_current$id,
                              speciesCount = sp_count,
                              range = sp_current$cell_count,
                              total_amount = sp_current$population,
                              niche0.position.mean = sp_current$niche_position_mean_0,
                              niche0.position.sd = sp_current$niche_position_sd_0,
                              niche0.breadth.mean = sp_current$niche_breadth_mean_0,
                              niche0.breadth.sd = sp_current$niche_breadth_sd_0,
                              niche0.sp.min = sp_current$niche_min_0,
                              niche0.sp.max = sp_current$niche_max_0)
    
      sp.summary <- round(sp.summary, 4)

      list.for.text <- c(list(current.time=paleoTime,
                              elapsed.time.total = difftime(Sys.time(), starttime_global),
                              elapsed.time.step = difftime(Sys.time(), starttime_timestep)),
                              as.list(sp.summary[1,]))
      text.update(list(species_range_niche=list.for.text))
    }
  }

  # update the dynamic plot FOR ALL SPECIES at specified frequency
  if (do.display) {
    if (currentStep %% image_frequency == 0) {  # only do the image at the specified frequency
      
      # get the data needed for plots

      env.ras <- rasterFromEnv(m, envOrig, 1)
      
      if (image_to_file) {
        display.to.file.start(image_dir=image.dir, time=image_counter, image_filename = paste("animation_multisp_", sep=""))
		# use time=image_counter for consecutively numbered images (eg for ffmpeg) or time=currentStep to number by the model time
      }

      model.params <- list(niche.breadth=round(sp.df$niche_breadth_mean_0[1]), niche.evolution.rate=niche.evolution.rate, dispersal = dispersal)

      # display.update.multispecies(list(current.time = paleoTime, model.params = model.params), env = env.ras, allDemes = all.demes, plot_demes_amount_position = TRUE)
	  
      display.update.multispecies(list(current.time = paleoTime, model.params = model.params), env = env.ras, allDemes = all.demes, plot_species_ranges = TRUE)	  

      display.update.multispecies(list(current.time = paleoTime, model.params = model.params, speciesCount = sp_count), env = env.ras, allDemes = all.demes, plot_richness = TRUE)

      if (do.display.genome) {
        display.update.multispecies(list(current.time = paleoTime, model.params = model.params, genome.columns=genome.columns),
                                    env = env.ras, 
                                    allDemes = all.demes,
                                    plot_genome_scatter = TRUE,
                                    plot_genome_map = TRUE)
      }
 
      if (image_to_file) {
        display.to.file.stop()
		image_counter <- image_counter + 1
      }
    }
  }
}
