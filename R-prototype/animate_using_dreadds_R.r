rm(list=ls())

library(dreadds)
library(raster)
library(data.table)

source("~/code/DREaD_ds/R-prototype/dynamicDisplay.r")

rasterFromEnv <- function(dreadModel, envTemplate, envNumber=1) {
  env         <- dreadModel$getEnv()[[envNumber]]
  env.extent  <- extent(envTemplate)
  env.ras     <- raster(env.extent, nrows=envTemplate@nrows, ncols=envTemplate@ncols)
  env.ras[]   <- env[]
  return(env.ras)
}

envOrig <- raster("/home/danr/code/DREaD_extras/habitat_sizes.asc")
environment.rows <- nrow(envOrig)

niche.evolution.rate    <- 0  # this should be extracted from the model arguments, but not yet possible
dispersal               <- 2  # this should be extracted from the model arguments, but not yet possible
gene.flow.max.distance  <- 20  # this should be extracted from the model arguments, but not yet possible
genome.column.name.format   <- "genetic_position_"

# set up the plotting environment
# display settings
do.display            <- TRUE
do.display.diff       <- TRUE
do.display.genome     <- TRUE
do.text.output        <- TRUE
image_to_file         <- TRUE
raster_to_file        <- FALSE
image_frequency       <- 20

run_interval          <- 5  # this is the number of steps to run before returning to R

# define and create the directory structure
base.dir              <- "~/simulations/multisize_test/"
output.dir            <- paste(base.dir, "output/", sep='')
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

# set up and start running the model here
nSteps <- 1000

starttime_global <- Sys.time()

m <- createDreadDS(
        config.file = "/home/danr/simulations/multisize_test/multisize.conf",
        output.dir = output.dir,
        iterations = nSteps
      )

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
if(sp_all_count > 1){browser()}
  # subset sp.df to only extant species
  sp.df <- sp.df[which(sp.df$extinction == -1 & sp.df$split == -1), ]
  sp_count <- nrow(sp.df)
  
  # get all demes
  all.demes <- rbindlist(m$getDemes())
  # BUT NEED TO DISTINGUISH SPECIES BEFORE COMBINING ALL OF THE DEMES

  for (sp_idx in 1:sp_count) {
    sp_current <- sp.df[sp_idx,]

    species_name <- paste("species_", sp_current$id, sep="")
    sp_current_demes <- m$getDemes()[[species_name]]
    setDT(sp_current_demes) #make sp_current_demes into a data.table for compatibility with existing code in dynamicDisplay.r
    
    # on first step, identify the genome position columns
    if (exists("genome.columns")==FALSE) {
      demes.df <- m$getDemes()[[species_name]]
      genome.columns        <- grep(genome.column.name.format, names(demes.df), value=FALSE)
      rm(demes.df)
    }    
   
    # calculate and print a species summary - probably drop this once development is done
    sp.summary <- data.frame(speciesID = sp_current$id,
                              speciesCount = sp_count,
                              range = sp_current$cell_count,
                              total_amount = sp_current$population,
                              niche0.position.mean = sp_current$niche_position_mean_0,
                              niche0.position.sd = sp_current$niche_position_sd_0,
                              niche0.breadth.mean = sp_current$niche_breadth_mean_0,
                              niche0.breadth.sd = sp_current$niche_breadth_sd_0,
                              niche0.sp.min = sp_current$niche_min_0,
                              niche0.sp.max = sp_current$niche_max_0
    )
    
    sp.summary <- round(sp.summary, 4)
    if (do.text.output) {
      list.for.text <- c(list(current.time=currentStep,
                              elapsed.time.total = difftime(Sys.time(), starttime_global),
                              elapsed.time.step = difftime(Sys.time(), starttime_timestep)),
                         as.list(sp.summary[1,]))
      text.update(list(species_range_niche=list.for.text))
    }

    # update the dynamic plot at specified frequency
    if (do.display) {
      if (currentStep %% image_frequency == 0) {  # only do the image at the specified frequency
        
        # get the data needed for plots
        sp_current_demes <- m$getDemes()[[species_name]]
        setDT(sp_current_demes) #make sp_current_demes into a data.table for compatibility with existing code in dynamicDisplay.r
        
        env.ras <- rasterFromEnv(m, envOrig, 1)
        
        if (image_to_file) {
          display.to.file.start(image.dir, currentStep, image_filename = paste("animation_sp", sp_current$id, "_", sep=""))
        }

        model.params <- list(niche.breadth=round(sp.summary$niche0.breadth.mean), niche.evolution.rate=niche.evolution.rate, dispersal = dispersal)
  
        display.update(list(env = env.ras, demes_amount_position=sp_current_demes, current.time=currentStep, model.params = model.params))
   
        if (do.display.diff) {
          display.update(list(env = env.ras, demes_amount_position_diff = sp_current_demes, current.time = currentStep))
        }
    
        if (do.display.genome) {
          display.update(list(env = env.ras,
                              demes_genecolour = sp_current_demes,
                              current.time = currentStep,
                              genome.columns = genome.columns))
        }
   
        if (image_to_file) {
          display.to.file.stop()
        }
        
        gc()  # trying to solve a crash after 100 - 200 steps
      }
    }
  }
}
