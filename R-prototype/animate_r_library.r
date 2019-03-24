rm(list=ls())

library(dreadds)
library(raster)

source("~/code/DREaD_ds/R-prototype/dynamicDisplay.r")

rasterFromEnv <- function(dreadModel, envTemplate, envNumber=1) {
  env       <- dreadModel$get_env()[[envNumber]]
  env.extent  <- extent(envTemplate)
  return(raster(env.extent, nrows=envTemplate@nrows, ncols=envTemplate@ncols))
}

envOrig <- raster("/home/danr/code/DREaD_extras/habitat_sizes.asc")
environment.rows <- nrow(envOrig)

niche.evolution.rate <- 0  # this should be extracted from the model arguments, but not yet possible
dispersal            <- 2  # this should be extracted from the model arguments, but not yet possible
# set up the plotting environment
# display settings
do.display            <- TRUE
do.display.diff       <- TRUE
do.display.genome     <- TRUE
do.text.output        <- TRUE
image_to_file         <- TRUE
raster_to_file        <- FALSE

# define and create the directory structure
base.dir              <- "~/simulations/multisize_test/"
output.dir            <- paste(base.dir, "output/", sep='')
image_dir             <- paste(output.dir, "images", sep='')

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

m <- dreadds_model(
        "/home/danr/simulations/multisize_test/multisize.conf",
        c('-o', output.dir)
      )

# step through the model 
for (currentTime in 1:nSteps) {
  
  starttime_timestep <- Sys.time()
  
  m$do_step()
  
  ##### species summary and plots ######
  
  # get the species data frame and number of species
  sp.df <- m$get_species()
  sp_count <- nrow(sp.df)
  
  for (sp_num in 1:sp_count) {
    sp_current <- sp.df[sp_num,]
    sp_current_demes <- m$get_demes()[[sp_num]]
    env.ras <- rasterFromEnv(m, envOrig, sp_num)
    
    # calculate and print a species summary - probably drop this once development is done
    sp.summary <- data.frame(speciesID = sp_current$id,
                              range = nrow(sp_current_demes),
                              total_amount = sum(sp_current_demes$amount,
                              niche0.position.mean = sp_current$niche_position_mean_0),
                              niche0.position.sd = sp_current$niche_position_sd_0,
                              niche1.breadth.mean = sp_current$niche_breadth_mean_0,
                              niche1.breadth.sd = sp_current$niche_breadth_sd_0,
                              niche1.sp.min = sp_current$niche_min_0,
                              niche1.sp.max = sp_current$niche_max_0
    )
    sp.summary <- round(sp.summary, 4)
    if (do.text.output) {
      list.for.text <- c(list(current.time=currentTime,
                              elapsed.time.total = difftime(Sys.time(), starttime_global),
                              elapsed.time.step = difftime(Sys.time(), starttime_timestep)),
                         as.list(sp.summary[1,]))
      text.update(list(species_range_niche=list.for.text))
    }

    # update the dynamic plot
    if (do.display) {
      if (image_to_file) {
        display.to.file.start(image_dir, currentTime, image_filename = "animation")
      }
   
      model.params <- list(niche.breadth=round(sp.summary$niche1.breadth.mean), niche.evolution.rate=niche.evolution.rate, dispersal = dispersal)
browser()
      display.update(list(env = env.ras, demes_amount_position=sp_current_demes, current.time=currentTime, model.params = model.params))
 
      if (do.display.diff) {
        display.update(list(env = env.ras, demes_amount_position_diff = sp_current_demes, current.time = currentTime))
      }
      
      if (do.display.genome) {
        display.update(list(env = env.ras,
                            demes_genecolour = sp_current_demes,
                            current.time = currentTime,
                            genome.columns = genome.columns))
      }
 
      if (image_to_file) {
        display.to.file.stop()
      }
    }
    
  }  

}