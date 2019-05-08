rm(list=ls())

library(dreadds)
library(raster)
library(data.table)
library(ape)

source("~/code/DREaD_ds/R-prototype/dynamicDisplay.r")

# turn the input arguments into local variables with relevant names
input.args	<- commandArgs(trailingOnly = TRUE)

if (length(input.args) > 0) {
  output.dir								<- input.args[1]
  config.file 							<- input.args[2]
  dispersal									<- as.numeric(input.args[3])
  niche.evolution.rate			<- as.numeric(input.args[4])
  gene.flow.max.distance   	<- as.numeric(input.args[5])
  cat("\nExternal arguments received in R\n\toutput.dir\t", output.dir, 
																					"\n\tconfig.file\t", config.file, 
																					"\n\tDispersal\t", dispersal, 
																					"\n\tniche evolution rate\t", niche.evolution.rate, 
																					"\n\tgene flow max distance\t", gene.flow.max.distance, "\n")
} else {
  output.dir							<- "~/simulation/AMT_paleoclim/movie_5my_d2_e2_s12/output4/"
  config.file 						<- "/short/ka2/simulation/AMT_paleoclim/scripts/AMT_paleoclim_narrowniche.conf"
	dispersal               <- 2     # this should be extracted from the model arguments, but not yet possible
	niche.evolution.rate    <- 0.02  # this should be extracted from the model arguments, but not yet possible
	gene.flow.max.distance  <- 12    # this should be extracted from the model arguments, but not yet possible
  cat("\nNo external arguments received in R. Hard coded arguments used instead.\n\toutput.dir\t", output.dir, "\n\tconfig.file\t", config.file, "\n\tDispersal\t", dispersal, "\n\tniche evolution rate\t", niche.evolution.rate, "\n\tgene flow max distance\t", gene.flow.max.distance, "\n")
}

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

envOrig						<- raster("~/simulation/input_data/AMT2.5min-grids/T5000/Dan AMT2.5min_bio PaleoClimDat - MinTemp.txt.tif")
#envOrig						<- raster("/home/805/dfr805/code/DREaD_extras/habitat_sizes.asc")
environment.rows	<- nrow(envOrig)

coastline.shape		<- shapefile("~/simulation/input_data/shapefiles/aus5m_coast_no_norfolk.shp")

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
image_frequency   <- 4

runInterval       <- 4  # this is the number of steps to run before returning to R
nSteps						<- 2000

# define and create the directory structure
base.dir              <- "~/simulation/AMT_paleoclim/movie_5my_d2_e2_2000/"
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

# set up a data frame of extant species through time, initially with 0 rows
diversity_time.df <- data.frame(time=0, currentStep=0, current_species=0, extinct_species=0)
diversity_time.df <- diversity_time.df[-1,]

starttime_global <- Sys.time()

m <- createDreadDS(
        config.file = config.file,
        output.dir = output.dir,
        iterations = nSteps,
				niche.evolution.rate = niche.evolution.rate,
				gene.flow.max.distance = gene.flow.max.distance)
				
# print the model version
try(cat("Running dreadds version:", m$version(), "\n"))
# this line is wrapped in a 'try()' because it would fail on earlier versions of dreadds which
# lack the m$version() method

# step through the model 
currentStep     <- 0
currentInterval <- 1  # start plotting with a single step

while (currentStep <= nSteps) {
  
  starttime_timestep <- Sys.time()
  
  # run the next model step(s)
  m$runSteps(currentInterval)
  currentStep <- currentStep + currentInterval  
  
  # calculate the time before present for map and text outputs
  paleoTime <- paleoTime.start - (paleoTime.step * currentStep)

  ##### species summary and plots ######
  
  # get the species data frame and number of species
  sp.df <- m$getSpecies()
  sp_all_count <- nrow(sp.df)  # all species ever (ie count of branches)

  # write the species data frame to file (keeping only the latest version)
  write.csv(sp.df, paste(output.dir, "species.csv", sep=""))  

  # count extinct species, then subset sp.df to only extant species
  extinct_sp_count <- length(which(sp.df$extinction > 0))
	sp.df <- sp.df[which(sp.df$extinction == -1 & sp.df$split == -1), ]
  sp_count <- nrow(sp.df)

  # update diversity_time.df
  diversity_time.df <- rbind(diversity_time.df, c(paleoTime, currentStep, sp_count, extinct_sp_count))
  names(diversity_time.df) <- c("time", "currentStep", "current_species", "extinct_species")
  
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
 
  ###########################################################
  # experimenting with the phylogeny
  multiple_species_exist <- FALSE
  
  if (sp_count > 1) {
    multiple_species_exist    <- TRUE
    
    tree.text        <- m$getPhylogeny()
    
    # fix tree text to work with read.tree, by stripping the outer set of brackets
    if ( (substr(tree.text, 1, 2) == "((") & (substring(tree.text, nchar(tree.text)-1, nchar(tree.text)) == ");") ) {
      tree.text <- substring(tree.text, 2)
      substr(tree.text, nchar(tree.text)-1, nchar(tree.text)) <- "; "
      tree.text <- substring(tree.text, 1, nchar(tree.text)-1)
    }  
    
    tree      	     <- read.tree(text = tree.text)
  }
  ###########################################################

  # update the dynamic plot FOR ALL SPECIES at specified frequency
  if (do.display) {
    if (currentStep==1 | currentStep %% image_frequency == 0) {  # only do the image at the specified frequency
      
      # get the data needed for plots

      env.ras <- rasterFromEnv(m, envOrig, 1)
      
      if (image_to_file) {
        display.to.file.start(image_dir=image.dir, time=image_counter, image_filename = paste("animation_multisp_", sep=""), plot_rows=2, plot_cols=3)
				# use time=image_counter for consecutively numbered images (eg for ffmpeg) or time=currentStep to number by the model time
      } else {
        display.to.screen.start(window_name = paste("plot_", image_counter, sep=""), plot_rows=2, plot_cols=2)
      }

      model.params <- list(niche.breadth=round(sp.df$niche_breadth_mean_0[1]), niche.evolution.rate=niche.evolution.rate, dispersal = dispersal, speciation.distance = gene.flow.max.distance)

			if (multiple_species_exist) {
			# plot the main header and phylogeny
				display.update.multispecies(list(current.time = paleoTime, model.params = model.params, tree = tree, diversity_time = diversity_time.df), includes.map = F, plot_tree = TRUE, plot_mainheader = TRUE)	  
			} else {
				# plot the main header and demes by amount and niche optimum (for environment 1)
				display.update.multispecies(list(current.time = paleoTime, model.params = model.params, shape = coastline.shape, diversity_time = diversity_time.df), env = env.ras, allDemes = all.demes, plot_demes_amount_position = TRUE, plot_mainheader = TRUE)
			}
				
      if (multiple_species_exist) {
        # plot species accumulation over time
        display.update.multispecies(list(current.time = paleoTime, model.params = model.params, diversity_time = diversity_time.df), includes.map = F, plot_species_over_time = TRUE)
      }
      
      #display.update.multispecies(list(current.time = paleoTime, model.params = model.params), env = env.ras, allDemes = all.demes, plot_species_ranges = TRUE)	  
      
      display.update.multispecies(list(current.time = paleoTime, model.params = model.params, speciesCount = sp_count, shape = coastline.shape), 
                                  env = env.ras, 
                                  allDemes = all.demes, 
                                  plot_richness = TRUE)
      
      if (do.display.genome) {
    
        if (!exists("genome.extremes")) {
          genome.extremes <- matrix(data=NA, nrow=2, ncol=length(genome.columns))
        }
        
        genome.extremes <- get.genome.extremes(all.demes, genome.columns, genome.extremes)
        display.update.multispecies(list(current.time=paleoTime, model.params=model.params, genome.columns=genome.columns, gene.flow.max.distance=gene.flow.max.distance, genome.extremes=genome.extremes, shape=coastline.shape),
                                    env = env.ras, 
                                    allDemes = all.demes,
                                    plot_genome_scatter = TRUE,
                                    plot_genome_map = TRUE)
      }

      if (image_to_file) {
        display.to.file.stop()
      } else {
        display.to.screen.stop()
      }
      image_counter <- image_counter + 1
    }
  }
  
  # allow the plot to begin at step one, then adjust to the requested stepping schedule
  if (runInterval==1 | currentStep >= runInterval) {
    currentInterval <- runInterval
  } else {
    currentInterval <- runInterval - currentStep # this allows the first two steps to be 1, runInterval - 1
  }
  
}
