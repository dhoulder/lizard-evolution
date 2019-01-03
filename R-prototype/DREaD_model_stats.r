
##### required libraries
require(data.table)
require(dplyr)
#require(gstat)
#require(raster)
library(yaml)

source("C:/Users/u3579238/Work/Software/dan-github/DREaD_ds/R-prototype/dataStructures.r")

###############################################################################

# range_though_time <- function(input.dir = "",
                     # input.prefix = "") {

input.dir <- "C:/Users/u3579238/Work/Simulation/alps250_test5_40/run_1"
input.prefix <- "out"
genome.dimensions <- 4

  starttime_global <- Sys.time()

  # first check for input files
  if (input.dir == "") {stop("A directory for model output files must be provided.")}

  input.filter <- paste(input.prefix, "[[:digit:]]*.csv", sep="")
  input.files  <- list.files(path=input.dir, pattern = input.filter, include.dirs = F)

  # get input file numbers to order correctly
  input.files.num <- sub(".csv", "", input.files) # remove the suffix
  if (input.prefix != "") {
    input.files.num <- sub(input.prefix, "", input.files.num)
  }

  input.file.index <- data.frame(index=1:length(input.files), number=as.numeric(input.files.num), filename=input.files)
  input.file.index <- input.file.index[order(input.file.index$number), ]
  rm(input.files, input.files.num)

  nfiles <- nrow(input.file.index)
  total.time <- nfiles

  input.file.path   <- paste(input.dir, "/", input.file.index$filename[1], sep="")
  first.input.file  <- read.csv(input.file.path)

  presence.cells <- 1:nrow(first.input.file)
  coords <- first.input.file[,2:3]
  rownum <- 1

  #demetable <- makeDemeTable(genome.dimensions = genome.dimensions, rowcount = 10000)
  genomeInitial <- as.list(rep(0, genome.dimensions))
  demetable.columncount <- ncol(demetable)

  # identify the genome position columns - to avoid repeating it within loops
  first.gene.col.idx    <- which(names(demetable)=="gene.pos1")
  genome.columns        <- first.gene.col.idx:(first.gene.col.idx + genome.dimensions - 1)

  input.genome.columns  <- grep("genetic", names(first.input.file))

  # for (i in 1:length(presence.cells)) {
  #   #cell <- presence.cells[i]
  # 
  #   new.row <- list(cellFromRowCol(env, first.input.file$row[i], first.input.file$column[i]), # cellID
  #                   1,                                    # speciesID
  #                   first.input.file$column[i],           # col
  #                   first.input.file$row[i],              # row
  #                   1,                                    # amount - need to sort out values!
  #                   first.input.file$niche_centre_0[i],   # niche1.position
  #                   first.input.file$niche_breadth_0[i],  # niche1.breadth
  #                   first.input.file$niche_centre_1[i],   # niche2.position
  #                   first.input.file$niche_breadth_1[i]   # niche2.breadth
  #   )
  # 
  #   new.row <- c(new.row, first.input.file[i, input.genome.columns]) # add the genome columns - 2 or more
  #   set(demetable, i=i, j=1:demetable.columncount, value=new.row)
  # }

  # loop through time steps in the simulation
  while(current.time <= total.time) {

    starttime_timestep <- Sys.time()

    # read in the input data for the current timestep
    input.file.path   <- paste(input.dir, "/", input.file.index$filename[current.time], sep="")
    current.input.file  <- read.csv(input.file.path)
    input.demecount <- nrow(current.input.file)

    # update demetable from the current.input.file
    demetable <- makeDemeTable(genome.dimensions = genome.dimensions, rowcount = input.demecount)
    demetable$cellID <- cellFromRowCol(env, current.input.file$row, current.input.file$column)
    demetable[, 2:9] <- current.input.file[, c(1,3,2,4,6,7,9,10)]
    demetable[, genome.columns] <- current.input.file[, input.genome.columns]

    # iterate through extant species
    for (i in species.tips) {

      #skips rows that have speciated (i.e., is an ancestral branch not an extant lineage)
      # this section shouldn't be needed if internal branches are correctly coded when speciation occurs
      extant.species <- which(!is.na(edgetable[,10]))
      if (length(extant.species) > 1) {
        if (edgetable[i, 2] %in% edgetable[extant.species, 1] ) {
          edgetable[i, 10] <- 2 #set this species as an internal branch
          next
        }
      }

      # # select species current iteration
      current.speciesID <- edgetable[i, 11]

      demetable.species <- demetable[demetable$speciesID==current.speciesID, ]

      env.table <- as.data.table(cbind(all.coords, env[]))
      names(env.table)[4] <- "env1"
      setkey(env.table, cellID)

      # # temporary test of within species genetic distances
      # genome.distances <- dist(demetable.species[, genome.columns, with=FALSE])
      # genome.distance.max     <- max(genome.distances)
      # genome.distance.median  <- median(genome.distances)

      # calculate and print a species summary - probably drop this once development is done
      sp.summary <- demetable.species[, .(range = .N,
                                          total_amount = sum(amount),
                                          niche1.position.mean = mean(niche1.position),
                                          niche1.position.sd = sd(niche1.position),
                                          niche1.breadth.mean = mean(niche1.breadth),
                                          niche1.breadth.sd = sd(niche1.breadth),
                                          niche1.sp.min = min(niche1.position - (niche1.breadth/2)),
                                          niche1.sp.max = max(niche1.position + (niche1.breadth/2))
                                          #gen.distance.max = genome.distance.max,
                                          #gen.distance.median = genome.distance.median
                                          )]
      sp.summary <- round(sp.summary, 4)
      if (do.text.output) {
        list.for.text <- c(list(current.time=current.time,
                                elapsed.time.total = difftime(Sys.time(), starttime_global),
                                elapsed.time.step = difftime(Sys.time(), starttime_timestep),
                                current.speciesID=as.integer(current.speciesID)),
                           as.list(sp.summary[1,]))
        text.update(list(species_range_niche=list.for.text))
      }

      # update the dynamic plot
      if (do.display) {
        if (image_to_file) {
          display.to.file.start(image_dir, current.time, image_filename = "animation")
        }

        niche.params <- list(niche.breadth=round(demetable.species[,max(niche1.breadth)],2), niche.evolution.rate=niche.evolution.rate, dispersal = dispersal)
        display.update(list(env = env, demes_amount_position=demetable.species, current.time=current.time, niche.params = niche.params))

        if (do.display.diff) {
          display.update(list(env = env, demes_amount_position_diff = demetable.species, current.time = current.time))
        }

        if (do.display.genome) {
          display.update(list(env = env,
                              demes_genecolour = demetable.species,
                              current.time = current.time,
                              genome.columns = genome.columns))
        }

        if (do.animate) {
          ani.record() # record the current frame
        }

        if (image_to_file) {
          display.to.file.stop()
        }
      }

      if (raster_to_file & current.time %in% generations_to_save) {
        current.species.ras <- demes_to_raster(demetable.species, current.speciesID, env)
        if (run_number_in_filename) {
          raster.filename <- paste(raster_dir, "run",  run_number, "species", current.speciesID, "_time", current.time, ".asc", sep="")
        } else {
          raster.filename <- paste(raster_dir, "species", current.speciesID, "_time", current.time, ".asc", sep="")
        }
        writeRaster(current.species.ras, raster.filename, overwrite=T)
      }

      gc(verbose=FALSE)

    } #end of looping through species
    
    current.time <- current.time + timestep.size

  }   # end of while loop

  # return()
# }
