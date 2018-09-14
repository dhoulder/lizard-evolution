#  SOME QUICK NOTES AND EDITS TO CONSIDER MODIFYING ROLLING STONE
#  TO SIMULATE EVOLUTION WITHIN SPECIES USING DEMES AS THE BASIC UNITS

# RollingStone has three basic data structures
# an environment raster
# a species raster stack
# an edge table which stores species / branches with their niche and their phylogenetic parent

# The intraspecific model would instead have
# an environment raster (initially same, but allow > 1 environmental dimension)
# a demes table which stores for each deme:
# location  as a cell index number (or numbers if it can occupy multiple cells)
# niche value (eg mean and breadth) for each environmental dimension
# genetic position in multiple dimensions - probably 3
# abundance and/or suitability in the cell
# species membership - to streamline gene flow calculations and generate species phylogeny

# a species table, updated at each timestep from the demes table, storing for each species (and branch above species level):
# its parent branch/edge
# status - a) current species, b) speciated (internal branch) or c) extinct
# timestep when it originated

# if helpful, the species table could also store or point to the cells / demes which comprise it, and record the overall niche limits of the species.  These would both be a form of indexing, with the authoritative values at the level of the deme

# a raster stack of distributions of current species would be an option - but again, it is a result not a primary authority for ranges.



######### Dynamic Range Evolution and Diversification model (DREaD) #########

# main simulator function of geographic diversification model

############# Arguments ###############

# dispersal = dispersal kernal parameter
# enviro.mode = environmental change mode - character "sine" or "linear"
# amp = amplitude of environmental change sine wave
# freq = frequency of environmnental change sine wave
# slope = slope of environmental change linear model
# breadth.ev.rate = rate of niche breadth evolution - not currently used
# niche.evolution.rate = rate of niche evolution - impact depends on the mode of niche evolution
# enviro.hetero = biniary variable. whether environmental change model varies spatially or operates in synchrony across the domain
# plot = whether to plot species ranges and phylogeny at end of simulation (plot results)
# animate = whethe to save a snapshot of species range at each time step (png) to turn into a GIF using otehr software (visualisation tool)
# timestep.size = size of each time step
# generateSummaryStatistics = logical. whether to generate the 32 summary statistics

# Arguments specific to the dynamic speciation model

# total.time = the amount of time to elapse while running the model.  This replaces the number of tips to determine how long the model runs
# niche.blocksize = the similarity of niches, that are treated as the same for calculating dispersal and 'abundance'
  # a larger number should lead to faster running, but could miss effects of small niche differences
  # each niche axis is scaled from 0 to 25

DREaD_ds <- function(total.time,
                      dispersal,
                      niche.evolution.rate,
                      breadth.ev.rate,
                      enviro.hetero,
                      enviro.mode,
                      amp = NA,
                      freq = NA,
                      slope = NA,
                      plot = FALSE,
                      timestep.size = 1,
                      generateSummaryStatistics = TRUE,
                      genome.dimensions = 3,
                      niche.blocksize = 0.1,
                      suitability.mode="block",
                      speciation.gene.distance,
                      environment.source,
                      initial.species.defined=NA,
                      environment.dimension = 100) {

  ##### required libraries
  #require(ape)
  #require(apTreeshape)
  require(data.table)
  require(dplyr)
  #require(ENMTools)
  #require(fossil)
  #require(geiger)
  require(gstat)
  #require(moments)
  #require(phyloclim)
  #require(phytools)
  require(raster)

  starttime_global <- Sys.time()

  # matrix describing the degree of environmental change across the landscape (for enviro.hetero=T)
  env.change.matrix <- matrix(rep(seq(from=0.01, to =1, by=1/environment.dimension), environment.dimension), ncol=environment.dimension, nrow=environment.dimension, byrow=F)

  # vector of parameters
  params <- data.frame(total.time=total.time, dispersal=dispersal, amp=amp, freq=freq,
              niche.evolution.rate=niche.evolution.rate, breath.ev.rate=breadth.ev.rate,
              slope=slope, enviro.mode=enviro.mode, enviro.hetero=enviro.hetero)

  # define the minimum occurrence amount (between 0 and 1).  Values below this will be treated as 0
  minimum.amount <- 0.001

  ####################### CHECK INPUT ARGUMENTS - just checking some at this point ########

  # slope parameter required only for linear climate models, amp and freq only for sine climate models
  if (enviro.mode == "linear") {
    if (is.na(slope)) {stop("When enviro.mode is 'linear', A value must be provided for slope.")}
  } else if (enviro.mode == "sine") {
    if (is.na(amp) | is.na(freq)) {stop("When enviro.mode is 'sine', values must be provided for amp and freq.")}
  }

  ######################### 1. Generate or load background environment ########################

  if (environment.source == "internal") {
    env <- generateEnv(grid.size = environment.dimension, original = T)
  } else {
    env <- raster(environment.source)
  }
  starting.env <- env

  # use the background environment to define all cells in the model
  all.coords <- rowColFromCell(env, 1:ncell(env))
  all.coords <- as.data.table(cbind(1:nrow(all.coords), all.coords))
  names(all.coords)[1] <- "cellID"

  if (do.display & !do.display.diff & !image_to_file) {
    display.update(list(env = env))
  }

  ###########################  2. Seed initial species and data structures #############################
  #initial.species <- seedSpecies(env, dispersal = dispersal, print.species.params = TRUE)

  initial.species <- seedSpecies_defined(env, dispersal = dispersal, initial.species.defined)
  initial.species.ras <- initial.species[[1]]

  # generate edgetable  -  edgetable is a matrix that stores information on each species'
  # phylogeny, niche position, niche breadth, speciation modes, range size each row is a
  # species

  edgetable <- makeEdgeTable(1000, dynamicSpeciation = TRUE)

  edgetable[1,5]  <- sum(initial.species.ras@data@values)
  edgetable[1,7]  <- initial.species[[2]]
  edgetable[1,8]  <- initial.species[[3]]
  edgetable[1,11] <- 1

  presence.cells <- which(initial.species.ras[] == 1)
  coords <- rowColFromCell(initial.species.ras, presence.cells)
  rownum <- 1

  demetable <- makeDemeTable(genome.dimensions = genome.dimensions, rowcount = 10000)
  genomeInitial <- as.list(rep(0, genome.dimensions))
  demetable.columncount <- ncol(demetable)

  # identify the genome position columns - to avoid repeating it within loops
  first.gene.col.idx    <- which(names(demetable)=="gene.pos1")
  genome.columns        <- first.gene.col.idx:(first.gene.col.idx + genome.dimensions - 1)

  for (i in 1:length(presence.cells)) {
    cell <- presence.cells[i]
    new.row <-  list(cell,                  # cellID
                  1,                     # speciesID
                  coords[i, 2],          # col
                  coords[i, 1],          # row
                  1,                     # amount - need to sort out values!
                  initial.species[[2]],  # niche1.position
                  initial.species[[3]],  # niche1.breadth
                  0,                     # niche2.position - need to sort out values if using
                  0                     # niche2.breadth  - need to sort out values if using
    )
    new.row <- c(new.row, genomeInitial) # add the genome columns - 2 or more

    set(demetable, i=i, j=1:demetable.columncount, value=new.row)
  }

  demetable.used.rows <- i

  # # species rasters is a list that hangs onto each species geographic range in the form of a raster
  # species.rasters <- vector('list', 10000)
  # species.rasters[[1]] <- initial.species[[1]]

  current.time <- 0
  tips <- 1

  extinct <- vector("logical", 1000)
  extinct.number <- 0

  # while loop propels the simulation. iterations repeat until the condition (number of species generated) is met

  while(current.time < total.time) {

  starttime_timestep <- Sys.time()

  ############################# 3. Environmental change ###########################

    # time changes
    current.time <- current.time + timestep.size

    # environment changes
    env <- enviroChange(start.env=starting.env, env=env, current.time=current.time, amp=amp, freq=freq, slope=slope,
                        model= enviro.mode, hetero=enviro.hetero, env.change.matrix=env.change.matrix,
                        env.dimension=environment.dimension)
    # species.tips is the row index of non-extinct lineages (rows in the edgetable)
    if(any(extinct==TRUE)){
      species.tips <- seq_along(which(!is.na(edgetable[,10])))[-which(extinct == TRUE)]
    } else {
      species.tips <- seq_along(which(edgetable[,10]==1))
    }

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

  ############################# 4. Dispersal ###########################
      demetable.species <- demetable[demetable$speciesID==current.speciesID, ]

      env.table <- as.data.table(cbind(all.coords, env[]))
      names(env.table)[4] <- "env1"
      setkey(env.table, cellID)

      # run the deme dispersal function
      demetable.species.overlap <- disperse_ds(demetable.species, env=env, env.table, dispersal.range=dispersal, suitability.mode=suitability.mode)

      # check for extinction here
      if (nrow(demetable.species.overlap) == 0) {
        edgetable <- extinction(edgetable, current.speciesID, current.time)
      }

      demetable.species <- combine.demes(demetable.species.overlap, genome.columns, speciation.gene.distance, minimum.amount, env.table, verbose=FALSE)

      # check for extinction here
      if (nrow(demetable.species) == 0) {
        edgetable <- extinction(edgetable, current.speciesID, current.time)
      }

  ############################  5. Evolution ######################

      # niche evolution for each deme
      demetable.species <- niche.evolution(demetable.species, env.table, niche.evolution.rate)

      # genetic drift for each deme
      demetable.species <- genetic_drift(demetable.species, timestep.size, genome.columns)

      # temporary test of within species genetic distances
      genome.distances <- dist(demetable.species[, genome.columns, with=FALSE])
      genome.distance.max     <- max(genome.distances)
      genome.distance.median  <- median(genome.distances)

      # calculate and print a species summary - probably drop this once development is done
      sp.summary <- demetable.species[, .(range = .N,
                                      total_amount = sum(amount),
                                      niche1.position.mean = mean(niche1.position),
                                      niche1.position.sd = sd(niche1.position),
                                      niche1.breadth.mean = mean(niche1.breadth),
                                      niche1.breadth.sd = sd(niche1.breadth),
                                      niche1.sp.min = min(niche1.position - (niche1.breadth/2)),
                                      niche1.sp.max = max(niche1.position + (niche1.breadth/2)),
                                      gen.distance.max = genome.distance.max,
                                      gen.distance.median = genome.distance.median)]
      sp.summary <- round(sp.summary, 4)
      if (do.text.output) {
        list.for.text <- c(list(current.time=current.time,
                                elapsed.time.total = difftime(Sys.time(), starttime_global),
                                elapsed.time.step = difftime(Sys.time(), starttime_timestep),
                                current.speciesID=as.integer(current.speciesID)),
                          as.list(sp.summary[1,]))
        text.update(list(species_range_niche=list.for.text))
      }

      # update this species in demetable
      demetable <- demetable[demetable$speciesID != current.speciesID, ]
      demetable <- rbind(demetable, demetable.species)

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
          raster.filename <- paste(raster_dir, "species", current.speciesID, "_run",  run.number, "_time", current.time, ".asc", sep="")
        } else {
          raster.filename <- paste(raster_dir, "species", current.speciesID, "_time", current.time, ".asc", sep="")
        }
        writeRaster(current.species.ras, raster.filename, overwrite=T)
      }

      gc(verbose=FALSE)

    } #end of looping through species

  #   # if all species are extinct reset to the starting values of the simulation
  #     if (all(edgetable[, 10][which(!edgetable[, 2] %in% edgetable[, 1])] =="extinct")) {
  #       initial.species <- seedSpecies(env, dispersal=dispersal)
  #       edgetable <- matrix(ncol=10, nrow=10000)
  #       edgetable[1,] <- c(0, 1, 0, NA, NA, 1, NA, NA, NA, "X")
  #       edgetable[1,5] <- sum(initial.species[[1]]@data@values, na.rm=T)
  #       edgetable[1,7] <- initial.species[[2]]
  #       edgetable[1,8] <- breadth
  #       species.rasters <- vector('list', 10000)
  #       species.rasters[[1]] <- initial.species[[1]]
  #       current.time <- 1
  #       tips <- 1
  #       extinct <- vector("logical", 10000)
  #       # Record mass extinction event. When 5 mass extinction events occur (extinction =5) in a row the simulation will sample new parameters and start again
  #       extinction <- extinction + 1
  #       extinct.number <- 0
  #       if(extinction == 5) {return(list(params, "too many extinctions with these parameters"))}
  #       next
  #     }
  #     extinct.species <- which(extinct==TRUE)
  #     extinct.number <- length(extinct.species)
  #     terminal.branches <- which(!edgetable[, 2] %in% edgetable[, 1])
  #     tips <- length(terminal.branches)
     }# end of while loop
  #
  # # species.rasters.2 are extant species rasters (extinct rasters are empty)
  # species.rasters.2 <- species.rasters[terminal.branches[!terminal.branches %in% extinct.species]]
  # # Build phylogeny
  # phy.table <- buildPhyInSim(edgetable, species.rasters.2@layers, extinct, tips)
  #
  # # plot species ranges
  # if (plot == TRUE ) {
  #   cols <- sample(rainbow(100,end=1, start=0, alpha=0.6 ), 100, replace = TRUE)
  #   par(mfrow = c(2, 1), mar=c(2,2,2,2))
  #   plot(env)
  #   for (i in 1:length(species.rasters.2)) {
  #     par(new=T)
  #     sp<-rasterToPolygons(species.rasters.2[[i]], dissolve=T)
  #     plot(sp, add=T, col=cols[i])
  #   }
  #   plot(phy.table[[2]]); axisPhylo()
  # }
  #
  # species.rasters.2 <- stack(species.rasters.2)
  # # record proportion of domain occupied by the whole clade
  # clade.area <- sum(stackApply(species.rasters.2, indices = rep(1, length(species.rasters.2@layers)), fun=mean)@data@values, na.rm=T)
  # #generate summary statistics (must have ENMTools loaded)
  # if(generateSummaryStatistics == TRUE) {
  #   summaryStats <- try(generateSummaryStatistics(phy.table[[1]], species.rasters.2, phy.table[[2]], drop.fossil = T))
  #   ARC.clade <- createENMToolsclade(species.rasters.2, phy.table[[1]], phy.table[[2]], env, drop.fossil = T)
  #   ARC <- try(enmtools.aoc(clade = ARC.clade,  nreps = 0, overlap.source = "range"))
  #   if(class(summaryStats) == "try-error") {
  #     summary.vector=NA
  #   } else {
  #     if(class(ARC)=="try-error"){
  #       summary.vector = c(summaryStats$analysis, NA, NA)
  #       ARC <- NA
  #     } else {
  #       summary.vector = c(summaryStats$analysis, ARC$coefficients[[1,1]], ARC$coefficients[[1,2]])
  #     }
  #   }
  #   results <- list("df" = phy.table[[1]], "phy" = phy.table[[2]], "rasters" = species.rasters.2, "env" = env, "params"=params, "ARC"=ARC, "summaryStats"=summaryStats, "summaryVector"= summary.vector , "extinct"= length(which(extinct==TRUE)), "clade.area"= clade.area)
  # } else {
  # results <- list("df" = phy.table[[1]], "phy" = phy.table[[2]], "rasters" = species.rasters.2, "env" = env, "params"=params, "extinct"= length(which(extinct==TRUE)), "clade.area"=clade.area)
  # }


# returns a list with the following elements 1) data frame with species information, 2) phylogeney, 3) environment (as was at end of simulation), 4) parameters used in that simulation
# 5) the Age-Range_correlation object crated by ENMTools, 6) a list of summary statistics (see generateSummaryStatistics function for details), 7) number of extinct lineages, 8) total area occupied by clade
  return()
}
