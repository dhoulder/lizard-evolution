# this file contains functions for dispersal, competition, niche evolution and drift


######### rangeDispersal #########

# Function extends a species range based on a dispersal kernal and the inherent niche breadth of the species

############# Arguments ###############

# env = current environmental layer
# position = niche position of species 
# breadth = niche breadth of species
# species.ras is the current species range in raster format (presence/absence)

disperse_ds <- function (demetable.species, 
                         env, 
                         env.table, 
                         dispersal.range,
                         suitability.mode="block"){
  
  # group demes by similarity of the niche - so for very similar niches the environment is calculated just once
  
  # this is a clunky grouping method - there may be more efficient functions in data.table
  demetable.species$niche1.breadth.group <- as.integer(round(demetable.species$niche1.breadth/niche.blocksize))
  demetable.species$niche1.position.group <- as.integer(round(demetable.species$niche1.position/niche.blocksize))
  demetable.species$niche2.breadth.group <- as.integer(round(demetable.species$niche2.breadth/niche.blocksize))
  demetable.species$niche2.position.group <- as.integer(round(demetable.species$niche2.position/niche.blocksize))
  
  # create groupos of demes with similar environemnt to calculate dispersal
  niche.groups <- demetable.species[, list(niche1.breadth.group, niche1.position.group, niche2.breadth.group, niche2.position.group, .N), by=list(niche1.breadth.group, niche1.position.group, niche2.breadth.group, niche2.position.group)]
  niche.groups <- niche.groups[, list(niche1.breadth.group, niche1.position.group, niche2.breadth.group, niche2.position.group, N)]
  
  # loop through the niche groups
  for (k in 1:nrow(niche.groups)) {
    niche <- niche.groups[k, ]
    niche.values    <- niche * niche.blocksize  # turns the rounded niche group into niche values
    
    demetable.nichegroup <- demetable.species[niche1.breadth.group == niche$niche1.breadth.group & 
                                    niche1.position.group== niche$niche1.position.group &
                                    niche2.breadth.group == niche$niche2.breadth.group &
                                    niche2.position.group== niche$niche2.position.group, 1:12]
    
    # filter env cells to those within the spatial limits (for the niche group) 
    bounds        <- as.list(demetable.species[, list(xmin=(min(x)-dispersal.range), xmax=(max(x)+dispersal.range), ymin=(min(y)-dispersal.range), ymax=(max(y)+dispersal.range))])
    bounds[bounds <0] <- 0
    #niche.extent  <- extent(bounds[[1]], bounds[[2]], bounds[[3]], bounds[[4]])

    env.table.dispersal <- env.table[col >= bounds$xmin & 
                                      col <= bounds$xmax &
                                      row >= bounds$ymin &
                                      row <= bounds$ymax, ]
    #env.dispersal <- crop(env, niche.extent)
browser()     
    # filter env cells to those within the niche limits
    niche1.min <- niche.values$niche1.position.group - niche.values$niche1.breadth /2
    niche1.max <- niche.values$niche1.position.group + niche.values$niche1.breadth /2    
    #niche2.min <- niche.values$niche2.position.group - niche.values$niche2.breadth /2
    #niche2.max <- niche.values$niche2.position.group + niche.values$niche2.breadth /2    
    
    env.table.dispersal <- env.table.dispersal[env1 >= niche1.min & env1 <= niche1.max, ]
    # dispersal.cells <- which(!is.na(env.dispersal[]))
    
    #apply niche suitability function to env.dispersal to give suitability for 0 to 1
    env.table.dispersal <- niche_suitability(env=env.table.dispersal, suitability.mode = suitability.mode, niche1.breadth = niche.values$niche1.breadth, niche1.position = niche.values$niche1.position)
    
    demetable.nichegroup.new <- demetable.nichegroup[0, ]
    row.pointer <- 0

#TEMPTEMPTEMP
    points(env.table.dispersal$col, env.table.dispersal$row, col="blue", pch=20, cex=0.8)
    points(demetable.species$x, demetable.species$y, col="red", pch=20, cex=0.6)

    
    # loop through each deme for the niche group
    for (d in 1:nrow(demetable.nichegroup)) {

      deme <-  demetable.nichegroup[d,]
      #TEMPTEMPTEMP
      points(deme$x, deme$y, col="black", pch=20, cex=1.5)
      
      # find the cells in dispersal distance
      deme.dest <- env.table.dispersal[col >= deme$x-dispersal.range & col <= deme$x+dispersal.range, ]
      deme.dest <- deme.dest[row >= deme$y-dispersal.range & row <= deme$y+dispersal.range, ]
      new.count <- nrow(deme.dest)
      
      new.amount  <- deme$amount * deme.dest$suitability  # should return a vector
      new.rows    <- (row.pointer+1):(row.pointer+new.count)
      
      if (length(new.amount) > 0) {
        # create dispersed demes
        demetable.nichegroup.new <- rbind(demetable.nichegroup.new,
                                          list(cellID=deme.dest$cellNum,
                                               x=deme.dest$col,
                                               y=deme.dest$row,
                                               amount=deme$amount * deme.dest$suitability), fill=T)
        demetable.nichegroup.new$speciesID[new.rows] <- deme$speciesID
        demetable.nichegroup.new$niche1.position[new.rows] <- deme$niche1.position
        demetable.nichegroup.new$niche1.breadth[new.rows]  <- deme$niche1.breadth
        demetable.nichegroup.new$niche2.position[new.rows] <- deme$niche2.position
        demetable.nichegroup.new$niche2.breadth[new.rows]  <- deme$niche2.breadth
        demetable.nichegroup.new$gene.pos1[new.rows] <- deme$gene.pos1
        demetable.nichegroup.new$gene.pos2[new.rows] <- deme$gene.pos2
        demetable.nichegroup.new$gene.pos3[new.rows] <- deme$gene.pos3
browser()         
        # apply distance function
        weights <- distance.weights(source = deme[1, c("x", "y")], 
                                    destinations = demetable.nichegroup.new[new.rows, c("x", "y")],
                                    dispersal.range = dispersal.range,
                                    distance.type = "euclidean", 
                                    distance.decay = "linear")
        
        demetable.nichegroup.new$amount[new.rows] <- demetable.nichegroup.new$amount[new.rows] * weights 

      }

      row.pointer <- row.pointer + new.count
      
      #TEMPTEMPTEMP
      points(deme$x, deme$y, col="white", pch=20, cex=1.5)
      
      print(deme)
    }
  }
  
  print(demetable.nichegroup)
  
  # loop through each of the grouped niches for the species
  
    # apply its niche function to env

    # loop through each deme for the niche group
  
      # apply dispersal based on the niche, and add new demes to a temporary demetable which is just for this species, with co-occurrence
  
      # record an initial amount/influence based on proximity and niche fit
  
  # now resolve all of the co-occurences based on
    # prior occupation
    # amount
    # genetic distance
  
    # to call an occurence of
      # a) combining with gene flow -> new single occurrence
      # b) competitive exclusion
      # c) co-occurrence ??
  
  # return modified edgetable
  
}

co_occur_ds <- function (demetable_sp){
  # determine what happens to populations which are in the same cell
  
  
}

niche_suitability <- function(env,
                              suitability.mode,
                              niche1.position, 
                              niche1.breadth, 
                              niche2.position=0, 
                              niche2.breadth=0) {
  # convert environment to suitability

  env.class <- class(env)  # allow env to be either a raster, or a data.table, and handle accordingly

  niche1.min <- niche1.position - (niche1.breadth /2)
  niche1.max <- niche1.position + (niche1.breadth /2)
  niche2.min <- niche2.position - (niche2.breadth /2)
  niche2.max <- niche2.position + (niche2.breadth /2)

  if (suitability.mode=="block") {
    # block suitability is 1 within the suitable range, zero elsewhere
    if (any(env.class=="data.table")) {
      env$suitability <- env$env1
      env[env1 < niche1.min | env1 > niche1.max, "suitability"] <- 0
      env[env1 > 0, "suitability"] <- 1
    } else if (env.class=="raster") {
      suitability[suitability < niche1.min | suitability > niche1.max] <- 0
      suitability[suitability > 0] <- 1
    }
    
  } else if (suitability.mode=="sine") {
    # sine suitability is 1 for the central niche position, declining to 0 at
    # the niche limits according to a sine curve
    
    # transform suitability to range from 0 (for minimum to pi for maximum)
    if (any(env.class=="data.table")) {
      env$suitability <- env$env1
      env$suitability <- env[, suitability - niche1.min]
      if (min(env$suitability) < 0) {
        env[suitability < 0, suitability] <- 0
      } 
      env$suitability <- env$suitability * (pi / niche1.breadth)
      if (max(env$suitability) > pi) {
        env[suitability > pi, suitability] <- 0
      }
      env$suitability <- sin(env$suitability)
    } else if (env.class=="raster") {
      suitability <- suitability - niche1.min
      suitability[suitability < 0] <- 0
      suitability <- suitability * (pi / niche1.breadth)
      suitability[suitability > pi] <- 0
      suitability <- sin(suitability)
    }

  }
  
  if (any(env.class=="data.table")) {
    return(env)
    } else {
    return(suitability)
  }
}

distances <- function(source, destinations, distance.type="euclidean") {
  # calculates a vector of distances based on:
  #   a single source location; and 
  #   one or more destination locations, provided as a two column data type 
  
  # for a single distance type this function is not needed, but useful to allow easy
  # switching between distance types - for example euclidean or cost based

  if (distance.type=="euclidean" | substr(distance.type, 1, 3) == "euc") {
    distances <- raster::pointDistance(source, destinations, lonlat = FALSE)
  }
  
  return(distances)
}

reduce.to.max     <- function(x, dispersal.range) {return(min(x, dispersal.range))}
set.to.zero  <- function(x, threshold) {
  if (x < threshold) {x <- 0}
  return(x)
}

distance.weights  <- function(source,
                              destinations, 
                              dispersal.range, 
                              distance.type="euclidean", 
                              distance.decay="linear") {
  # calculates a weighting (0 to 1) which adjusts, by distance, the amount of a
  # deme arriving in a cell
  
  # calculation is based on:
  #   a single source location; and 
  #   one or more destination locations, provided as a two column data type
  #   a parameter for the type of distance measure - at present only euclidean diatnce is available
  #   a distance decay function which defines how amount relates to distance
  
  distances <- distances(source, destinations, distance.type)
  
  if (distance.decay=="linear") {
    min.weight <- 0.2  # this means that at the maximum distance, the weight is 0.2, rather than declining to 0
    
    #distances         <- sapply(distances, FUN=reduce.to.max, dispersal.range)
    distance.weights  <- (dispersal.range - (distances * (1-min.weight))) / dispersal.range
    distance.weights  <- sapply(distance.weights, FUN=set.to.zero, threshold=min.weight)
  }

  return(distance.weights)
}