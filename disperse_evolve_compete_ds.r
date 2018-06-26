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
browser()    
    # filter env cells to those within the spatial limits (for the niche group) 
    bounds        <- as.list(demetable.species[, list(xmin=(min(x)-dispersal.range), xmax=(max(x)+dispersal.range), ymin=(min(y)-dispersal.range), ymax=(max(y)+dispersal.range))])
    bounds[bounds <0] <- 0
    #niche.extent  <- extent(bounds[[1]], bounds[[2]], bounds[[3]], bounds[[4]])

    env.table.dispersal <- env.table[col >= bounds$xmin & 
                                      col <= bounds$xmax &
                                      row >= bounds$ymin &
                                      row <= bounds$ymax, ]
    #env.dispersal <- crop(env, niche.extent)
   
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

      # find the cells in dispersal distance
      deme.dest <- env.table.dispersal[(col >= deme$x-dispersal.range 
                                       & col <= deme$x+dispersal.range
                                       & row >= deme$x-dispersal.range
                                       & row <= deme$y+dispersal.range), ]
      #deme.dest <- deme.dest[row >= deme$y-dispersal.range & row <= deme$y+dispersal.range, ]
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
        demetable.nichegroup.new$originCell[new.rows] <- deme$cellID
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

    }

    demetable.nichegroup.new <- demetable.nichegroup.new[amount > 0, ]
    
    # combine the niche group demetable rows
    if (exists("demetable.species.new")) {
      demetable.species.new <- rbind(demetable.species.new, demetable.nichegroup.new)
    } else {
      demetable.species.new <- demetable.nichegroup.new
    }
  }

  #TEMPTEMPTEMP
  newdemes <- unique(demetable.species.new[, c("x","y")])
  points(newdemes$x, newdemes$y, col="white", pch=20, cex=1.5)
  
  return(demetable.species.new)

}

combine.demes <- function (demetable.species.overlap, genomeDimensions, speciation.gene.distance){
  # determine what happens to populations which are in the same cell
  
  demetable.species <- demetable.species.overlap[0, -"originCell"]

  deme.cells <- unique(demetable.species.overlap$cellID)
  for (cell in deme.cells) {
    demetable.cell <- demetable.species.overlap[cellID==cell, ]
    deme.count     <- demetable.species.overlap[cellID==cell, .N]
   
    if (deme.count==1) {
      demetable.species <- rbindlist(list(demetable.species, demetable.cell[, -"originCell"]))
    } else {
      
      cat("For species", demetable.species.overlap[1, speciesID], "in cell", cell, "there are", deme.count, "origin demes to combine\n")
      
      # determine the primary deme, which the others may join with
      # if the current cell is one of the sources, it's the primary cell
      deme.primary <- which(demetable.cell$originCell==cell)
      
      if (length(deme.primary) == 0) {
        # otherwise pick a cell at random, wit probability weighted by amount
        deme.primary <- sample(x=deme.count, size=1, prob=demetable.cell$amount, replace = T)
      }
       
      # now loop through each of the source demes and combine them with the primary deme, as determined above
      demetable.cell$gene.flow <- is.geneflow(demetable.cell, deme.primary, dimensions=3, speciation.gene.distance)
      demetable.cell <- demetable.cell[gene.flow==TRUE,]
      
      if (nrow(demetable.cell) > 0) { # make sure that there are still demes with gene flow - better check why not - should be gene flow to self
        demetable.species.new <- demetable.cell[, list(amount=mean(amount),
                                  niche1.position=weighted.mean(niche1.position, amount), 
                                  niche1.breadth  =weighted.mean(niche1.breadth, amount),
                                  niche2.position =weighted.mean(niche2.position, amount), 
                                  niche2.breadth  =weighted.mean(niche2.breadth, amount),
                                  gene.pos1       =weighted.mean(gene.pos1, amount),
                                  gene.pos2       =weighted.mean(gene.pos1, amount)), 
                                  by=list(cellID, speciesID, x, y) ]
        if (genomeDimensions > 2) {
          for (g in 3:genomeDimensions) {
            column.name <- paste("gene.pos", g, sep="")
     
            #cat("\nDebugging:", weighted.mean(demetable.cell[[column.name]], demetable.cell[["amount"]])  ,"\n")
            demetable.species.new[[column.name]] <- weighted.mean(demetable.cell[[column.name]], demetable.cell[["amount"]])
          }
        }
        
        demetable.species <- rbindlist(list(demetable.species, demetable.species.new))
      }
    }

  }
browser()   
  return(demetable.species)
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

gene.distances <- function(source, destinations, dimensions=3) {
  # calculates a vector of distances based on:
  #   a single source location; and 
  #   one or more destination locations, provided as an n column data, where
  #   n is the number of dimensions

  d <- vector(mode="numeric", length=nrow(destinations))
  
  for (i in 1:dimensions) {
    source.coord <- as.numeric(source[1, j=i, with=FALSE])
    d <- d + (destinations[, j=i, with=FALSE] - source.coord) ^ 2
  }
  distances <- sqrt(d)

  return(distances)
}

prob.geneflow <- function(gene.dist, threshold=0.001, zero_flow_dist=5) {
  A <- 14
  B <- 0.5
  
  # first a fairly standard logistic function
  logistic_result <- (1 / (1 + exp(((gene.dist / zero_flow_dist) - B) * A)))
  
  # then expand on the y axis (while remaining centred at 0.5)
  # so that the values of 0 and 1 are actually reached
  rescale_y <- (1 / (1 - (2 * threshold))) 
  prob <-  (logistic_result * rescale_y) - threshold
  return(prob)
}

is.geneflow <- function(demetable.cell, deme.primary, dimensions=3, speciation.gene.distance) {
  # this function measures the genetic distance between the primary deme and the
  # other origin demes and applies a distance function to determine if there is
  # gene flow.  The result is returned as TRUE or FALSE for each origin deme.

  first.gene.col.idx    <- which(names(demetable.cell)=="gene.pos1")
  gene.cols             <- first.gene.col.idx:(first.gene.col.idx + dimensions - 1)
  source.coords         <- demetable.cell[deme.primary, gene.cols, with=FALSE]
  destination.coords    <- demetable.cell[, gene.cols, with=FALSE]
  g.distances           <- gene.distances(source.coords, destination.coords, dimensions)
  demetable.cell$probs  <- prob.geneflow(g.distances, zero_flow_dist=speciation.gene.distance)
  demetable.cell[probs>1, "probs"] <- 1

  # determine geneflow using the probablities
  origin.deme.count     <- nrow(demetable.cell)
  for (k in 1:origin.deme.count) {
    demetable.cell$geneflow[k] <- sample(x=c(1,0), size=1, replace=TRUE, prob=c(demetable.cell$probs[k],
                                                                                demetable.cell$probs[k]))
  }
  demetable.cell$geneflow <- demetable.cell$geneflow == 1
  return(demetable.cell$geneflow)

}