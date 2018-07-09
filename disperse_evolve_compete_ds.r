# this file contains functions for dispersal, competition, niche evolution and drift

# env = current environmental layer
# position = niche position of species
# breadth = niche breadth of species

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

  # get maximum x and y for simulation area
  domain.max <- env.table[,.(xmax=max(col), ymax=max(row))]

  # create an empty demetable.species.new
  demetable.species.new <- cbind(demetable.species[0, 1:12], data.table(originCell=0)[0])

  # loop through the niche groups
  for (k in 1:nrow(niche.groups)) {
    niche <- niche.groups[k, ]
    niche.values    <- niche * niche.blocksize  # turns the rounded niche group into niche values

    demetable.nichegroup <- demetable.species[niche1.breadth.group == niche$niche1.breadth.group &
                                    niche1.position.group== niche$niche1.position.group &
                                    niche2.breadth.group == niche$niche2.breadth.group &
                                    niche2.position.group== niche$niche2.position.group, 1:12]

    # filter env cells to those within the spatial limits (for the niche group)
    bounds        <- as.list(demetable.nichegroup[, list(xmin=(min(x)-dispersal.range), xmax=(max(x)+dispersal.range), ymin=(min(y)-dispersal.range), ymax=(max(y)+dispersal.range))])

    bounds[bounds < 0] <- 0
    if(bounds$xmax > domain.max$xmax) {bounds$xmax <- domain.max$xmax}
    if(bounds$ymax > domain.max$ymax) {bounds$ymax <- domain.max$ymax}

    env.table.dispersal <- env.table[col >= bounds$xmin &
                                      col <= bounds$xmax &
                                      row >= bounds$ymin &
                                      row <= bounds$ymax]

    # filter env cells to those within the niche limits
    niche1.min <- niche.values$niche1.position.group - niche.values$niche1.breadth / 2
    niche1.max <- niche.values$niche1.position.group + niche.values$niche1.breadth / 2
    #niche2.min <- niche.values$niche2.position.group - niche.values$niche2.breadth / 2
    #niche2.max <- niche.values$niche2.position.group + niche.values$niche2.breadth / 2

    env.table.dispersal <- env.table.dispersal[env1 >= niche1.min & env1 <= niche1.max, ]  # need to add 2nd niche dimension

    if (env.table.dispersal[, .N] > 0) {

      #apply niche suitability function to env.dispersal to give suitability for 0 to 1
      env.table.dispersal$suitability <- niche_suitability(env = env.table.dispersal$env1,
                                                           niche1.breadth = niche.values$niche1.breadth,
                                                           niche1.position = niche.values$niche1.position,
                                                           suitability.mode = suitability.mode)

      #remove dispersal destinations which are within the bounding rectangle, but too far away
      include.in.dispersal <- close.enough(demetable.nichegroup[, c("x", "y")],
                                           env.table.dispersal[, c("col", "row")],
                                           dispersal.range)
      env.table.dispersal <- env.table.dispersal[include.in.dispersal]
    }

    # process this nichegroup only if there are cells to disperse to
    if (nrow(env.table.dispersal) > 0) {

      demetable.nichegroup.new <- demetable.nichegroup[0, ]
      row.pointer <- 0

      # loop through each deme for the niche group
      for (d in 1:nrow(demetable.nichegroup)) {

        deme <-  demetable.nichegroup[d,]

        # remove once dispersal is working correctly
        # if (do.display) {
        #   display.update(list(one_deme=deme))
        # }

        # find the cells in dispersal distance
        deme.dest <- env.table.dispersal[(col >= deme$x-dispersal.range
                                         & col <= deme$x+dispersal.range
                                         & row >= deme$y-dispersal.range
                                         & row <= deme$y+dispersal.range)]
        new.count <- nrow(deme.dest)
        new.amount  <- deme$amount * deme.dest$suitability  # should return a vector
        new.rows    <- (row.pointer + 1):(row.pointer + new.count)

        if (length(new.amount) > 0) {
          # create dispersed demes

          demetable.nichegroup.new <- rbind(demetable.nichegroup.new,
                                            list(cellID = deme.dest$cellID,
                                                 x =      deme.dest$col,
                                                 y =      deme.dest$row,
                                                 amount = new.amount), fill = TRUE)

          set(demetable.nichegroup.new,  # assign values which are the same for each new row in the niche group
              i=new.rows,
              j=c("speciesID", "niche1.position", "niche1.breadth", "niche2.position",
                  "niche2.breadth", "gene.pos1", "gene.pos2", "gene.pos3", "originCell"),
              value=list(deme$speciesID, deme$niche1.position, deme$niche1.breadth, deme$niche2.position,
                         deme$niche2.breadth, deme$gene.pos1, deme$gene.pos2, deme$gene.pos3, deme$cellID)
              )

          # apply distance function
          weights <- distance.weights(source = deme[1, c("x", "y")],
                                      destinations = demetable.nichegroup.new[new.rows, c("x", "y")],
                                      dispersal.range = dispersal.range,
                                      distance.type = "euclidean",
                                      distance.decay = "linear")

          demetable.nichegroup.new$amount[new.rows] <- demetable.nichegroup.new$amount[new.rows] * weights
          new.count <- length(which(demetable.nichegroup.new$amount[new.rows] > 0))
          demetable.nichegroup.new <- demetable.nichegroup.new[amount > 0]

        }

        row.pointer <- row.pointer + new.count

      }

      # combine the niche group demetable rows
      demetable.species.new <- rbind(demetable.species.new, demetable.nichegroup.new)
    }
  }

  return(demetable.species.new)

}

combine.demes <- function (demetable.species.overlap,
                           genomeDimensions,
                           speciation.gene.distance,
                           minimum.amount,
                           env.table,
                           verbose=FALSE) {
  # determine what happens to populations which are in the same cell

  demetable.species <- demetable.species.overlap[0, -"originCell"]

  deme.cells <- unique(demetable.species.overlap$cellID)
  for (cell in deme.cells) {
    demetable.cell <- demetable.species.overlap[cellID==cell]
    deme.count     <- demetable.species.overlap[cellID==cell, .N]

    if (deme.count==1) {
      demetable.species <- rbindlist(list(demetable.species, demetable.cell[, -"originCell"]))
    } else {

      if (verbose) {
        cat("For species", demetable.species.overlap[1, speciesID], "in cell", cell, "there are", deme.count, "origin demes to combine\n")
      }
      # determine the primary deme, which the others may join with
      # if the current cell is one of the sources, it's the primary cell
      deme.primary <- which(demetable.cell$originCell==cell)

      if (length(deme.primary) == 0) {
        # otherwise pick a cell at random, with probability weighted by amount
        deme.primary <- sample(x=deme.count, size=1, prob=demetable.cell$amount, replace = T)
      }

      # now loop through each of the source demes and combine them with the primary deme, as determined above
      demetable.cell <- demetable.cell[is.geneflow(demetable.cell, deme.primary, dimensions=3, speciation.gene.distance)]
      set(demetable.cell, j="gene.flow", value=TRUE)

      if (nrow(demetable.cell) > 0) { # make sure that there are still demes with gene flow - better check why not - should be gene flow to self
        demetable.species.new <- demetable.cell[, list(amount=1,  # this amount is just a place filler. Real amount will be the suitability for the new niche
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

            demetable.species.new[[column.name]] <- weighted.mean(demetable.cell[[column.name]], demetable.cell[["amount"]])
          }
        }

        demetable.species <- rbindlist(list(demetable.species, demetable.species.new))
      }
    }

  }

  # update amount based on environment and new niche
if (nrow(demetable.species) > 0) {
  setkey(demetable.species, cellID)

  env.demes <- merge(env.table, demetable.species, by="cellID", nomatch=0)
  demetable.species$amount <- env.demes[, .(amount=niche_suitability(env1, niche1.position, niche1.breadth, niche2.position, niche2.breadth, suitability.mode)), by=cellID][, 2]
  }

  # remove occurrences below the amount threshold, after combinations are complete
  demetable.species <- subset(demetable.species, amount >= minimum.amount)

  return(demetable.species)
}

niche_suitability <- function(env,
                               niche1.position,
                               niche1.breadth,
                               niche2.position=0,
                               niche2.breadth=0,
                               suitability.mode="sine") {
  # convert environment to suitability

  niche1.min <- niche1.position - (niche1.breadth /2)
  niche1.max <- niche1.position + (niche1.breadth /2)
  niche2.min <- niche2.position - (niche2.breadth /2)
  niche2.max <- niche2.position + (niche2.breadth /2)

  env.class <- class(env)  # allow env to be a numeric vector, a raster or a data.table, but handle it as a vector

  if (env.class=="raster") {
    suitability <- env[]
  } else if (any(env.class=="data.table")) {
    suitability <- env[,1] # the first column must be env1.
    # Note that this does not yet handle a second environment dimension
  } else {
    suitability <- env
  }

  if (suitability.mode=="block") {
    # block suitability is 1 within the suitable range, zero elsewhere
    suitability[suitability < niche1.min | suitability > niche1.max] <- 0
    suitability[suitability > 0] <- 1

  } else if (suitability.mode=="sine") {
    # sine suitability is 1 for the central niche position, declining to 0 at
    # the niche limits according to a sine curve

    # transform suitability to range from 0 (for minimum to pi for maximum)
    suitability <- suitability - niche1.min
    suitability[suitability < 0] <- 0
    suitability <- suitability * (pi / niche1.breadth)
    suitability[suitability > pi] <- 0
    suitability <- sin(suitability)
  }

  return(suitability)

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

    distance.weights  <- round((dispersal.range - (distances * (1-min.weight))) / dispersal.range, 5)
    # rounding is to prevent unintended results due to numerical precision issues
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
    d <- d + ((destinations[, j=i, with=FALSE] - source.coord) ^ 2)
  }

  return(sqrt(d))
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
  # this function calculates the genetic distance between the primary deme and the
  # other origin demes and applies a distance function to determine if there is
  # gene flow.  The result is returned as TRUE or FALSE for each origin deme.

  first.gene.col.idx    <- which(names(demetable.cell)=="gene.pos1")
  gene.cols             <- first.gene.col.idx:(first.gene.col.idx + dimensions - 1)
  source.coords         <- demetable.cell[deme.primary, gene.cols, with=FALSE]
  destination.coords    <- demetable.cell[, gene.cols, with=FALSE]

  set(demetable.cell, j="geneflow", value=FALSE)

  # deal quickly with those that have the same location
  for (dest.row in 1:nrow(destination.coords)) {
    if (all(destination.coords[dest.row]==source.coords)) {
      set(demetable.cell, dest.row, "geneflow", TRUE)
    }
  }

  dest.rows <- which(!demetable.cell$geneflow)

  if (length(dest.rows) > 0) {  # do this loop where the destination gene position differs from the origin
    g.distances   <- gene.distances(source.coords, destination.coords[dest.rows], dimensions)
    set(demetable.cell, i=dest.rows, j="probs", value=prob.geneflow(g.distances, zero_flow_dist=speciation.gene.distance))
    demetable.cell <- demetable.cell[(probs>1), probs := 1]

    # determine geneflow using the probablities
    for (k in dest.rows) {
      demetable.cell$geneflow[k] <- (sample(x=c(1,0),
                                            size=1,
                                            replace=TRUE,
                                            prob=c(demetable.cell$probs[k], (1-demetable.cell$probs[k])))
                                     ==1)
    }
  }

  return(demetable.cell$geneflow)

}

close.enough <- function(dispersal.origins, dispersal.destinations, dispersal.range) {
  # this function returns the index of dispersal destinations which are close enough to include
  # in the dispersal calculation

  # dispersal.origins should be a 2 column data.table of x,y
  # dispersal.destinations should be a 3 column data.table of index,x,y with index order maintained
  #   to return a subset which is close enough, as an index vector
  # dispersal.range is the maximum dispersal distance

  # check for valid rows
  if (nrow(dispersal.destinations) > 0 & nrow(dispersal.origins) > 0) {

    dispersal.destinations <- cbind(1:nrow(dispersal.destinations), dispersal.destinations, include=0)
    names(dispersal.destinations)[1:3] <- c("rownum", "x", "y")

    on.origin <- dispersal.destinations[dispersal.origins, on=c(x="x", y="y"), nomatch=0][, 1]
    dispersal.destinations$include[on.origin$rownum] <- 1

    not.on.origin <- which(dispersal.destinations$include==0)

    for (destination.index in not.on.origin) {
      xy <- dispersal.destinations[destination.index, 2:3]
      dist.from.dest <- pointDistance(xy, dispersal.origins, lonlat = FALSE)
      if(min(dist.from.dest) <= dispersal.range) {
        set(dispersal.destinations, destination.index, "include", 1)
      }
    }

  } else {
    cat("\nFunction close.enough() requires at least one dispersal.destinations row and one dispersal.origins row.\n")
browser()
  }

  return(which(dispersal.destinations$include==1))

}

niche.evolution <- function(demetable.species,
                            env.table,
                            niche.evolution.rate) {

  # this function adjusts the niche limits of each deme towards the local environment
  for (d in 1:nrow(demetable.species)) {
    deme <- demetable.species[d]

    niche1.min <- deme$niche1.position - (deme$niche1.breadth / 2)
    niche1.max <- deme$niche1.position + (deme$niche1.breadth / 2)

    env1 <- env.table$env1[env.table$cellID == deme$cellID]

    # move the niche min and max towards the local environment
    niche1.min.new <- niche1.min + ((env1 - niche1.min) * niche.evolution.rate)
    niche1.max.new <- niche1.max - ((niche1.max - env1) * niche.evolution.rate)

    demetable.species[d, "niche1.breadth"]  <- niche1.max.new - niche1.min.new
    demetable.species[d, "niche1.position"] <- (niche1.max.new + niche1.min.new) / 2
  }
  return(demetable.species)
}

extinction <- function(edgetable, speciesID, current.time) {
  # this function takes the steps required when a species is found to have no demes
  # initially it makes the most minimal changes, but more could be done here

  cat("\nSpecies", speciesID, "has become extinct at time", current.time, "\n")
  browser()
}