# this file contains functions for dispersal, competition, niche evolution and drift


######### rangeDispersal #########

# Function extends a species range based on a dispersal kernal and the inherent niche breadth of the species

############# Arguments ###############

# env = current environmental layer
# position = niche position of species 
# breadth = niche breadth of species
# species.ras is the current species range in raster format (presence/absence)

disperse_ds <- function (demetable.species, env, dispersal.range){
  
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
    demetable.nichegroup <- demetable.species[niche1.breadth.group == niche$niche1.breadth.group & 
                                    niche1.position.group== niche$niche1.position.group &
                                    niche2.breadth.group == niche$niche2.breadth.group &
                                    niche2.position.group== niche$niche2.position.group, 1:9]
    
    # filter env cells to those within the spatial limits (for the niche group) 
    bounds        <- as.list(demetable.species[, list(xmin=(min(x)-dispersal.range), xmax=(max(x)+dispersal.range), ymin=(min(y)-dispersal.range), ymax=(max(y)+dispersal.range))])
    niche.extent  <- extent(bounds[[1]], bounds[[2]], bounds[[3]], bounds[[4]])
    
    env.dispersal <- crop(env, niche.extent)
    
    # filter env cells to those within the niche limits
    niche.values    <- niche * niche.blocksize  # turns the rounded niche group into niche values
    niche1.min <- niche.values$niche1.position.group - niche.values$niche1.breadth /2
    niche1.max <- niche.values$niche1.position.group + niche.values$niche1.breadth /2    
    niche2.min <- niche.values$niche2.position.group - niche.values$niche2.breadth /2
    niche2.max <- niche.values$niche2.position.group + niche.values$niche2.breadth /2    
    
    env.dispersal[env.dispersal[] < niche1.min | env.dispersal[] > niche1.max] <- NA
    
    #apply niche suitability function to env.dispersal to give suitability for 0 to 1
    
    dispersal.cells <- which(!is.na(env.dispersal[]))
    suitability <- niche_suitability(env, suitability.mode = "sine", niche1.breadth = niche.values$niche1.breadth, niche1.position = niche.values$niche1.position)
      
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
  # convert environment raster to a suitability raster

  niche1.min <- niche1.position - (niche1.breadth /2)
  niche1.max <- niche1.position + (niche1.breadth /2)
  niche2.min <- niche2.position - (niche2.breadth /2)
  niche2.max <- niche2.position + (niche2.breadth /2)
  suitability <- env

  if (suitability.mode=="block") {
    
    suitability[suitability < niche1.min | suitability > niche1.max] <- 0
    suitability[suitability > 0] <- 1
    
  } else if (suitability.mode=="sine") {
    
    # transform suitability to range from 0 (for minimum to pi for maximum)
    suitability <- suitability - niche1.min
    suitability[suitability < 0] <- 0
    suitability <- suitability * (pi / niche1.breadth)
    suitability[suitability > pi] <- 0
    suitability <- sin(suitability)
  }
  
  return(suitability)
}