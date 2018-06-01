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