# this file contains functions for dispersal, competition, niche evolution and drift


######### rangeDispersal #########

# Function extends a species range based on a dispersal kernal and the inherent niche breadth of the species

############# Arguments ###############

# env = current environmental layer
# position = niche position of species 
# breadth = niche breadth of species
# species.ras is the current species range in raster format (presence/absence)

disperse_ds <- function (demetable_sp, env, dispersal.range){
  
  # group demes by similarity of the niche - so for very similar niches the environment is calculated just once
  
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