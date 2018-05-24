######### nicheEvolution #########

# evolves both the niche position and nice bredath of a species according to modified Brownian motion model

############# Arguments ###############

# env = current environmental layer
# position = niche position of species 
# breadth = niche breadth of species
# pos.ev.rate = niche position evolution rate
# breadth ev. rate = niche breadth evolution rate
# min.breadth = minimum niche breadth permitted
# max.breadth = maximum niche breadth permitted (set to 10 which means a species can occpuy 80% of environmental conditions in the domain at the widest niche breadth)
# species.raster is the current species range in raster format (presence/absence)

nicheEvolution <- function(position, breadth, pos.ev.rate, breadth.ev.rate, species.raster, env, min.breadth=0.1, max.breadth=10){
  orig.position <- position
  orig.breadth <- breadth
  env2 <- envRasterFromDistribution(species.raster, env)
  #new position and breadth adds a random deviate to original niche values
  position <- position + rnorm(1, mean = 0, sd = pos.ev.rate)
  breadth  <- breadth + rnorm(1, mean = 0, sd = breadth.ev.rate)
  # condition that if new breadth extends beyond bounds then breadth = the max or min
  if(breadth < min.breadth) breadth <- min.breadth
  if(breadth > max.breadth) breadth <- max.breadth
  # set niche bounds
  upper <- position + breadth
  lower <- position - breadth
  # if all the environmental values within the species range are greater than their niche maximum or lower than their niche minimum need to re-evolve niche
  # if species have very small ranges, small niche breadths there may be very little room for the niche to shift without removing all suitable habitat from within the species range
  # for this reaon we include the condition to break after 100 times
  if(all(env2@data@values > upper | env2@data@values < lower, na.rm=T)){
    loops <- 1
    while(all(env2@data@values > upper | env2@data@values < lower, na.rm=T)) {
      position <- orig.position + rnorm(1, mean = 0, sd = pos.ev.rate)
      breadth  <- orig.breadth + rnorm(1, mean = 0, sd = breadth.ev.rate)
      if(breadth < min.breadth) breadth <- min.breadth
      if(breadth > max.breadth) breadth <- max.breadth
      upper <- position + breadth
      lower <- position - breadth 
      loops <- loops + 1
      if(loops == 100) {
        position <- orig.position
        breadth <- orig.breadth
        upper <- position + breadth;
        lower <- position - breadth
        break
      }
    }
  }
  # remove values from within species range that do no longer fall within their new niche breadth
  species.raster[which(env2@data@values > upper | env2@data@values < lower)] <- NA
  # return new species range, niche position, and niche breadth
  return(c(species.raster, position, breadth))  
} 