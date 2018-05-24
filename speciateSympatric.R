######### speciateSympatric #########

# Function divides a species into two daughter species based on an sympatric model of speciation

############# Arguments ###############

# env = current environmental layer
# position = niche position of species 
# breadth = niche breadth of species
# phylo.sig = phylogenetic signal of species niches at speciation
# species.raster is the current species range in raster format (presence/absence)

speciateSympatric <- function (position, breadth, species.raster, env, phylo.sig){
  #obtain x and y coordinates of current species range
  ras <- species.raster
  occupied.space <- rasterToPoints(ras)
  #sample 2 x and 2 y coordinates to form the range boundaries of the new daughter species geographic range from within the original species range
  xcoords <- sample(unique(occupied.space[,1]), 2, replace=F)
  ycoords <- sample(unique(occupied.space[,2]), 2, replace=F)
  # conditions to make sure the north boundary is south of the south boundary etc.
  new.range <- matrix(, ncol=2, nrow=2)
  if (xcoords[1] > xcoords[2]) {
    new.range[1,2] <- xcoords[1]
    new.range[1,1] <- xcoords[2]
  } else {
    new.range[1,2] <- xcoords[2]
    new.range[1,1] <- xcoords[1] 
  }
  if(ycoords[1] > ycoords[2]) {
    new.range[2,2] <- ycoords[1]
    new.range[2,1] <- ycoords[2]
  } else {
    new.range[2,2] <- ycoords[2]
    new.range[2,1] <- ycoords[1]  
  }
  # new species range
  new.range<-extent(new.range) 
  new.ras <- crop(ras, new.range, snap="out")
  # check that the new species range falls within the range of the roiginal species
  # resample species boundaries untill have a completely sympatric daughter species
  while(any(is.na(new.ras@data@values))){
    xcoords <- sample(unique(occupied.space[,1]), 2, replace=F)
    ycoords <- sample(unique(occupied.space[,2]), 2, replace=F)
    new.range <- matrix(, ncol=2, nrow=2)
    if (xcoords[1] > xcoords[2]) {
      new.range[1,2] <- xcoords[1]
      new.range[1,1] <- xcoords[2]
    } else {
      new.range[1,2] <- xcoords[2]
      new.range[1,1] <- xcoords[1] 
    }
    if(ycoords[1] > ycoords[2]) {
      new.range[2,2] <- ycoords[1]
      new.range[2,1] <- ycoords[2]
    } else {
      new.range[2,2] <- ycoords[2]
      new.range[2,1] <- ycoords[1]  
    }
    new.range<-extent(new.range) 
    new.ras <- crop(ras, new.range, snap="out")  
  }
  new.ras<-extend(new.ras, env)
  # new daughter species niche and range is recentered based on the degree of phylogenetic signal in the niche (phylo.sig)
  recenteredsp2 <- nicheRecenter(position, breadth, new.ras, env, phylo.sig)
  # returns new daughter species' range, niche breadth and niche position (other daughter species retains the values of the parent species)
  result <- list(positions = recenteredsp2[[1]], breadths = recenteredsp2[[2]], species.rasters=recenteredsp2[[3]])
  return(result)
} 

