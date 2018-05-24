######### speciateAllopatric #########

# Function divides a species into two daughter species based on an allopatric model of speciation

############# Arguments ###############

# env = current environmental layer
# position = niche position of species 
# breadth = niche breadth of species
# phylo.sig = phylogenetic signal of species niches at speciation
# species.raster is the current species range in raster format (presence/absence)

speciateAllopatric <- function (position, breadth, species.raster, env, phylo.sig){
  # get coordinates of species range
  ras1 <- species.raster
  coords <- xyFromCell(ras1, Which(ras1 == 1, cells=T))
  x <- range(coords[,1])
  y <- range(coords[,2])
  #  range bisection occurs across x or y axis of species range (latitudinaly or longitudinaly)
  x.or.y <- sample(c("x", "y"),1)
  # select two points along axis to draw line of bisection
  if(x.or.y == "x"){
    x <- runif(2, min=x[1], max=x[2])
  } else {
    y <- runif(2, min=y[1], max=y[2])
  }
  # make points spatial coordinates
  xy <- cbind(x,y)
  xy.sp = sp::SpatialPoints(xy)
  # create spatial line linking two spatial coordinates
  spl <- SpatialLines(list(Lines(Line(xy.sp), ID="a")))
  # find coordinates of the grid cells that fall along line
  cells <- extract(ras1, spl, cellnumbers=T)
  xy <- xyFromCell(ras1, cells[[1]][,1])
  #determine which cells fall on one side of the bisection line
  if(x.or.y =="x"){
    if(x[1] > x[2]){
      sp1.coords <- list()
      for (i in 1:length(xy[,1])){
        sp1.coords[[i]] <- coords[which(coords[,1] > xy[i,1] & coords[,2] > xy[i,2]),]  
      }
    } else {
      sp1.coords <- list()
      for (i in 1:length(xy[,1])){
        sp1.coords[[i]] <- coords[which(coords[,1] > xy[i,1] & coords[,2] < xy[i,2]),]  
      }
    }
  } else {
    if(y[1] > y[2]){
      sp1.coords <- list()
      for (i in 1:length(xy[,1])){
        sp1.coords[[i]] <- coords[which(coords[,1] > xy[i,1] & coords[,2] > xy[i,2]),]  
      }
    } else {
      sp1.coords <- list()
      for (i in 1:length(xy[,1])){
        sp1.coords[[i]] <- coords[which(coords[,1] > xy[i,1] & coords[,2] < xy[i,2]),]  
      }
    } 
  }
  # these cells form the basis for the first daughter species and the second daughter species occupies the rest of the cells
  sp1.df <- do.call("rbind", sp1.coords)
  sp1.cells <- extract(ras1, sp1.df, cellnumbers=T)[,1]
  ras1[sp1.cells] <- 2
  sp1 <- ras1
  sp1@data@values[which(sp1@data@values==1)] <- NA
  sp1@data@values[which(sp1@data@values == 2)] <- 1
  sp2 <- ras1
  sp2@data@values[which(sp2@data@values==2)] <- NA
  #each species niche and range is recentered based on the degree of phylogenetic signal in the niche (phylo.sig)
  recenteredsp1 <- nicheRecenter(position, breadth, sp1, env, phylo.sig)
  recenteredsp2 <- nicheRecenter(position, breadth, sp2, env, phylo.sig)
  # returns two daughter species' ranges, niche breadths and niche positions
  result <- list(positions = c(recenteredsp1[[1]],recenteredsp2[[1]]), breadths = c(recenteredsp1[[2]],recenteredsp2[[2]]), species.rasters=c(recenteredsp1[[3]],recenteredsp2[[3]]))
  return(result)
}
