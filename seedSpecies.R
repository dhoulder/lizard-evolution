######### seedSpecies #########

# Function generate the first species in the simulator (the ancestor of the simulated clade)
# this species has a geographic range in raster format with presence = 1 and absence = NA
# the species also has a niche position and niche breadth

############# Arguments ###############

# env = current environmental layer
# dispersal = dispersal kernal of the species 
# resolution of the environmental layer. defaults to 1 as this is the resolution of the simulated environment

seedSpecies <- function(env, dispersal, resolution=1) {
# sample a single cell from the domain  
  cell<-sample(1:length(env), 1) 
#checks the values of the cell is not NA (may be the case with imported rasters)  
  while(is.na(env[cell])==TRUE){cell<-sample(1:length(env), 1) }
# starting niche position is the value of the selected cell
    starting.cell <- xyFromCell(env, cell)
    initial.position <- extract(env, starting.cell)
# niche breadth is sampled from a uniform distribution (0.1 to 10)    
    initial.breadth <- runif(1, min=0.1, max=10)
# range boundaries are sampled from a uniform distribution (from the value of the dispersal kernal to 25)
    initial.extent <- runif(4, min=dispersal, max=25)
# range boundaries are drawn around the starting cell
    initial.extent <- extent(c(starting.cell[[1]] - initial.extent[[1]], starting.cell[[1]] + initial.extent[[2]], starting.cell[[2]] - initial.extent[[3]], starting.cell[[2]] + initial.extent[[4]])) 
# conditions to ensure that the range does not extend beyond the boundaries of the domain
    if (initial.extent[1] < env@extent[1])  {initial.extent[1] <- env@extent[1]}
    if (initial.extent[2] > env@extent[2]) {initial.extent[2] <- env@extent[2]}
    if (initial.extent[3] < env@extent[3]) {initial.extent[3] <- env@extent[3]}
    if (initial.extent[4] > env@extent[4]) {initial.extent[4] <- env@extent[4]}
# remove cells from the range that do not fall within the species niche    
initial.ras <- crop(env, initial.extent)
upper <- initial.position + initial.breadth
lower <- initial.position - initial.breadth
initial.ras[initial.ras@data@values > upper] <- NA
initial.ras[initial.ras@data@values < lower]  <- NA
initial.ras[!is.na(initial.ras@data@values)] <- 1
initial.ras<-extend(initial.ras, env)
# return species range, niche position, and niche breadth
return(c(initial.ras, initial.position, initial.breadth))
}