rm(list=ls())
gc()

scripts <- c("disperse_evolve_compete_ds.r","seedSpecies.R","environmentalChange.R", "DREaD_ds.R",
             "generateSummaryStatistics.R", "helperFunctions.R", "dataStructures.r")

lapply(scripts, source)

#required.packages <- (c("raster","gstat","ape","phytools","geiger","phyloclim","ggplot2","gridExtra","moments",
#                        "apTreeshape","parallel", "doSNOW", "rgeos", "data.table", "fossil", "ENMTools"))
#lapply(required.packages, require, character.only=T)


#sample parameters
total.time            <- 100
D                     <- runif(1, 1, 10)
ENVa                  <- runif(1, 0.25, 2)
ENVf                  <- runif(1, 0.25, 2)
NEp                   <- runif(1, 0.005, 2)
NEb                   <- runif(1, 0.0025, 1)
PS                    <- runif(1, 0.25, 1)
m                     <- runif(1, 50, 250)
niche.blocksize       <- 0.1
suitability.mode      <- "sine"
speciation.gene.distance <- 10  # this parameter will need to be set with the drift rate
niche.evolution.rate  <- 0.2  # this (so far) is just the proportion of the gap between the min (or max)
                          # and the local environment is reduced each timestep
do.display            <- FALSE
do.text.output        <- TRUE

source("dynamicDisplay.r")

if (do.display) {
  my.colours <- display.initialise()
}

# Run model
simulation.1 <- DREaD_ds(total.time=total.time, dispersal=D, amp=ENVa, freq=ENVf,
                  niche.ev.rate=NEp, breadth.ev.rate=NEb, phylo.sig=PS, Me=m,
                  enviro.hetero=T, geo.mode="dispersal", enviro.mode="sine", suitability.mode=suitability.mode,
                  speciation.gene.distance=speciation.gene.distance, niche.evolution.rate=niche.evolution.rate)
