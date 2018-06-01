
scripts <- c("rangeDispersal.R", "nicheEvolution.R","speciateAllopatric.R","speciateSympatric.R",
             "speciateParapatric.R","speciateDispersal.R","seedSpecies.R","environmentalChange.R",
             "nicheRecenter.R","DREaD_ds.R","generateSummaryStatistics.R", "helperFunctions.R", "dataStructures.r")

lapply(scripts, source)

#required.packages <- (c("raster","gstat","ape","phytools","geiger","phyloclim","ggplot2","gridExtra","moments",
#                        "apTreeshape","parallel", "doSNOW", "rgeos", "data.table", "fossil", "ENMTools"))
#lapply(required.packages, require, character.only=T)


#sample parameters
total.time      <- 100
D               <- runif(1, 1, 10)
ENVa            <- runif(1, 0.25, 2)
ENVf            <- runif(1, 0.25, 2)
NEp             <- runif(1, 0.005, 2)
NEb             <- runif(1, 0.0025, 1)
PS              <- runif(1, 0.25, 1)
m               <- runif(1, 50, 250)
niche.blocksize <- 0.1

# Run model
simulation.1 <- DREaD_ds(total.time=total.time, dispersal=D, amp=ENVa, freq=ENVf,
                  niche.ev.rate=NEp, breadth.ev.rate=NEb, phylo.sig=PS, Me=m,
                  enviro.hetero=T, geo.mode="dispersal", enviro.mode="sine")
