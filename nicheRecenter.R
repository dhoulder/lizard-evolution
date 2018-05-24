######### nicheRecenter #########

# Function shifts the niche towards the mean values of the species range at speciation based on the degree of phylogenetic signal in species ranges
# high phylogenetic signal and species will maintain a niche more similar to their parents ncihe
# low phylogenetic signal means species will move further towards the mean of their new range

############# Arguments ###############

# env = current environmental layer
# position = niche position of species 
# breadth = niche breadth of species
# phylo.sig = phylogenetic signal of species niches at speciation
# species.raster is the current species range in raster format (presence/absence)

nicheRecenter <- function (position, breadth, species.raster, env, phylo.sig) {
# new psotion moves towards the mean of the new geographic range  
new.pos <- position - ((position - mean(env@data@values[Which(species.raster, cells=T)], na.rm=T)) * phylo.sig)  
# breadth shifts as per nicheEvolution function
new.breadth <- breadth + rnorm(1, mean = 0, sd = phylo.sig)
env.raster <- envRasterFromDistribution(species.raster, env)
env.raster@data@values[which(env.raster@data@values > new.pos+new.breadth | env.raster@data@values < new.pos-new.breadth)] <- NA
new.range <- replaceNonNAValues(env.raster, 1)
return(c(new.pos, new.breadth, new.range)) 
}
