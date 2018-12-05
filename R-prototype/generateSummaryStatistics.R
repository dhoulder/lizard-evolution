######### generateSummaryStatistics #########

# Function to generate the summary statsitics required for data analysis and ABC


############# Arguments ###############

# df = edgetable from simulation
# rasters = species rasters list
# phy = phylogeny
# drop.fossil = logical to drop extinct taxa from analysis


generateSummaryStatistics <- function (df, rasters, phy, drop.fossil=T) {

  if(drop.fossil==T){ phy <- drop.fossil(phy) }
  if(class(rasters)=="RasterStack") {rasters <- unstack(rasters)}
  # reduce data.frame to just extant species to match phylogeny
  df.dropped <- df[which(df[,2] %in% phy$tip.label),]
  df.dropped <- df.dropped[order(as.numeric(df.dropped[,"rasterMatch"]), decreasing=FALSE),]
  # ID sister species pairs
  sisters <- findSisters(phy, solve.polytomies=T)
  sisters.rasters <- list()
  df<-df.dropped
  # match species rasters numbers in the df to their position in rasters
  df[,"rasterMatch"]<-c(1:length(rasters))

# Find sister species's raster numbers and get divergence times of sister species
  sisters.split.time<-c()
  root.age.true <- max(findMRCA(phy, type=c("height")))
  # fix root age to 1 so branch lengths are in relative time for comparison between replicates
  phy$edge.length<-phy$edge.length/root.age.true
  root.age<-1
  for ( i in 1:length(sisters)) {
    raster.nums<-c()  
    sisters.split.time[i]<-(root.age-findMRCA(phy, tips=paste(sisters[[i]]), type=c("height")))
    for (j in 1:length(sisters[[i]])) {
      raster.nums[j] <- df[which(df[,2] == sisters[[i]][j]), "rasterMatch"]
    }
    sisters.rasters[[length(sisters.rasters) + 1]] <- raster.nums
  }

# Get range overlap of sister pairs
# range overlap is calculated as the overlap of the larger ranged sister over the smaller ranged sister 
  sisters.overlap<-list()
  for (i in 1:length(sisters.rasters)) {
    sis.1 <- rasters[[as.numeric(sisters.rasters[[i]][1])]]
    sis.2 <- rasters[[as.numeric(sisters.rasters[[i]][2])]]
    sis.1.area <- sum(sis.1@data@values, na.rm=T)
    sis.2.area <- sum(sis.2@data@values, na.rm=T)
    sis.overlap<-sis.1+sis.2
    if(any(sis.overlap@data@values==2, na.rm=T)) {
      if (sis.1.area <= sis.2.area) {
        range.overlap <- length(which(sis.overlap@data@values == 2))/ sis.1.area
      }
      if (sis.1.area > sis.2.area) {
       range.overlap <- length(which(sis.overlap@data@values == 2))/ sis.2.area
      }
    } else { range.overlap <- 0 }
    sisters.overlap[[length(sisters.overlap) + 1]] <- range.overlap
  }

# Get range size asymmetry and range isolation of sister pairs
# range size asymmetry is calculated as the size of the smaller ranged species / the sum of both species range sizes
# range isolation is calculated as the minimum distance between the sister species' ranges divided by the total possible distance between species (a constant for comparison with empirical data)  
  sister.asymmetry <- c()
  range.distance <- c()
  for (i in 1:length(sisters.rasters)) {
    sis.1 <- rasters[[as.numeric(sisters.rasters[[i]][1])]]
    sis.2 <- rasters[[as.numeric(sisters.rasters[[i]][2])]]
    sis.1.size <- sum(sis.1@data@values, na.rm=T)
    sis.2.size <- sum(sis.2@data@values, na.rm=T)
    p1 <- data.frame(rasterToPoints(sis.1))[,1:2]
    p2 <- data.frame(rasterToPoints(sis.2))[,1:2]
    range.distance[i]<- min(pointDistance(p1, p2, longlat=F))/140.0071
    if(sis.1.size < sis.2.size) {
    sister.asymmetry[i] <-sis.1.size / (sis.1.size + sis.2.size)
  } else {
    sister.asymmetry[i] <-sis.2.size / (sis.2.size + sis.1.size)
  }
}
# calculate mean and SD of range distances
  range.distance.av <- mean(range.distance, na.rm=T)
  range.distance.sd <- sd(range.distance, na.rm=T)
  
# Linear regression of asymmetry, overlap, and distance against node ages for sister species
  lm.asym <- lm(sister.asymmetry~sisters.split.time)
  lm.overlap <- lm(unlist(sisters.overlap)~sisters.split.time)
  lm.distance <- lm(range.distance~sisters.split.time)
# record the slope and intercept of these models  
  overlap.slope<-lm.overlap$coefficients[[2]]
  overlap.intercept<-lm.overlap$coefficients[[1]]
  distance.slope<-lm.distance$coefficients[[2]]
  distance.intercept<-lm.distance$coefficients[[1]]
  asym.slope<-lm.asym$coefficients[[2]]
  asym.intercept<-lm.asym$coefficients[[1]]
  
# Calculate the mean range overlap between each sister species and the members of it's outgroup
  # TOmean = average of the range overlap between sister species - this mean OG overlap across all sister species pairs
  # TOSD = the SD of these values

  df_OG<-data.frame(model = NA, TO = NA)
  OG.rasters <- list()
  OG<-findOG(phy, sisters)
  for ( i in 1:length(OG)) {
    raster.nums.OG<-c()  
    for (j in 1:length(OG[[i]])) {
      raster.nums.OG[j] <- df[which(df[,2] == OG[[i]][j]), "rasterMatch"]
    }
    OG.rasters[[length(OG.rasters) + 1]] <- raster.nums.OG
  }  
  OG.overlap <- list()
  outgroup.overlap <- list()
  TO <- c()
  for (i in 1:length(OG.rasters)) {
    OG.overlap.sis2 <- c()
    OG.overlap.sis1 <- c()
    for (j in 1:length(OG.rasters[[i]])) {
      sis.1 <- rasters[[as.numeric(sisters.rasters[[i]][1])]] 
      sis.2 <- rasters[[as.numeric(sisters.rasters[[i]][2])]] 
      og <- rasters[[as.numeric(OG.rasters[[i]][j])]] 
      sis1.overlap<-sis.1+og
      sis2.overlap<-sis.2+og
      sis.overlap<-sis.1+sis.2
      if(any(sis1.overlap@data@values == 2, na.rm=T)) {
        if (sum(sis.1@data@values, na.rm=T) < sum(og@data@values, na.rm=T)) {
          range.overlap <- length(which(sis1.overlap@data@values == 2))/ sum(sis.1@data@values, na.rm=T)
        }
        if (sum(sis.1@data@values, na.rm=T) > sum(og@data@values, na.rm=T)) {
          range.overlap <- length(which(sis1.overlap@data@values == 2))/ sum(og@data@values, na.rm=T)
        }
      } else { range.overlap <- 0 }
      OG.overlap.sis1[j] <- range.overlap
      if(any(sis1.overlap@data@values == 2, na.rm=T)) {
        if (sum(sis.2@data@values, na.rm=T) <= sum(og@data@values, na.rm=T)) {
          range.overlap <- length(which(sis2.overlap@data@values == 2))/ sum(sis.2@data@values, na.rm=T)
        }
        if (sum(sis.2@data@values, na.rm=T) > sum(og@data@values, na.rm=T)) {
          range.overlap <- length(which(sis2.overlap@data@values == 2))/ sum(og@data@values, na.rm=T)
        }
      } else { range.overlap <- 0 }
      OG.overlap.sis2[j] <- range.overlap
      og.overlap.list<-list(OG.overlap.sis1, OG.overlap.sis2)
    }  
    OG.overlap[[i]]<-og.overlap.list
    outgroup.overlap[[i]]<-mean(mean(og.overlap.list[[1]], mean(og.overlap.list[[2]])))
    TO[i]<- sisters.overlap[[i]] - outgroup.overlap[[i]]
  }
  TO.total<- mean(TO, na.rm=T)
  TO.var <- sd(TO, na.rm=T)

# Species' range sizes
  species.range.sizes<-c()
  for (i in 1:length(rasters)) {
    species.range.sizes[i] <- length(which(!is.na(rasters[[i]]@data@values)))
  }
# standardise species range sizes
  species.range.sizes.stand <- species.range.sizes/max(species.range.sizes)
# range size mean, SD, and skew
  range.skew <- skewness(species.range.sizes.stand)
  range.av <- mean(species.range.sizes.stand, na.rm=T)
  range.sd <- sd(species.range.sizes.stand, na.rm=T)

# Phylogenetic imbalance and tempo metrics

  phy2<-as.treeshape(phy)
  beta <- maxlik.betasplit(phy2, up = 10, remove.outgroup = FALSE, confidence.interval = "none", conf.level = 0.95, size.bootstrap = 100)$max_lik
  collessI<-colless(phy2)
  sackinI<-sackin(phy2)
  ltt<-ltt(phy, plot=FALSE)
  gamma<-ltt$gamma
  
# sister species overlap metrics'
  
  symp0 <- length(which(sisters.overlap == 0))/length(sisters.overlap)
  symp50 <- length(which(sisters.overlap >= 0.5))/length(sisters.overlap)
  symp75 <- length(which(sisters.overlap >= 0.75))/length(sisters.overlap)
  symp90 <- length(which(sisters.overlap >= 0.9))/length(sisters.overlap)
  symp100 <- length(which(sisters.overlap == 1))/length(sisters.overlap)
  av.overlap <- mean(unlist(sisters.overlap), na.rm=T)
  asym <- mean(sister.asymmetry, na.rm=T)
  sisterpairs <- length(sisters)
  bimodality100<-((symp0*sisterpairs) *  (symp100*sisterpairs)) / ((sisterpairs/2)*(sisterpairs/2))
  bimodality90<-((symp0*sisterpairs) *  (symp90*sisterpairs)) / ((sisterpairs/2)*(sisterpairs/2))
  bimodality75<-((symp0*sisterpairs) *  (symp75*sisterpairs)) / ((sisterpairs/2)*(sisterpairs/2))
  bimodality50<-((symp0*sisterpairs) *  (symp50*sisterpairs)) / ((sisterpairs/2)*(sisterpairs/2))
  overlap<-unlist(sisters.overlap)
  kurt<- kurtosis(overlap)
  skew <- skewness(overlap)

  # this is the vector of summary statistics
  analysis <- data.frame(RO0 = symp0, RO50 = symp50, RO75 = symp75, RO90 = symp90, RO100 = symp100, ROmean = av.overlap,
                         asymmean = asym, bimod100 = bimodality100, bimod90 = bimodality90, bimod75 = bimodality75, bimod50 = bimodality50,
                         ROkurt = kurt, ROskew = skew, TOmean = TO.total, TOsd = TO.var,
                         RSskew= range.skew, RSmean = range.av, RSsd = range.sd, RDmean = range.distance.av, RDsd = range.distance.sd,
                         ROslope= overlap.slope, ROintercept = overlap.intercept, asymslope = asym.slope, asymintercept = asym.intercept, 
                         RDslope = distance.slope, RDintercept = distance.intercept,
                         Beta = beta, CI = collessI, SI = sackinI, gamma = gamma
                        )
# results return the vector of summary statistics as well as the measured values of sister species range overlap, range size asymmetry, TO, divergence times, and the species rasters
  Results <- list(sister.overlap = overlap,
                  sisters.split.time = sisters.split.time,
                  sister.rasters = sisters.rasters,
                  sister.asymmetry = sister.asymmetry,
                  Tip.Overlap = TO,
                  analysis = analysis)
return(Results)
}

# function to identify sister species pairs
findSisters <- function(phy,solve.polytomies=FALSE,include.polytomies=TRUE){
  require(phytools)
  n <- length(phy$tip.label)
  nb.node <- phy$Nnode
  if (n < 4)
    stop("tree has fewer than 4 tips")
  if(solve.polytomies){
    phy=multi2di(phy)
    include.polytomies=FALSE
  }
  if(include.polytomies){
    sister.ancestors <- which(tabulate(phy$edge[, 1][phy$edge[, 2] <= n])>=2)
    cherries <- vector("list",length(sister.ancestors))
    for(i in 1:length(cherries)){
      x <- getDescendants(phy, sister.ancestors[i])
      cherries[[i]] <- phy$tip.label[x]
      if(length(which(is.na(cherries[[i]])))>0) cherries[[i]]<-NULL  #remove internal nodes
    }
  }
  if(!include.polytomies){
    sister.ancestors <- which(tabulate(phy$edge[, 1][phy$edge[, 2] <= n])==2)
    cherries <- vector("list",length(sister.ancestors))
    for(i in 1:length(cherries)){
      x <- getDescendants(phy, sister.ancestors[i])
      cherries[[i]] <- c(phy$tip.label[x[1]], phy$tip.label[x[2]])
    }
  }
  cherries[!sapply(cherries,is.null)]
  return(cherries)
}

# function to find the species in the outgroup of a sister species pair
findOG <- function(phy,sistips){
  require(phytools)
  mrca <- mrca(phy)
  og <- list()
  for(i in 1:length(sistips)){
    sistip1 <- sistips[[i]][1]
    sistip2 <- sistips[[i]][2]
    curr.mrca <- mrca[which(rownames(mrca)==sistip1),which(colnames(mrca)==sistip2)]
    curr.anc <- phy$edge[,1][which(phy$edge[,2]==curr.mrca)]
    curr.tips <- phy$tip.label[getDescendants(phy,curr.anc)[getDescendants(phy,curr.anc)<Ntip(phy)+1]]
    curr.tips1 <- c(sistip1,sistip2)
    og[[i]] <- setdiff(curr.tips,curr.tips1)
  }
  return(og)
}

# this function creates a ENMTools clade object from the simulated data - used for calculating an Age-Range-Correlation with the ENMTools package
createENMToolsclade <- function(rasters, df, phy, env, drop.fossil=T) {
  if(drop.fossil==T){phy<-drop.fossil(phy)}
  root.age.true <- max(findMRCA(phy, type=c("height")))
  phy$edge.length<-phy$edge.length/root.age.true
  
  root.age<-1
  if(class(rasters)=="RasterStack") {rasters <- unstack(rasters)}
  df.dropped <- df[which(df[,2] %in% phy$tip.label),]
  df.dropped <- df.dropped[order(as.numeric(df.dropped[,10]), decreasing=FALSE),]
  rasters.dropped <- rasters[as.numeric(df.dropped[, "rasterMatch"])]
  species<-list()
  for (i in 1:length(phy$tip.label)) {
    a <- enmtools.species()
    a$species.name <- as.character(df.dropped[i, 2])
    a$range<-rasters.dropped[[i]]
    species[[i]] <- a 
  }
  test.clade <- enmtools.clade(species=species, tree=phy)
  return(test.clade)
}

