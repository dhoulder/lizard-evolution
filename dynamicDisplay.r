# functions for displaying information in runtime such as a map and summary stats

display.initialise <- function() {
  if (!output_to_file) {
    x11(15,11)
  }

  stored_par <- par(mfcol=c(1,1))

  my.colours.function <- colorRampPalette(colors = c("blue", "yellow",  "red"))
  my.colours    <- my.colours.function(250)
  return(list(my.colours, stored_par))
}

display.initialise.double <- function() {

  if (!output_to_file) {
    x11(18,9) # open on-screen display
  }

  stored_par <- par(mfcol=c(1,2))

  my.colours.function <- colorRampPalette(colors = c("blue", "yellow",  "red"))
  my.colours    <- my.colours.function(250)

  my.diffcolours.function <- colorRampPalette(colors = c("red", "white", "blue"))
  my.diffcolours    <- my.diffcolours.function(250)

  return(list(my.colours, my.diffcolours, stored_par))
}

display.initialise.2by2 <- function(output_to_file = F) {

  if (!output_to_file) {
    x11(10,10) # open on-screen display
  }

  stored_par <- par(mfcol=c(2,2))

  my.colours.function <- colorRampPalette(colors = c("blue", "yellow",  "red"))
  my.colours    <- my.colours.function(250)

  my.diffcolours.function <- colorRampPalette(colors = c("red", "white", "blue"))
  my.diffcolours    <- my.diffcolours.function(250)

  return(list(my.colours, my.diffcolours, stored_par))
}

display.update <- function(plotItems) {
  # elements is a list of named components to include in the display

  # this function relies on the data to be plotted, being in scope, rather than passed as argument
  # this can be revised if it is a problem

  dot.size.scaler <- 0.8  # 1 is good for a 100 x 100 plot (4 x4), smaller for higher resolution

  if (length(plotItems[["env"]]) > 0) {

    if (length(plotItems[["current.time"]]) > 0) {
      main.header <- paste("Time:", plotItems[["current.time"]])
    } else {
      main.header <- ""
    }

    if (length(plotItems[["niche.params"]]) > 0) {
      main.header <- paste(main.header, "\nNiche breadth:", plotItems[["niche.params"]][[1]],
                           "\tNiche evol rate:", plotItems[["niche.params"]][[2]],
                           "\tDispersal:", plotItems[["niche.params"]][[3]])
    }

    if (length(plotItems[["demes_genecolour"]]) > 0) {
      plot(plotItems[["env"]], main=main.header, col="white")   # trying a blank environment map to highlight gene colours
    } else {
      plot(plotItems[["env"]], main=main.header, col=my.colours)
    }

  }

  if (length(plotItems[["demes_amount_position"]]) > 0) {
    demetable.species <- plotItems[["demes_amount_position"]]
    these.colours <- my.colours[round(demetable.species$niche1.position*10)]
    these.sizes   <- sqrt(demetable.species$amount) * 2 * dot.size.scaler

    #points(demetable.species$col, demetable.species$row, bg=these.colours, pch=21, fg="black", cex=these.sizes)
    points(demetable.species$col, (environment.dimension-demetable.species$row), col=these.colours, pch=19, cex=these.sizes) # trying these settings for html animation
  }

  if (length(plotItems[["demes_amount_position_diff"]]) > 0) {
    demetable.species <- plotItems[["demes_amount_position_diff"]]

    env <- plotItems[["env"]]
    set(demetable.species, j="env", value=env[demetable.species$cellID])
    demetable.species[, diff:=niche1.position-env]
    breadth.max <- demetable.species[, max(niche1.breadth)]

    colour.nums <- round(demetable.species$diff * (248 / breadth.max)) + 124
    colour.nums[colour.nums < 1] <- 1
    colour.nums[colour.nums >250] <- 250
    these.colours <- my.coloursdiff[colour.nums]
    #points(demetable.species$x, demetable.species$y, bg=these.colours, pch=21, fg="black", cex=1)
    points(demetable.species$col, (environment.dimension-demetable.species$row), col=these.colours, pch=19, cex=dot.size.scaler)
  }

  if (length(plotItems[["demes_genecolour"]]) > 0) {

    demetable.species <- plotItems[["demes_genecolour"]]
    env <- plotItems[["env"]]
    genome.columns <- plotItems[["genome.columns"]]

    # call genome.colours function to turn gene positions into R, G & B
    deme.colours <- genome.colour(demetable.species, genome.columns)

    these.colours <- rgb(red = deme.colours[,1], green = deme.colours[,2], blue = deme.colours[,3])
    points(demetable.species$col, (environment.dimension-demetable.species$row), col=these.colours, pch=19, cex=dot.size.scaler)

    ########################################################################################
    these.sizes   <- sqrt(demetable.species$amount) * 1.5

    # give the plots a standard extent, to see the dispersion increasing.  But allow the extent to increase when needed
    plot.limit  <- speciation.gene.distance * 0.6
    plot.limits <- as.numeric(demetable.species[, .(min(gene.pos1, (plot.limit * -1)), max(gene.pos1, plot.limit), min(gene.pos2, (plot.limit * -1)), max(gene.pos2, plot.limit))])

    plot(demetable.species$gene.pos1, demetable.species$gene.pos2, col=these.colours, cex=these.sizes,
         pch=20, xlab="Genome axis 1", ylab="Genome axis 2", xlim=plot.limits[1:2], ylim=plot.limits[3:4])

    # add a weighted genome mean for the species, to the plot
    means <- genome.mean(demetable.species, genome.columns)
    points(means[1], means[2], pch=3, cex=1.5)
    ########################################################################################

  }

  if (length(plotItems[["one_deme"]] > 0)) {
    deme <- plotItems[["one_deme"]]
    points(deme$col, deme$row, col="black", pch=0)
  }

  return()
}

text.update <- function(textItems) {

  if (length(textItems[["species_range_niche"]]) > 0) {

    sp.summary <- textItems[["species_range_niche"]]
    geneflow.prob <- prob.geneflow(sp.summary$gen.distance.max, zero_flow_dist=speciation.gene.distance)

    cat("\nTime:", sp.summary$current.time, "\tSpecies:", sp.summary$current.speciesID,
        "\nRange:", sp.summary$range,
        "\tTotal amount:", sp.summary$total_amount,
        "\n\n\tNiche",
        "\nNiche1 position\tMean:", sp.summary$niche1.position.mean, "\tSD:", sp.summary$niche1.position.sd,
        "\nNiche1 breadth\tMean:", sp.summary$niche1.breadth.mean, "\tSD:", sp.summary$niche1.breadth.sd,
        "\nNiche1 species\tMin: ", sp.summary$niche1.sp.min, "\tMax:", sp.summary$niche1.sp.max,
        "\nTotal niche breadth:", (sp.summary$niche1.sp.max - sp.summary$niche1.sp.min),
        "\n\n\tGenetic divergence",
        "\nMaximum:", sp.summary$gen.distance.max,
        "\tMedian:", sp.summary$gen.distance.median,
        "\ngeneflow prob at max distance:", paste(round(geneflow.prob * 100, 2), "%", sep=""),
        "\nThis step:", sp.summary$elapsed.time.step, "\tTotal time elapsed:", sp.summary$elapsed.time.total,
        "\n*****************************************\n")
  }

  return()
}

genome.colour <- function(demetable.species, genome.columns) {

  # this function uses the gene.pos columns to generate red, green, blue colours (from 0 to 1)
  #  for each deme based on genetic position.
  # the simplest approach is for 3 gene dismensions to translate to RGB
  # a more general approach uses an ordination

  # the no ordination method
  max.dist <- 25  # this should be the genetic distance for maximum colour intensity in any one dimension
  genome.dimensions <- length(genome.columns)

  # this is an inelegant method, but I can't effectively reference columns in a data.table via a variable
  genomes.species.df <- as.data.frame(demetable.species[, ..genome.columns])

  for (col in 1:genome.dimensions) {
    col.range <- range(genomes.species.df[, col])
    mid.range <- mean(col.range)
    col.span  <- col.range[2] - col.range[1]

    # if the range of values in a dimesion is greater than max.dist, expand max.dist so that all values
    # fit within the range of colours intensities
    max.dist <- (max(max.dist, col.span))

    max.range <- mid.range + (max.dist/2)
    min.range <- mid.range - (max.dist/2)
    new.col.name <- paste("genome.colour", col, sep = '')
    genomes.species.df[, new.col.name] <- (genomes.species.df[, col] - min.range) / max.dist
    genomes.species.df[which(genomes.species.df[new.col.name] < 0), new.col.name] <- 0
    genomes.species.df[which(genomes.species.df[new.col.name] > 1), new.col.name] <- 1
  }

  return(genomes.species.df[, (genome.dimensions + 1):(genome.dimensions * 2)])  # return a data frame of just the colour columns

}

genome.mean <- function(demetable.species, genome.columns) {

  # this function calculates the mean genome position for the species, weighted by deme amount

  genome.dimensions <- length(genome.columns)
  genome.mean.pos   <- genome.columns # a lazy way to get a vector of the right length

  for (col in 1:genome.dimensions) {
    genome.col <- as.data.frame(demetable.species[, genome.columns[col], with=F])
    genome.mean.pos[col] <- weighted.mean(genome.col[,1], demetable.species$amount)
  }

  return(genome.mean.pos)
}

display.to.file.start <- function(image_dir, time, image_filename = "plot") {
  image.width = 960
  image.width = 960

  output.filename <- paste(image_dir, image_filename, time, ".png", sep='')
  png(output.filename, width = 960, height=960)
  #par(my.par)
  par(mfcol=c(2,2))
}

display.to.file.stop <- function() {
  dev.off()
}

