# functions for displaying information in runtime such as a map and summary stats

library(data.table)

display.initialise <- function() {
  if (!image_to_file) {
    x11(15,11)
  }

  stored_par <- par(mfcol=c(1,1))

  my.colours.function <- colorRampPalette(colors = c("blue", "yellow",  "red"))
  my.colours    <- my.colours.function(250)
  return(list(my.colours, stored_par))
}

display.initialise.double <- function() {

  if (!image_to_file) {
    x11(18,9) # open on-screen display
  }

  stored_par <- par(mfcol=c(1,2))

  my.colours.function <- colorRampPalette(colors = c("blue", "yellow",  "red"))
  my.colours    <- my.colours.function(250)

  my.diffcolours.function <- colorRampPalette(colors = c("red", "white", "blue"))
  my.diffcolours    <- my.diffcolours.function(250)

  return(list(my.colours, my.diffcolours, stored_par))
}

display.initialise.2by2 <- function(image_to_file = F) {

  if (!image_to_file) {
    x11(10,10) # open on-screen display
  }

  stored_par <- par(mfcol=c(2,2))

  my.colours.function <- colorRampPalette(colors = c("blue", "yellow",  "red"))
  my.colours    <- my.colours.function(250)

  my.diffcolours.function <- colorRampPalette(colors = c("red", "white", "blue"))
  my.diffcolours    <- my.diffcolours.function(250)

  return(list(my.colours, my.diffcolours, stored_par))
}

display.initialise.colours <- function() {
  my.colours.function <- colorRampPalette(colors = c("blue", "yellow",  "red"))
  my.colours    <- my.colours.function(250)
  
  my.diffcolours.function <- colorRampPalette(colors = c("red", "white", "blue"))
  my.diffcolours    <- my.diffcolours.function(250)
  
  return(list(my.colours, my.diffcolours))
}

display.update <- function(plotItems, plot_env=T, plot_genome_scatter=T, plot_genome_map=T) {
  # plotItems is a list of named components to include in the display

  dot.size.scaler <- 0.8  # 1 is good for a 100 x 100 plot (4 x4), smaller for higher resolution

  if (length(plotItems[["env"]]) > 0) {
    
    # to maintain a consistent scale of environment colours, fix two pixels to the extremes of the range
    # currently this is the initial range of values (0 to 100) plus the amplitude of cyclical variation
    # 10, so the range should be -10 to 110
    plotItems[["env"]][1] <- -10
    plotItems[["env"]][2] <- 110

    if (length(plotItems[["current.time"]]) > 0) {
      main.header <- paste("Time:", plotItems[["current.time"]])
    } else {
      main.header <- ""
    }

    if (length(plotItems[["model.params"]]) > 0) {
      main.header <- paste(main.header, "\nNiche breadth:", plotItems[["model.params"]][[1]],
                           "\tNiche evol rate:", plotItems[["model.params"]][[2]],
                           "\tDispersal:", plotItems[["model.params"]][[3]])
    }

    if ( (length(plotItems[["demes_genecolour"]]) > 0) | (!plot_env) ) {
      # plot the environment in colours to get the legend, then overplot in white
      #plot(plotItems[["env"]], main=main.header, col=my.colours)
      plot(plotItems[["env"]], main=main.header, col="white")   # trying a blank environment map to highlight gene colours
    } else {
      plot(plotItems[["env"]], main=main.header, col=my.colours)
    }
  }
	
	# convert row / column coordinates to continuous coordinates, if they don't match the environment raster
	# for example, if the raster is in degrees, then row and column numbers won't plot correctly
	if (length(plotItems[["demes_amount_position"]]) > 0 | 
			length(plotItems[["demes_amount_position_diff"]]) > 0 |
			length(plotItems[["demes_genecolour"]]) > 0){
		
		# check if the environment layer has coordinates matching the row and column numbers
		env_temp <- plotItems[["env"]]
		env.has.rowcol.coords <- (env_temp@ncols == env_temp@extent@xmax)
		
		# check the range of environment raster values
		min.env 	<- min(env_temp[], na.rm=T)
		range.env <- max(env_temp[], na.rm=T) - min.env
	}

  if (length(plotItems[["demes_amount_position"]]) > 0) {
    demetable.species <- plotItems[["demes_amount_position"]]
		
		# assign colours to niche0.position, based on the 250 colours defined above in display.initialise.2by2()
		colour.count   <- 250
    colour.indices <- round((demetable.species$niche0.position - min.env) * colour.count / range.env)
		these.colours  <- my.colours[colour.indices]
    these.sizes    <- sqrt(demetable.species$amount) * 2 * dot.size.scaler

    # replace the row and column values with x, y if needed
		if (env.has.rowcol.coords) {
			demetable.species$row <- environment.rows - demetable.species$row  # where row numbers are used for the y value, this converts
																																							# to standard y values where y=0 as at the bottom, not top
		} else  {
		demetable.species$col <- xFromCol(env_temp, col=demetable.species$col)		# replace row and column with x and y values
		demetable.species$row <- yFromRow(env_temp, row=demetable.species$row)
		}

    points(demetable.species$col, demetable.species$row, col=these.colours, pch=19, cex=these.sizes)
		
  }

  if (length(plotItems[["demes_amount_position_diff"]]) > 0) {
    demetable.species <- plotItems[["demes_amount_position_diff"]]

    #env.ras <- plotItems[["env"]]
    #demetable.species[, env:= env.ras[cellFromRowCol(env.ras, demetable.species$row, demetable.species$col)]]
    demetable.species[, diff:=niche_centre_0 - env_0]
    breadth.max <- demetable.species[, max(niche_breadth_0)]

    colour.nums <- round(demetable.species$diff * (248 / breadth.max)) + 124
    colour.nums[colour.nums < 1] <- 1
    colour.nums[colour.nums >250] <- 250
    these.colours <- my.coloursdiff[colour.nums]
		
    # replace the row and column values with x, y if needed
		if (env.has.rowcol.coords) {
			demetable.species$row <- environment.rows - demetable.species$row  # where row numbers are used for the y value, this converts
																																							# to standard y values where y=0 as at the bottom, not top
		} else  {
		demetable.species$col <- xFromCol(env_temp, col=demetable.species$col)		# replace row and column with x and y values
		demetable.species$row <- yFromRow(env_temp, row=demetable.species$row)
		}
		
    #points(demetable.species$x, demetable.species$y, bg=these.colours, pch=21, fg="black", cex=1)
    points(demetable.species$col, demetable.species$row, col=these.colours, pch=19, cex=dot.size.scaler)
  }

  if (length(plotItems[["demes_genecolour"]]) > 0) {

    demetable.species <- plotItems[["demes_genecolour"]]
    env <- plotItems[["env"]]
    genome.columns <- plotItems[["genome.columns"]]

    # call genome.colours function to turn gene positions into R, G & B
    deme.colours <- genome.colour(demetable.species, genome.columns)

    these.colours <- rgb(red = deme.colours[,1], green = deme.colours[,2], blue = deme.colours[,3])
		
    if (plot_genome_map) {
      # replace the row and column values with x, y if needed
  		if (env.has.rowcol.coords) {
  			demetable.species$row <- environment.rows - demetable.species$row  # where row numbers are used for the y value, this converts
  																																							# to standard y values where y=0 as at the bottom, not top
  		} else  {
  		demetable.species$col <- xFromCol(env_temp, col=demetable.species$col)		# replace row and column with x and y values
  		demetable.species$row <- yFromRow(env_temp, row=demetable.species$row)
  		}
  		
      points(demetable.species$col, demetable.species$row, col=these.colours, pch=19, cex=dot.size.scaler)
    }
    ########################################################################################
    if (plot_genome_scatter) {
      these.sizes   <- sqrt(demetable.species$amount) * 1.5

      # give the plots a standard extent, to see the dispersion increasing.  But allow the extent to increase when needed
      plot.limit  <- gene.flow.max.distance * 0.6
      plot.limits <- as.numeric(demetable.species[, .(min(genetic_position_0, (plot.limit * -1)), max(genetic_position_0, plot.limit), min(genetic_position_1, (plot.limit * -1)), max(genetic_position_1, plot.limit))])
  
      genome_scatter_x <- demetable.species[, genetic_position_0]
      genome_scatter_y <- demetable.species[, genetic_position_1]
      plot(genome_scatter_x, genome_scatter_y, col=these.colours, cex=these.sizes,
           pch=20, xlab="Genome axis 1", ylab="Genome axis 2", xlim=plot.limits[1:2], ylim=plot.limits[3:4])
  
      # add a weighted genome mean for the species, to the plot
      means <- genome.mean(demetable.species, genome.columns)
      points(means[1], means[2], pch=3, cex=1.5)
    }
    ########################################################################################

  }

  if (length(plotItems[["one_deme"]] > 0)) {
    deme <- plotItems[["one_deme"]]
    points(deme$col, deme$row, col="black", pch=0)
  }

  return(1)
}

text.update <- function(textItems) {

  if (length(textItems[["species_range_niche"]]) > 0) {

    sp.summary <- textItems[["species_range_niche"]]

    gen.distance.available <- (length(sp.summary$gen.distance.max) > 0)

    if (gen.distance.available) {
      geneflow.prob <- prob.geneflow(sp.summary$gen.distance.max, zero_flow_dist=gene.flow.max.distance)
    }

    cat("\nTime:", sp.summary$current.time, "\tSpecies:", sp.summary$speciesID, "of", sp.summary$speciesCount, "extant species",
        "\nRange:", sp.summary$range,
        "\tTotal amount:", sp.summary$total_amount,
        "\n\n\tNiche",
        "\nNiche0 position\tMean:", sp.summary$niche0.position.mean, "\tSD:", sp.summary$niche0.position.sd,
        "\nNiche0 breadth\tMean:", sp.summary$niche0.breadth.mean, "\tSD:", sp.summary$niche0.breadth.sd,
        "\nNiche0 species\tMin: ", sp.summary$niche0.sp.min, "\tMax:", sp.summary$niche0.sp.max,
        "\nTotal niche breadth:", (sp.summary$niche0.sp.max - sp.summary$niche0.sp.min))

    if (gen.distance.available) {
      cat("\n\n\tGenetic divergence",
        "\nMaximum:", sp.summary$gen.distance.max,
        "\tMedian:", sp.summary$gen.distance.median,
        "\ngeneflow prob at max distance:", paste(round(geneflow.prob * 100, 2), "%", sep=""))
    }

    cat("\n\nThis step:", sp.summary$elapsed.time.step, "\tTotal time elapsed:", sp.summary$elapsed.time.total,
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

display.to.file.start <- function(image_dir, time, image_filename = "plot", plot_rows=2, plot_cols=2) {
  image.width	<- 1600
  image.height	<- 1200

  output.filename <- paste(image_dir, image_filename, time, ".png", sep='')
  png(output.filename, width = image.width, height=image.height)

  par(mfcol=c(plot_rows, plot_cols))
}

display.to.file.stop <- function() {
  dev.off()
}

display.update.multispecies <- function(plotItems, env, allDemes, plot_env=T, plot_demes_amount_position=F, plot_species_ranges=F, plot_richness=F, plot_genome_scatter=F, plot_genome_map=F) {
  # plotItems is a list of named components to include in the display

  dot.size.scaler <- 0.9  # 1 is good for a 100 x 100 plot (4 x4), smaller for higher resolution

  # to maintain a consistent scale of environment colours, fix two pixels to the extremes of the range
  # currently this is the initial range of values (0 to 100) plus the amplitude of cyclical variation
  # 10, so the range should be -10 to 110
  #plotItems[["env"]][1] <- -10
  #plotItems[["env"]][2] <- 110
  
  if (length(plotItems[["current.time"]]) > 0) {
    main.header <- paste("Time:", plotItems[["current.time"]])
  } else {
    main.header <- ""
  }
  
  if (length(plotItems[["model.params"]]) > 0) {
    main.header <- paste(main.header, "\nNiche breadth:", plotItems[["model.params"]][[1]],
                         "\tNiche evol rate:", plotItems[["model.params"]][[2]],
                         "\tDispersal:", plotItems[["model.params"]][[3]])
  }
  # convert row / column coordinates to continuous coordinates, if they don't match the environment raster
  # for example, if the raster is in degrees, then row and column numbers won't plot correctly

  # check if the environment layer has coordinates matching the row and column numbers
  env.has.rowcol.coords <- (env@ncols == env@extent@xmax)
  
  # check the range of environment raster values
  min.env 	<- min(env[], na.rm=T)
  range.env <- max(env[], na.rm=T) - min.env

  # replace the row and column values with x, y if needed
  if (env.has.rowcol.coords) {
    allDemes$row <- environment.rows - allDemes$row  # where row numbers are used for the y value, this converts
    # to standard y values where y=0 as at the bottom, not top
  } else  {
    allDemes$column <- xFromCol(env, col=allDemes$col)		# replace row and column with x and y values
    allDemes$row <- yFromRow(env, row=allDemes$row)
  }

  these.symbols  <- getDiscreteSymbols(allDemes$species_name)

  if (plot_demes_amount_position) {

    # assign colours to niche0.position, based on the 250 colours defined above in display.initialise.2by2()
    colour.count   <- 250
    colour.indices <- round((allDemes$niche_centre_0 - min.env) * colour.count / range.env)
    these.colours  <- my.colours[colour.indices]
    these.sizes    <- sqrt(allDemes$amount) * 2 * dot.size.scaler
    
    plot(env, main=main.header, col="white")  
    points(allDemes$col, allDemes$row, col=these.colours, pch=these.symbols, cex=these.sizes)
  }
  
  if (plot_species_ranges) {
    these.sizes    <- sqrt(allDemes$amount) * 2 * dot.size.scaler
      
	plot(env, main=main.header, col="white", legend=FALSE)
	points(allDemes$col, allDemes$row, col=allDemes$species, pch=19, cex=these.sizes)
  }
  
  if (plot_richness) {

    #richness.dt <- allDemes[, .(cellRichness=length(species_name)), by=.(row, column)]
    richness.dt <- allDemes[, .(spRichness=.N), by=.(column, row)]
    
    # create a richness raster of the same size as env
    env.extent   <- extent(env)
    richness.ras <- env
    richness.ras[] <- 0
    if (env.has.rowcol.coords) {
      richness.dt[, cellIndex:=cellFromRowCol(richness.ras,(environment.rows - row), column)]
    } else {
      richness.dt[, cellIndex:=cellFromXY(richness.ras, cbind(column, row))]
    }

    richness.ras[richness.dt$cellIndex] <- richness.dt$spRichness
  
	main.txt <- paste("Species richness\t\tTotal species:", plotItems[["speciesCount"]])
    plot(richness.ras, main=main.txt)
  }

  if (plot_genome_scatter | plot_genome_map) {
    plot(env, main=main.header, col="white", legend=FALSE)   # a blank environment map to highlight gene colours
    
    genome.columns <- plotItems[["genome.columns"]]
    
    # call genome.colours function to turn gene positions into R, G & B
    deme.colours <- genome.colour(allDemes, genome.columns)
    
    these.colours <- rgb(red = deme.colours[,1], green = deme.colours[,2], blue = deme.colours[,3])
    
    if (plot_genome_map) {
      points(allDemes$col, allDemes$row, col=these.colours, pch=these.symbols, cex=dot.size.scaler)
    }

    if (plot_genome_scatter) {
      these.sizes   <- sqrt(allDemes$amount) * 1.5
      
      # give the plots a standard extent, to see the dispersion increasing.  But allow the extent to increase when needed
      plot.limit  <- gene.flow.max.distance * 0.8
      plot.limits <- as.numeric(allDemes[, .(min(genetic_position_0, (plot.limit * -1)), max(genetic_position_0, plot.limit), min(genetic_position_1, (plot.limit * -1)), max(genetic_position_1, plot.limit))])
      
      genome_scatter_x <- allDemes$genetic_position_0
      genome_scatter_y <- allDemes$genetic_position_1
      plot(genome_scatter_x, genome_scatter_y, col=these.colours, cex=these.sizes,
           pch=these.symbols, xlab="Genome axis 1", ylab="Genome axis 2", xlim=plot.limits[1:2], ylim=plot.limits[3:4])
      
      # add a weighted genome mean for each species, to the plot
      means <- genome.mean(allDemes, genome.columns)
      points(means[1], means[2], pch=3, cex=1.5)
    }  
  }

  return(1)
}

getDiscreteSymbols <- function(pointValues) {
  # pointValues should be a vector of class integer or factor, consecutively from 1
  
  if (! (class(pointValues) %in% c("factor", "integer"))) {exit}
  
  if (class(pointValues) == "factor") {
    pointValues <- as.integer(pointValues)
  }
  
  symbols <- 20 - pointValues
  return(symbols)
}

add.alpha <- function(col, alpha=1){
  # add an alpha value to a colour
  if(missing(col)) {stop("Please provide a vector of colours.")}
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}