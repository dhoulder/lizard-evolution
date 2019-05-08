# functions for displaying information in runtime such as a map and summary stats

library(data.table)


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
    gene.flow.max.distance <- plotItems[["gene.flow.max.distance"]]

    # call genome.colours function to turn gene positions into R, G & B
    deme.colours <- genome.colour(demetable.species, genome.columns, gene.flow.max.distance)

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

    cat("\nCurrent time:", sp.summary$current.time, "\tSpecies:", sp.summary$speciesID, "of", sp.summary$speciesCount, "extant species",
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

genome.colour <- function(demetable, genome.columns, gene.flow.max.distance, genome.extremes=NA) {

  # this function uses the gene.pos columns, or genome.extremes if provided, to generate red, green, blue colours (from 0 to 1)
  #  for each deme based on genetic position.
  # the simplest approach is for 3 gene dismensions to translate to RGB
  # a more general approach uses an ordination

  # the no ordination method
  max.dist <- gene.flow.max.distance * 1.1  # this should be the genetic distance for maximum colour intensity in any one dimension
  genome.dimensions <- length(genome.columns)

  if (is.matrix(genome.extremes)) {
    use.extremes.matrix <- TRUE
  } else {
    use.extremes.matrix <- FALSE
  }
  
  # this is an inelegant method, but I can't effectively reference columns in a data.table via a variable
  genomes.species.df <- as.data.frame(demetable[, ..genome.columns])

  for (col in 1:genome.dimensions) {
    if (use.extremes.matrix) {
      mid.range <- mean(genome.extremes[, col])
      col.span <- genome.extremes[1, col] - genome.extremes[2, col]
    } else {
      col.range <- range(genomes.species.df[, col])
      mid.range <- mean(col.range)
      col.span  <- col.range[2] - col.range[1]
    }
    
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

genome.mean <- function(demetable, genome.columns, by.species = FALSE) {

  # this function calculates the mean genome position for the species, weighted by deme amount
  if (by.species) {
    sp_means <- demetable[, lapply(.SD, weighted.mean, amount), by=species_name, .SDcols=(genome.columns)]
  } else {
    sp_means <- demetable[, lapply(.SD, weighted.mean, amount), .SDcols=(genome.columns)]
  }

  return(sp_means)
}

get.genome.extremes <- function(demetable, genome.columns, extremes=NA) {
  
  # this function calculates the extreme values for each genome dimension
  genome.dimensions <- length(genome.columns)  
  
  # test if a valid matrix of extreme values has been provided
  update.extremes <- FALSE
  if (!is.na(extremes)) {
    if (is.matrix(extremes)) {
      if (dim(extremes) == c(2, genome.dimensions)) {update.extremes <- TRUE}
    }
  }

  if (!update.extremes) {
    extremes <- matrix(data=NA, nrow=2, ncol=genome.dimensions)
  }
  
  for (i in 1:genome.dimensions) {
    k <- genome.columns[i]
    extremes[1, i] <- max(extremes[1, i], max(demetable[, ..k]), na.rm = TRUE)
    extremes[2, i] <- min(extremes[2, i], min(demetable[, ..k]), na.rm = TRUE)
  }
  
  return(extremes)
}

display.to.file.start <- function(image_dir, time, image_filename = "plot", plot_rows=2, plot_cols=3) {
  image.width	<- 2000
  image.height	<- 1200

  output.filename <- paste(image_dir, image_filename, time, ".png", sep='')
  png(output.filename, width = image.width, height=image.height)
  
  if (plot_rows==2 & plot_cols==3) {
    par(mar=c(5,5,4,2), cex.main=1.3, cex.lab=2, cex.axis=1.5)
    layout(mat = matrix(c(1, 2, 3, 4, 5, 6), ncol = plot_cols, byrow = TRUE))   # Widths of the two columns)
  } else if (plot_rows==2 & plot_cols==2) {
    par(mar=c(2,2,4,1), cex.main=0.8, cex.lab=0.8)
    layout(mat = matrix(c(1, 3, 2, 4), ncol = plot_cols))   # Widths of the two columns)
  } else {
    par(mfcol=c(plot_rows, plot_cols))
  }
}

display.to.file.stop <- function() {
  dev.off()
}

display.to.screen.start <- function(window_name = "plot", plot_rows=2, plot_cols=2) {
  
  x11(width=10, height=7, title=window_name)
  par(mar=c(2,2,4,1), cex.lab=1)
  layout(mat = matrix(c(1, 3, 2, 4), ncol = plot_cols))   # Widths of the two columns)
}

display.to.screen.stop <- function() {
  dev.off()
}

display.update.multispecies <- function(plotItems, 
                                        env=NA, 
                                        allDemes=NA, 
																				includes.map=TRUE,
                                        plot_env=T, 
                                        plot_demes_amount_position=F, 
                                        plot_species_ranges=F, 
                                        plot_richness=F, 
                                        plot_genome_scatter=F, 
                                        plot_genome_map=F, 
                                        plot_tree=F, 
                                        plot_mainheader=F,
                                        plot_species_over_time=F) {
  # plotItems is a list of named components which can be included in the display

  dot.size.scaler   <- 0.9  # 1 is good for a 100 x 100 plot (4 x4), smaller for higher resolution
  pane.header.size  <- 2

  # to maintain a consistent scale of environment colours, fix two pixels to the extremes of the range
  # currently this is the initial range of values (0 to 100) plus the amplitude of cyclical variation
  # 10, so the range should be -10 to 110
  #plotItems[["env"]][1] <- -10
  #plotItems[["env"]][2] <- 110
  
  if (plot_mainheader) {
		main.header <- ""
		
		if (length(plotItems[["model.params"]]) > 0) {
			main.header <- paste(main.header, "\n\nNiche breadth (1):", plotItems[["model.params"]][[1]],
													 "\nNiche evol rate:", plotItems[["model.params"]][[2]],
													 "\nDispersal:", plotItems[["model.params"]][[3]],
							 "\nSpeciation distance:", plotItems[["model.params"]][[4]])
		}		
		
		if (length(plotItems[["current.time"]]) > 0) {
			main.header <- paste(main.header, "\n\n\nCurrent time:", plotItems[["current.time"]])
		}
		
		if (length(plotItems[["diversity_time"]]) > 0) {
			diversity_time.df <- plotItems[["diversity_time"]]
			last_row <- nrow(diversity_time.df)
			main.header <- paste(main.header,
														"\n\nCurrent species:\t", diversity_time.df[last_row, "current_species"],
														"\nExtinct species: \t", diversity_time.df[last_row, "extinct_species"])
		}
 
    plot.new()
    text(0, 0.6, main.header, cex=2, font=2, adj= c(0, 0.5))
  }
   
  # convert row / column coordinates to continuous coordinates, if they don't match the environment raster
  # for example, if the raster is in degrees, then row and column numbers won't plot correctly

  if(includes.map) {
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
		
		# check if there is a shapefile to plot on the map
		if (length(plotItems[["shape"]]) > 0) {
			shape <- plotItems[["shape"]]
			have.shape <- TRUE
		} else {
			have.shape <- FALSE
		}

		#these.symbols  <- getDiscreteSymbols(allDemes$species_name)  
		# this WAS used to plot a different symbol per species, but they were hard to distinguish and caused plotting errors - perhaps due to going beyond the range of symbol numbers....
  }

  if (plot_demes_amount_position) {
    # assign colours to niche0.position, based on the 250 colours defined above in display.initialise.2by2()
    colour.count   <- 250
    colour.indices <- round((allDemes$niche_centre_0 - min.env) * colour.count / range.env)
    these.colours  <- my.colours[colour.indices]
    these.sizes    <- sqrt(allDemes$amount) * 2 * dot.size.scaler
    
    plot(env, main="Demes with environment", cex.main=pane.header.size, cex.axis=1.3) # a first plot to ensure the correct legend
    plot(env, col="white", add=T, legend=F)  
		points(allDemes$col, allDemes$row, col=these.colours, pch=19, cex=these.sizes)
		if (have.shape) {
			plot(shape, add=T)
		}
  }
  
  if (plot_species_ranges) {
    these.sizes    <- sqrt(allDemes$amount) * 2 * dot.size.scaler
      
		plot(env, main="species ranges", col="white", legend=FALSE, cex.main=pane.header.size)
		points(allDemes$col, allDemes$row, col=allDemes$species, pch=19, cex=these.sizes)
		if (have.shape) {
			plot(shape, add=T)
		}	
  }
  
  if (plot_tree) {
    plot(plotItems[["tree"]], type="phylogram", main="Phylogeny", cex.main=pane.header.size, root.edge = TRUE)
    add.scale.bar(col="blue", lwd=1.25)
  }
  
  if(plot_species_over_time) {
    diversity_time.df <- plotItems[["diversity_time"]]
		xrange <- range(diversity_time.df$time)
    plot(diversity_time.df$time, log2(diversity_time.df$current_species), 
			xlim = xrange[2:1], ylim=c(0, max(0, log2(max(diversity_time.df$current_species)))), 
			type="b", 
			pch=19, 
			col="blue", 
			xlab= "Time", 
			ylab="Current species log2", 
			main="Species through time (log scale)", cex.main=pane.header.size)
  }
  
  if (plot_genome_scatter | plot_genome_map) {
    	  
		if (length(plotItems[["genome.extremes"]]) > 0) {
			genome.extremes <- plotItems[["genome.extremes"]]
		}
    
    genome.columns <- plotItems[["genome.columns"]]
    gene.flow.max.distance <- plotItems[["gene.flow.max.distance"]]
    
    # call genome.colours function to turn gene positions into R, G & B
    deme.colours <- genome.colour(allDemes, genome.columns, gene.flow.max.distance, genome.extremes)
    these.colours <- rgb(red = deme.colours[,1], green = deme.colours[,2], blue = deme.colours[,3])
    
    if (plot_genome_map) {
      plot(env, main="Genomic divergence - map", col="white", legend=FALSE, cex.main=pane.header.size)   # a blank environment map to highlight gene colours
			points(allDemes$col, allDemes$row, col=these.colours, pch=19, cex=dot.size.scaler)
			if (have.shape) {
				plot(shape, add=T)
			}
    }
    
    if (plot_genome_scatter) {
      these.sizes   <- sqrt(allDemes$amount) * 1.5

			if (is.matrix(genome.extremes)) {
				use.extremes.matrix <- TRUE
			} else {
				use.extremes.matrix <- FALSE
			}
				
			# give the plots a standard extent, to see the dispersion increasing.  But allow the extent to increase when needed
			plot.limit  <- gene.flow.max.distance * 0.8
				
			if (use.extremes.matrix) {
				plot.limits <- c(min(genome.extremes[2, 1], (plot.limit * -1)), max(genome.extremes[1, 1], plot.limit), min(genome.extremes[2, 2], (plot.limit * -1)), max(genome.extremes[1, 2], plot.limit))
			} else {
				plot.limits <- as.numeric(allDemes[, .(min(genetic_position_0, (plot.limit * -1)), max(genetic_position_0, plot.limit), min(genetic_position_1, (plot.limit * -1)), max(genetic_position_1, plot.limit))])
			}
		
      genome_scatter_x <- allDemes$genetic_position_0
      genome_scatter_y <- allDemes$genetic_position_1
			
			# subsample the points for the genome scatterplot, where the numbers are very large
			demecount <- length(genome_scatter_x)
			max_points <- 150000
			if (demecount > max_points) {
				sample_indices <- sample(1:demecount, max_points)
				genome_scatter_x <- genome_scatter_x[sample_indices]
				genome_scatter_y <- genome_scatter_y[sample_indices]
				these.colours <- these.colours[sample_indices]
				these.sizes <- these.sizes[sample_indices]
				#these.symbols <- these.symbols[sample_indices]
			}
			
			plot(genome_scatter_x, genome_scatter_y, col=these.colours, cex=these.sizes, main="Genomic divergence - scatter",
           pch=19, xlab="Genome axis 1", ylab="Genome axis 2", xlim=plot.limits[1:2], ylim=plot.limits[3:4], cex.main=pane.header.size)
     
      # add a weighted genome mean for each species, to the plot
      means <- genome.mean(allDemes, genome.columns, by.species = TRUE)
      points(means$genetic_position_0, means$genetic_position_1, pch=3, cex=1.8)
    }  
  }
  
  if (plot_richness) {

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
    plot(richness.ras, main=main.txt, cex.main=pane.header.size)
		if (have.shape) {
			plot(shape, add=T)
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