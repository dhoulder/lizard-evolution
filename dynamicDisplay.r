# functions for displaying information in runtime such as a map and summary stats

display.initialise <- function() {
  windows(15,11)

  my.colours.function <- colorRampPalette(colors = c("blue", "yellow",  "red"))
  my.colours    <- my.colours.function(250)
  return(my.colours)
}

display.initialise.double <- function() {
  windows(18,9)
  par(mfcol=c(1,2))

  my.colours.function <- colorRampPalette(colors = c("blue", "yellow",  "red"))
  my.colours    <- my.colours.function(250)

  my.diffcolours.function <- colorRampPalette(colors = c("red", "white", "blue"))
  my.diffcolours    <- my.diffcolours.function(250)

  return(list(my.colours, my.diffcolours))
}

display.update <- function(plotItems) {
  # elements is a list of named components to include in the display

  # this function relies on the data to be plotted being in scope, rather than passed as argument
  # this can be revised if it is a problem

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

    plot(plotItems[["env"]], main=main.header, col=my.colours)
    #plot(plotItems[["env"]], main=main.header, col="white")   # trying a blank environment map to highlight diff values
  }

  if (length(plotItems[["demes_amount_position"]]) > 0) {
    demetable.species <- plotItems[["demes_amount_position"]]
    these.colours <- my.colours[round(demetable.species$niche1.position*10)]
    these.sizes   <- sqrt(demetable.species$amount) * 2
    points(demetable.species$x, demetable.species$y, bg=these.colours, pch=21, fg="black", cex=these.sizes)
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
    #these.sizes   <- sqrt(demetable.species$amount) * 2
    points(demetable.species$x, demetable.species$y, bg=these.colours, pch=21, fg="black", cex=1)#these.sizes)
  }

  if (length(plotItems[["one_deme"]] > 0)) {
    deme <- plotItems[["one_deme"]]
    points(deme$x, deme$y, col="black", pch=0)
  }

  return()
}

text.update <- function(textItems) {

  if (length(textItems[["species_range_niche"]]) > 0) {

    sp.summary <- textItems[["species_range_niche"]]

    cat("\nTime:", sp.summary$current.time, "\tSpecies:", sp.summary$current.speciesID, "\tRange:", sp.summary$range,
        "\nNiche1 position\tMean:", sp.summary$niche1.position.mean, "\tSD:", sp.summary$niche1.position.sd,
        "\nNiche1 breadth\tMean:", sp.summary$niche1.breadth.mean, "\tSD:", sp.summary$niche1.breadth.sd,
        "\nNiche1 species\tMin: ", sp.summary$niche1.sp.min, "\tMax:", sp.summary$niche1.sp.max,
        "\nTotal breadth:", (sp.summary$niche1.sp.max - sp.summary$niche1.sp.min),
        "\tTotal amount:", sp.summary$total_amount, "\n")
  }

  return()
}