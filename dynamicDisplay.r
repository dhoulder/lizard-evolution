# functions for displaying infoirmation in runtime such as 
# a map and summary stats

display.initialise <- function() {
  windows(15,11)
  
  my.colours.function <- colorRampPalette(colors = c("red", "yellow", "blue"))
  my.colours    <- my.colours.function(250)
  return(my.colours)
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
    plot(plotItems[["env"]], main=main.header)
  }
  
  if (length(plotItems[["demes_amount_position"]]) > 0) {
    demetable.species <- plotItems[["demes_amount_position"]]
    these.colours <- my.colours[round(demetable.species$niche1.position*10)]
    these.sizes   <- sqrt(demetable.species$amount) * 2
    points(demetable.species$x, demetable.species$y, bg=these.colours, pch=21, fg="black", cex=these.sizes)
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
        "\tTotal breadth:", (sp.summary$niche1.sp.max - sp.summary$niche1.sp.min), "\n")
  } 
  
  return()
}