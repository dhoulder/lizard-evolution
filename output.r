# this file should contain functions to return DREaD results

demes_to_raster <- function(demeTable, speciesID, env) {

  cell_amount_table <- demeTable[speciesID==speciesID, .(cellID, amount)]
  raster.new <- env
  raster.new[] <- 0
  raster.new[cell_amount_table$cellID] <- cell_amount_table$amount
  return(raster.new)
}
