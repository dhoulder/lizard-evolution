require(data.table)
require(dplyr)
require(raster)
#library(yaml)

# get external arguments passed to this script
input.args	        <- commandArgs(trailingOnly = TRUE)
run_num        	    <- as.numeric(input.args[1])
working_dir    	    <- as.numeric(input.args[2])
timesteps				    <- as.numeric(input.args[3])
all_times_filename  <- as.numeric(input.args[4])
raster_template     <- as.numeric(input.args[5])

# TEMP PARAMETERS AND STEPS FOR DEV
rm(list=ls())
run_num=1
working_dir="C:/Users/u3579238/Work/Simulation/scratch/"
timesteps=100
all_times_filename="combined_all_times.csv"
raster_template="C:/Users/u3579238/Work/Software/dan-github/DREaD_extras/realAlps250_rescaled_100_bio01.asc"

setwd(working_dir)

# create the raster template
template.ras <- raster(raster_template)
template.ras[] <- 0

# columns to skip (to save time / memory)
skip_cols <- c("niche_breadth_0", "niche_breadth_1")

# read in the csv file of all occurrences and timesteps
all_occs <- fread(all_times_filename, drop=skip_cols)

# get the final value per cell and write to csv or raster
final_occs <- all_occs[timestep==timesteps]
cells <- cellFromRowCol(template.ras, rownr = final_occs$row, colnr = final_occs$column)
final_amounts.ras <- template.ras
final_amounts.ras[cells] <- final_occs$amount
raster_filename <- paste("cellFinal_", run_num, ".asc", sep='')
writeRaster(final_amounts.ras, raster_filename, overwrite=T)

# get the mean value per cell and write to csv or raster
# don't use mean, but sum and divide by timesteps because this way empty cells are treated as 0
cell_sums <- all_occs[, .(amount_sum=sum(amount)), by=list(row, column)]
cell_sums[, amount_mean:=amount_sum / timesteps]

cells <- cellFromRowCol(template.ras, rownr = cell_sums$row, colnr = cell_sums$column)
cell_mean.ras <- template.ras
cell_mean.ras[cells] <- cell_sums$amount_mean
raster_filename <- paste("cellSums_", run_num, ".asc", sep='')
writeRaster(cell_mean.ras, raster_filename, overwrite=T)

# calculate environment mismatch as
# a) absolute difference between env_0 and niche_centre_0,  env_1 and niche_centre_1
# b) actual difference, niche_centre_0 - env_0, niche_centre_1 - env_1
all_occs[,env_0_diff:=env_0 - niche_centre_0]
all_occs[,env_1_diff:=env_1 - niche_centre_1]
all_occs[,env_0_absdiff:=abs(env_0 - niche_centre_0)]
all_occs[,env_1_absdiff:=abs(env_1 - niche_centre_1)]

cell_envdiff <- all_occs[, .(env_0_diff=mean(env_0_diff), 
                             env_1_diff=mean(env_1_diff), 
                             env_0_absdiff=mean(env_0_absdiff),
                             env_1_absdiff=mean(env_1_absdiff)),
                         by=list(row, column)]

cells <- cellFromRowCol(template.ras, rownr = cell_envdiff$row, colnr = cell_envdiff$column)
cell_env_0_diff.ras    <- template.ras
cell_env_1_diff.ras    <- template.ras
cell_env_0_absdiff.ras <- template.ras
cell_env_1_absdiff.ras <- template.ras

cell_env_0_diff.ras[cells]   <- cell_envdiff$env_0_diff
cell_env_1_diff.ras[cells]   <- cell_envdiff$env_1_diff
cell_env_0_absdiff.ras[cells]   <- cell_envdiff$env_0_absdiff
cell_env_1_absdiff.ras[cells]   <- cell_envdiff$env_1_absdiff

raster_filename <- paste("env_0_diff_", run_num, ".asc", sep='')
writeRaster(cell_env_0_diff.ras, raster_filename, overwrite=T)
raster_filename <- paste("env_1_diff_", run_num, ".asc", sep='')
writeRaster(cell_env_1_diff.ras, raster_filename, overwrite=T)
raster_filename <- paste("env_0_absdiff_", run_num, ".asc", sep='')
writeRaster(cell_env_0_absdiff.ras, raster_filename, overwrite=T)
raster_filename <- paste("env_1_absdiff_", run_num, ".asc", sep='')
writeRaster(cell_env_1_absdiff.ras, raster_filename, overwrite=T)