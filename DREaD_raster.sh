#!/bin/csh
#PBS -P ka2
#PBS -N Al_r200
#PBS -q normalbw
#PBS -l walltime=30:00:00
#PBS -l mem=16GB
#PBS -l ncpus=1
#PBS -j oe

echo "Run number:" $run_num
echo "Dispersal:"  $dispersal
echo "Timesteps:"  $timesteps
echo "Directory:"  $dir

module purge
module load R/3.4.3

cd ~/code/DREaD_ds/

Rscript --vanilla ~/code/DREaD_ds/Dread_save_raster_raijin.r "$dispersal" "$timesteps" "$run_num" "$dir" > $dir/output"$run_num".txt
