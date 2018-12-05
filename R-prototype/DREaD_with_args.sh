#!/bin/csh
#PBS -P xe5
#PBS -N r2d4a1f40
#PBS -q normalbw
#PBS -l walltime=48:00:00
#PBS -l mem=16GB
#PBS -l ncpus=1
#PBS -j oe

echo "Run number:" $run_num

module purge
module load R/3.4.3

cd /short/ka2/dfr805/simulation/test_runs/real225_nicherate0.02_disp4_amp1_freq40_200gens_2r

Rscript --vanilla ~/code/DREaD_ds/Dread_save_raster_raijin_5temp.r $run_num > output$run_num.txt
