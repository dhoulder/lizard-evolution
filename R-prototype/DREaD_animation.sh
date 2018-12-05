#!/bin/csh
#PBS -P ka2
#PBS -N Al_an_250
#PBS -q normalbw
#PBS -l walltime=48:00:00
#PBS -l mem=24GB
#PBS -l ncpus=1
#PBS -j oe

echo "Run number:" $run_num
echo "Dispersal:"  $dispersal
echo "Timesteps:"	 $timesteps
echo "Directory:"  $dir

module purge
module load R/3.4.3

cd ~/code/DREaD_ds/

Rscript --vanilla ~/code/DREaD_ds/DREaD_save_images_raijin.r "$dispersal" "$timesteps" "$dir" > $dir/output.txt
