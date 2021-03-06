#!/bin/bash
shopt -s globstar

home_dir=/short/ka2/dfr805/simulation/test_runs/

cd $home_dir

timesteps=200;

#for D in 1 1.5 2 3 5
for D in 1
do

	for run_num in {1..100}
	do

		echo $D  $timesteps $run_num;
	
		#new_dir=$home_dir'Alps225_nicherate0.02_envStatic_disp'$D'_'$timesteps'gens_r'
		new_dir=$home_dir'Alps225_nicherate0.02_amp1.5_freq50_disp'$D'_'$timesteps'gens_r'
		rm -r $new_dir/
		mkdir -p $new_dir
		
		qsub -v dispersal=$D,timesteps=$timesteps,run_num=$run_num,dir=$new_dir -N "Adv"$D"_"$run_num ~/code/DREaD_ds/DREaD_raster.sh
	
	done

done
