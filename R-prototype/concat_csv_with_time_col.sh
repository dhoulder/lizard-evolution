#!/bin/bash

# this script conmbined the csv files created for single timesteps into a single combined file for the full run

OutFileName=$1           			# Set the output name
TempFileName="outTemp.csv"
Timestep_colname="timestep,"
i=0                                       # Set a counter

for filename in out*.csv
do 
	if [ "$filename"  != "$OutFileName" ] ;      # Avoid recursion 
	then 
 
		# add a timestep column to the file (later add a check to ensure this is only added once
		time_from_filename=$(echo $filename | tr -dc '0-9')
		#echo $time_from_filename
		#awk -F, -v OFS=, '{$1='$time_from_filename';$2=k}1' $filename > outnew.csv
		awk -F, -v OFS=, '{$1='$time_from_filename' FS $1;}1' $filename > $TempFileName

		if [[ $i -eq 0 ]]
		then
		
			header_orig=$(head -1 $filename)
			echo "$Timestep_colname $header_orig" > $OutFileName 

		fi
		tail -n +2  $TempFileName >> $OutFileName # Append from the 2nd line each file

	fi
		i=$(( $i + 1 ))                        		# Increase the counter
done

rm $TempFileName
