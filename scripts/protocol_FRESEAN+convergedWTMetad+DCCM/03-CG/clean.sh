#!/bin/bash

wd=`pwd`

echo "This will delete your output data in ${wd}"
read -p "Do you want to proceed? (y/n) " yn

case $yn in 
	y ) echo ok, we will proceed;;
	n ) echo exiting...;
		exit;;
	* ) echo invalid response;
		exit 1;;
esac

rm *.mtop mtop.out trjconv*.out coarse.inp coarse.out ref.gro ref.pdb ref-cg.gro *.trr slurm-*.out 
