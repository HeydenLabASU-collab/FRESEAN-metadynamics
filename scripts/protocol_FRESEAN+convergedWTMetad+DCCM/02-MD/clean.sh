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

rm grompp.out mdrun.out mdout.mdp sample-NPT.edr sample-NPT.log sample-NPT.tpr sample-NPT.trr slurm-*.out
