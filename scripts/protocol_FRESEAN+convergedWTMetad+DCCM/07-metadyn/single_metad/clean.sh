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

rm *.tpr *.trr *.edr *.log *.cpt *.out *.hills *.fes mdout.mdp plumed-mode-input.pdb
