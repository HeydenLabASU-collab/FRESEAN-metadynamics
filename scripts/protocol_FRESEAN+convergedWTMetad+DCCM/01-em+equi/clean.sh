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

rm -r em equi slurm-*.out
