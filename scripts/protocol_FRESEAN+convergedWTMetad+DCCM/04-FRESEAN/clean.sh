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

rm *.mmat *mmat.dat *xyz covar.inp extract.inp *.out
