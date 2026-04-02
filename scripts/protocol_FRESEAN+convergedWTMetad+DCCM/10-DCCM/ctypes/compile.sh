#!/bin/bash

gcc -O3 -fopenmp -fpic -c unwrap.c
gcc -shared -lgomp unwrap.o -o libunwrap.so

gcc -O3 -fopenmp -fpic -c align.c
gcc -shared -lgomp align.o -o libalign.so