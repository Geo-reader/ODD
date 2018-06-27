#!/bin/sh
# Make wave modeling programs
all:
	mpicc -o elas Elastic_2m2_de2D.c -lm -fopenmp
	nohup mpirun -f mpd.hosts -np 4 ./elas &

