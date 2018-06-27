#!/bin/sh
# Make wave modeling programs
all:
	mpicc -o ac Acoustic_2m2_de2D.c -lm -fopenmp
	nohup mpirun -f mpd.hosts -np 4 ./ac &

