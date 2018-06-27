#!/bin/sh
# Make wave modeling programs
all:
	mpicc -o ac Acoustic_2m4_de2D.c -lm -fopenmp
	nohup mpirun -f mpd.hosts -np 4 ./ac &

