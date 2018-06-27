#!/bin/sh
# Make wave modeling programs
all:
	mpicc -o ac Acoustic_2m4_de1D_f.c -lm -fopenmp
	nohup mpirun -f mpd.hosts -np 16 ./ac &

