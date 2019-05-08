#!/bin/sh
# Make wave modeling programs
# if the Inter compiler installed, the instruction is :
all:
	mpicc -o elas Elastic_2m2_de2D.c -lm -fopenmp
	nohup mpirun -f mpd.hosts -np 4 ./elas &

## if the Inter compiler has not installed, the instruction is :
# 	$(MPICC) $(srcfile) -I$(IncludeDir) -I$(lib) -L$(MpiLibDir) -o $(exe) -lm -fopenmp
#run:
#	nohup $(MPIRUN) -np $(num_process) -machinefile $(MachineFile) $(exe) &
# Note taht you should replace the folder path with your own envirment path
