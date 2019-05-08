#!/bin/sh
# Make wave modeling programs
all:
	mpicc -o ac Acoustic_2m4_de1D_s.c -lm -fopenmp
	nohup mpirun -f mpd.hosts -np 16 ./ac &
	
## if the Inter compiler has not installed, the instruction is :
# 	$(MPICC) $(srcfile) -I$(IncludeDir) -I$(lib) -L$(MpiLibDir) -o $(exe) -lm -fopenmp
#run:
#	nohup $(MPIRUN) -np $(num_process) -machinefile $(MachineFile) $(exe) &
# Note taht you should replace the folder path with your own envirment path	

