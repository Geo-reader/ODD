# ODD
Our package aims to provide a set of 3D parallel simulation schemes based on the optimal decomposition strategies for acoustic and elastic wave equations.
---
title: Manual of code package
date: 2018-06-20
author: Ning Wang
mathjax: true
---


## Overview of the optimal domain decomposition strategies

The domain decomposition strategies include decomposition direction and decomposition dimension. As a matter of fact, domain decomposition strategies can significantly affect the computational performance. More specifically, decomposition direction determines the cache hit ratio during register addressing, while decomposition dimension controls the load balance and communication among each node.

In this paper, we consider three domain decomposition schemes: 1D decomposition perpendicular to the fastest dimension (denoted as 1D_f), 1D decomposition perpendicular to the slowest dimension (denoted as 1D_s), 2D decomposition perpendicular to the slowest and the second slowest dimensions. The slowest dimension is the dimension that data are continuously stored. The comparison of 1D_f decomposition method and 1D_s decomposition method is designed to illustrate the effect of decomposition direction on computational efficiency; and the comparison of 1D_s decomposition method and 2D decomposition method is designed to illustrate the effect of decomposition dimension on computational efficiency.

Conclusions about the impact of domain decomposition strategy on parallel computational performance are drawn via theoretical analysis and verification of four metrics including computation time, speedup ratio, strong scaling property and memory usage. They can be summarized as follows. First, decomposition direction has a dramatic influence on compute time, speedup ratio and strong scaling property. It should be given priority to decompose the model perpendicular to the slowest dimension, in which data are continuously saved. Second, the impact of decomposition dimension is more evident in the fine-grained parallelism. Two-dimensional domain decomposition performs better when the model is finely decomposed. Third, compared to the 1D domain decomposition, 2D domain decomposition can make more efficient use of memory. In general, 2D domain decomposition of the model is a better choice in most cases.

## The architecture of code package 

Our package aims to provide a set of 3D parallel simulation schemes based on the optimal decomposition strategies for acoustic and elastic wave equations. Three scenarios based on optimized 2D decomposition scheme are included in the `Optimal version` of our package.

1. **Acoustic_2m4_de2D** is used to run the time-space-domain high-order (2M, 4) acoustic simulation based on the optimized 2D decomposition scheme. 
   
2. **Acoustic_2m2_de2D** is used to run the traditional (2M, 2) acoustic simulation based on the optimized 2D decomposition scheme.
   
3. **Elastic_2m2_de2D** is used to run the traditional (2M, 2) elastic simulation based on the optimized 2D decomposition scheme.

If you want to verify the conclusions about domain decomposition strategies, we also provide a `Verification version`, which is consisted of following three scenarios.  

1. **Acoustic_2m4_de2D** is used to run the time-space-domain high-order (2M, 4) acoustic simulation based on the optimized 2D decomposition scheme.

2. **Acoustic_2m4_de1D_f** is used to run the time-space-domain high-order (2M, 4) acoustic simulation. The model is decomposed perpendicular to the fastest dimension. 

3. **Acoustic_2m4_de1D_s** is used to run the time-space-domain high-order (2M, 4) acoustic simulation. The model is decomposed perpendicular to the slowest dimension.

Next, we take the **Acoustic_2m4_de2D** as an example to outline the architect of our code. 

-   `input`: scripts for generating and decomposing models of velocity and density:
    - `gene_models.m`: generating velocity and density models;
    - `split.cpp`: decomposing velocity and density models;
    - `Makefile`: execution script of split.cpp.
-   `output`: scripts for storing generated sub- seismograms and splicing them:
    - `splice.cpp`: splicing sub-seismograms into a whole seismogram;
    - `Makefile`: execution script of splice.cpp;
    - `plot_rec.m`: plot the whole seismogram.
-   `Myfunctions.h`: header file;
-   `Acoustic_2m2_de2D.c`: forward code file;
-   `Makefile`: excution script;
-   `geometry.txt`: parameters of geometry;
-   `Parameters3d.txt`: some parameters of forward modeling.

## Prerequisites

our package is developed under `Linux` system, which should be equipped with the following environments:

- MPI environment (for example, mpich2-1.0.2.p1 version);
- matlab.

## How to run this package

We use **Acoustic_2m4_de2D** as an example to demonstrate how to run this package. If you want to quick test this package, you can generate a small 3D model with the size of 101^3 in `Optimal version` of our package. 

- Step 1: Run the matlab file `gene_models.m` in `./input` to generate velocity and density models;
- Step 2: Run the `split.cpp` in `./input` to decompose the velocity and the density models. The `Makefile` is written as  
``` bash
#!/bin/sh

# Decomposing velocity and density models

all:
	icc -w -o spl split.cpp -lm
	./spl
```
- Step 3: Run the `Acoustic_2m4_de2D.c`. The `Makefile` is written as
``` bash
#!/bin/sh

# Make wave modeling programs

all:
	mpicc -o ac Acoustic_2m2_de2D.c -lm -fopenmp
	nohup mpirun -f mpd.hosts -np 4 ./ac &
```
- Step 4: When the `Acoustic_2m4_de2D.c` is executed, view the time-consuming by the command line: `vi nohup`

- Step 5: Run the `splice.cpp` in `./output` to splice the sub-seismograms. The `Makefile` is written as  
``` bash
#!/bin/sh

# Splicing sub- seismograms into a whole seismogram

all:
	icc -w -o spli splice.cpp -lm
	./spli
```
- Step 6: Run the `plot_rec.m` in `./output` to plot the seismic record.

## Some points should be explained

1. Our code exploits hybrid MPI and OpenMp, which means the parallelism between nodes is completed by MPI, and the parallelism within a node is completed by OpenMP. It is strongly recommended to set the number of compute nodes equal to the number of sub-domains. For example, in the `Optimal version`, 4 compute nodes are required if the models are decomposed into 2 sub-domains at x- and y- directions, respectively.

2. For quick test our optimal code, we initially set up small 101^3 3D models calculated by 4 compute nodes in `Optimal version`(models are decomposed into 2 sub-domains at x- and y- directions, respectively). In order to highlight the difference in computational efficiency between three decomposition schemes, we initially set up 512^3 3D models calculated by 16 compute nodes in `Verification version`.

3. If you want to run our code to simulate other models, parameters in `Parameters3d.txt`, `geometry.txt`, `split.cpp`, `splice.cpp`, `gene_models.m` and `plot_rec.m` should also make corresponding changes.


4. Due to differences in hardware configuration, the wall-clock time of your cluster may be different from the wall-clock time in this paper. But the conclusions (comparison of three domain decomposition schemes) can be verified.

5. When you test the compute time, to avoid interference, please make sure that no other process is simultaneously executed on your compute nodes.

## Contact me

I am Ning Wang, a PhD candidate from China University of Petroleum, Beijing. If you have any question about this coda package, please feel free to contact me by wangn2311@163.com.

## Copyright

 
Copyright (C) 2018  China University of Petroleum, Beijing (Ning Wang)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
