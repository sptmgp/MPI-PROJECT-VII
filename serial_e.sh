#!/bin/bash

#BSUB -e jer-%J-pc_serial.err.log
#BSUB -o jer-%J-pc_serial.out.log
#BSUB -J jer-pc_serial.job
#BSUB -q normal
#BSUB -n 1
module purge 
# module load openmpi/1.8.8
module load numa/2.0.11
module load gcc/5.4.0 
#module load cuda/8.0
cd /home/HPC-YT-JAN-MAY-2109/BRYAN-GOMEZ/project_collision
g++  -std=c++11 particle_collisions_serial.cpp -o serial
./serial 2
