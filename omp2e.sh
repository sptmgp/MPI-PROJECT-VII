#BSUB -e jer-%J-2pc_omp.err.log
#BSUB -o jer-%J-2pc_omp.out.log
#BSUB -J    jer-2pc_omp.job
#BSUB -q normal
#BSUB -n 2
##BSUB -R span[ptile=10]
module purge 
module load numa/2.0.11
#module load cuda/8.0
cd /home/HPC-YT-JAN-MAY-2109/BRYAN-GOMEZ/project_collision

export OMP_NUM_THREADS=2

g++ -std=c++11 -fopenmp particle_collisions_omp.cpp -o omp
./omp 2
