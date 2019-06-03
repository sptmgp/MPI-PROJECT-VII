#BSUB -e jer-%J-4pc_mpi.err.log
#BSUB -o jer-%J-4pc_mpi.out.log
#BSUB -J    jer-4pc_mpi.job
#BSUB -q normal
#BSUB -n 4
##BSUB -R span[ptile=10]
module purge 
module load openmpi/1.8.8
module load numa/2.0.11
#module load cuda/8.0
cd /home/HPC-YT-JAN-MAY-2109/BRYAN-GOMEZ/project_collision
mpirun -np 4 mpi 2
