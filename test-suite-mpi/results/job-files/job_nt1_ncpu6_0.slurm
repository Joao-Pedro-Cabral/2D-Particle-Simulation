#!/usr/bin/env bash
#SBATCH --job-name=g24
#SBATCH --output=../mpi_outputs/output_nt1_ncpu6_0.txt
#SBATCH --error=../mpi_outputs/output_nt1_ncpu6_0.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --exclusive
#SBATCH --exclude=lab5p[1-20]
#SBATCH --ntasks-per-node=1
srun parsim-mpi 5893 0.05 3 10 10
