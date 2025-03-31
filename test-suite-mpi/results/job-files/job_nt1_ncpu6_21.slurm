#!/usr/bin/env bash
#SBATCH --job-name=g24
#SBATCH --output=../mpi_outputs/output_nt1_ncpu6_21.txt
#SBATCH --error=../mpi_outputs/output_nt1_ncpu6_21.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --exclusive
#SBATCH --exclude=lab5p[1-20]
#SBATCH --ntasks-per-node=1
srun parsim-mpi 1 1000 3 10000 10000
