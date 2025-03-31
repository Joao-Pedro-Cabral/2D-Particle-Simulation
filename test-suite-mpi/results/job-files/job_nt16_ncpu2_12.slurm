#!/usr/bin/env bash
#SBATCH --job-name=g24
#SBATCH --output=../mpi_outputs/output_nt16_ncpu2_12.txt
#SBATCH --error=../mpi_outputs/output_nt16_ncpu2_12.txt
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --exclude=lab5p[1-20]
#SBATCH --ntasks-per-node=1
srun parsim-mpi 1 5000 100 1000000 4
