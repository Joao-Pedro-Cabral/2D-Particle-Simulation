#!/usr/bin/env bash
#SBATCH --job-name=g24
#SBATCH --output=../mpi_outputs/output_nt1_ncpu8_24.txt
#SBATCH --error=../mpi_outputs/output_nt1_ncpu8_24.txt
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --exclude=lab5p[1-20]
#SBATCH --ntasks-per-node=1
srun parsim-mpi 3 5000 50 1000000 300
