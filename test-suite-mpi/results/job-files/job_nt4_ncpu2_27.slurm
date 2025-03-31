#!/usr/bin/env bash
#SBATCH --job-name=g24
#SBATCH --output=../mpi_outputs/output_nt4_ncpu2_27.txt
#SBATCH --error=../mpi_outputs/output_nt4_ncpu2_27.txt
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --exclude=lab5p[1-20]
#SBATCH --ntasks-per-node=1
srun parsim-mpi 3 5000 50 1000000 500
