#!/usr/bin/env bash
#SBATCH --job-name=g24
#SBATCH --output=../mpi_outputs/output_nt4_ncpu2_6.txt
#SBATCH --error=../mpi_outputs/output_nt4_ncpu2_6.txt
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --exclusive
#SBATCH --exclude=lab5p[1-20]
#SBATCH --ntasks-per-node=1
srun parsim-mpi 12 100 5 10000 10000
