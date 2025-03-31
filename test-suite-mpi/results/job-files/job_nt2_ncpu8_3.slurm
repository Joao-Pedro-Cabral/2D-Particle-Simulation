#!/usr/bin/env bash
#SBATCH --job-name=g24
#SBATCH --output=../mpi_outputs/output_nt2_ncpu8_3.txt
#SBATCH --error=../mpi_outputs/output_nt2_ncpu8_3.txt
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --exclusive
#SBATCH --exclude=lab5p[1-20]
#SBATCH --ntasks-per-node=1
srun parsim-mpi 8555 0.05 3 10 10
