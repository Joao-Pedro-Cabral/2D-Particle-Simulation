for f in *.slurm; do
    sbatch "$f" 
done 