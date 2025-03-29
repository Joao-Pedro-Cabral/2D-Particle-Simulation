#!/usr/bin/env bash

# Generate and submit MPI jobs using Slurm
TESTS=(
    "5893 0.05 3 10 10" "0.002 0.035" "2"
    "8555 0.05 3 10 10" "0.016 0.049" "1"
    "12 100 5 10000 10000" "76.732 61.943" "2209"
    "-11 3500 20 500000 10" "1984.878 1625.992" "35"
    "1 5000 100 1000000 4" "3936.506 131.472" "4"
    "1 5000 100 1000000 100" "3899.787 156.291" "163"
    "1 5000 20 1000000 10" "3918.912 143.364" "19"
    "1 1000 3 10000 10000" "287.788 261.446" "31"
    "3 5000 50 1000000 300" "3819.032 25.659" "469"
    "3 5000 50 1000000 500" "3738.436 58.743" "804"
    "-1 1000 30 100000 1000" "575.878 370.663" "1203"
    "32 100 6 10000 100000" "29.530 30.082" "4448"
)

NTASKS=(1 2 4 8 16 32)
CPUS_PER_TASK=(1 2 4 6 8)

rm -rf results/mpi_outputs
mkdir -p results/mpi_outputs

cd ../mpi || exit 1
make profile || exit 1

for ntasks in "${NTASKS[@]}"; do
    for cpus_per_task in "${CPUS_PER_TASK[@]}"; do
        # number of machine used must be atmost 32
        # if ((ntasks*cpus_per_task > 32)); then
        #     continue
        # fi
        for ((i=0; i<${#TESTS[@]}; i+=3)); do
            params="${TESTS[i]}"
            # If ncside < ntasks then jump the test
            ncside=$(echo "$params" | awk '{print $3}')
            if (( ncside < ntasks )); then
                continue
            fi
            job_name="g24" 
            output_file="../test-suite-mpi/results/mpi_outputs/output_nt${ntasks}_ncpu${cpus_per_task}.txt"
            #stderr_file="../test-suite-mpi/results/mpi_outputs/output_nt${ntasks}_ncpu${cpus_per_task}.txt"
            job_script="job.slurm"

            cat <<EOF > "$job_script"
#!/usr/bin/env bash
#SBATCH --job-name=$job_name
#SBATCH --output=$output_file
#SBATCH --error=$output_file
#SBATCH --ntasks=$ntasks
#SBATCH --cpus-per-task=$cpus_per_task
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
echo "$ ./parsim-mpi ${params}" >> $output_file
srun parsim-mpi $params
EOF

            sbatch "$job_script"
        done
    done
done

echo "All MPI jobs submitted."
