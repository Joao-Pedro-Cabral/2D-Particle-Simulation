#!/usr/bin/env bash

# Test parameters and expected results (params coords collisions)
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

# Thread counts to test
NTASKS=(1 2 4 8 16)
CPUS_PER_TASK=(1 2 4 6)

declare -A serial_times
declare -A mpi_times

function check_output() {
    local cmd=$1
    local params=$2
    local expected_coords=$3
    local expected_collisions=$4
    local version=$5
    local ntasks=$6
    local cpus_per_task=$7
    local output_file=$8

    echo "$cmd $params (ntasks=$ntasks, cpus-per-task=$cpus_per_task)" >> ../output-mpi.txt

    if [ "$version" == "serial" ]; then
        exec 3>&1 4>&2
        exec 2>temp.err
        output=$($cmd $params)
        exec 2>&4 4>&-
        execution_time=$(cat temp.err)
        rm temp.err
    else 
        output=$(cat "$output_file")
        execution_time=$(grep 's$' "$output_file" | tail -n 1 | cut -d's' -f1)
    fi

    echo "$output" >> ../output-mpi.txt
    echo "$execution_time" >> ../output-mpi.txt
    coords=$(echo "$output" | head -n 1)
    collisions=$(echo "$output" | head -n 2 | tail -n 1)
    
    if [ "$version" == "serial" ]; then
        serial_times[$params]=$execution_time
    else
        key="${params}_${ntasks}_${cpus_per_task}"
        mpi_times[$key]=$execution_time
        serial_time=${serial_times[$params]}
        speedup=$(echo "scale=2; $serial_time / $execution_time" | bc)
        echo "Speedup for test [$params] with $ntasks tasks and $cpus_per_task CPUs per task: $speedup" >> ../output-mpi.txt
    fi

    # Check if output matches expected values
    if [ "$coords" != "$expected_coords" ] || [ "$collisions" != "$expected_collisions" ]; then
        echo "Test failed for $cmd $params "
        echo "Expected coordinates: $expected_coords, got: $coords"
        echo "Expected collisions: $expected_collisions, got: $collisions"
        exit 1
    fi
}

echo -e "========== Serial version ==========\n" > output-mpi.txt
cd serial
make profile

for ((i=0; i<${#TESTS[@]}; i+=3)); do
    params="${TESTS[i]}"
    expected_coords="${TESTS[i+1]}"
    expected_collisions="${TESTS[i+2]}"
    check_output "./parsim" "$params" "$expected_coords" "$expected_colllision" "serial" "1"
done

echo -e "\n\n========== MPI version ==========" >> ../output-mpi.txt
cd ../mpi
make profile
for ntasks in "${NTASKS[@]}"; do
    for cpus_per_task in "${CPUS_PER_TASK[@]}"; do
        for ((i=0; i<${#TESTS[@]}; i+=3)); do
            params="${TESTS[i]}"
            expected_coords="${TESTS[i+1]}"
            expected_collisions="${TESTS[i+2]}"
            job_name="mpi_ntasks_${ntasks}_cpus_${cpus_per_task}_test_${i}"
            output_file="../output-mpi-${job_name}.txt"

            # Create a Slurm job script
            job_script="job_${job_name}.slurm"
            echo "#!/bin/bash" > $job_script
            echo "#SBATCH --job-name=$job_name" >> $job_script
            echo "#SBATCH --output=$output_file" >> $job_script
            echo "#SBATCH --ntasks=$ntasks" >> $job_script
            echo "#SBATCH --cpus-per-task=$cpus_per_task" >> $job_script
            echo "#SBATCH --exclusive" >> $job_script
            echo "#SBATCH --ntasks-per-node=1" >> $job_script
            echo "" >> $job_script
            echo "srun ./parsim-mpi $params" >> $job_script

            # Submit the job and wait for it to complete
            job_id=$(sbatch $job_script | awk '{print $4}')
            while squeue -j $job_id > /dev/null 2>&1; do
                sleep 1
            done

            # Check the output
            check_output "$params" "$expected_coords" "$expected_collisions" "mpi" "$ntasks" "$cpus_per_task" "$output_file"
        done
    done
done

# Print summary of all speedups
echo -e "\n========== Speedup Summary ==========\n" >> ../output-mpi.txt
for ((i=0; i<${#TESTS[@]}; i+=3)); do
    params="${TESTS[i]}"
    echo "Test [$params]:" >> ../output-mpi.txt
    for ntasks in "${NTASKS[@]}"; do
        for cpus_per_task in "${CPUS_PER_TASK[@]}"; do
            key="${params}_${ntasks}_${cpus_per_task}"
            if [ -n "${mpi_times[$key]}" ]; then
                speedup=$(echo "scale=2; ${serial_times[$params]} / ${mpi_times[$key]}" | bc)
                echo "  $ntasks tasks, $cpus_per_task CPUs per task: $speedup" >> ../output-mpi.txt
            fi
        done
    done
done