#!/bin/bash

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
)

# Thread counts to test
THREADS=(1 2 4 8)

declare -A serial_times
declare -A omp_times

function check_output() {
    local cmd=$1
    local params=$2
    local expected_coords=$3
    local expected_collisions=$4
    local version=$5
    local threads=$6

    echo "$ $cmd $params (threads: $threads)" >> ../output.txt
    
    # Execution time is redirected in stdrout in parsim.c
    # output=$($cmd $params)
    # echo "$output" >> ../output.txt
    # coords=$(echo "$output" | head -n 1)
    # collisions=$(echo "$output" | head -n 2 | tail -n 1)
    # execution_time=$(echo "$output" | head -n 3 | tail -n 1 | cut -d's' -f1)

    # Execution time is redirected in stderr in parsim.c
    exec 3>&1 4>&2
    exec 2>temp.err
    output=$($cmd $params)
    exec 2>&4 4>&-
    execution_time=$(cat temp.err)
    rm temp.err
    echo "$output" >> ../output.txt
    echo "$execution_time" >> ../output.txt
    coords=$(echo "$output" | head -n 1)
    collisions=$(echo "$output" | head -n 2 | tail -n 1)
    execution_time=$(echo "$execution_time" | cut -d's' -f1)

    

    # Store execution time based on version
    if [ "$version" == "serial" ]; then
        serial_times[$params]=$execution_time
    else
        omp_times["${params}_${threads}"]=$execution_time
        # Calculate and print speedup
        serial_time=${serial_times[$params]}
        speedup=$(echo "scale=2; $serial_time / $execution_time" | bc)
        echo "Speedup for test [$params] with $threads threads: $speedup" >> ../output.txt
    fi

    # Check if output matches expected values
    if [ "$coords" != "$expected_coords" ] || [ "$collisions" != "$expected_collisions" ]; then
        echo "Test failed for $cmd $params "
        echo "Expected coordinates: $expected_coords, got: $coords"
        echo "Expected collisions: $expected_collisions, got: $collisions"
        exit 1
    fi
}

echo -e "========== Serial version ==========\n" > output.txt
cd serial
make profile

for ((i=0; i<${#TESTS[@]}; i+=3)); do
    check_output "./parsim" "${TESTS[i]}" "${TESTS[i+1]}" "${TESTS[i+2]}" "serial" "1"
done

echo -e "\n\n========== OMP version ==========" >> ../output.txt
cd ../omp
make profile

for threads in "${THREADS[@]}"; do
    echo -e "\nRunning with $threads threads:" >> ../output.txt
    export OMP_NUM_THREADS=$threads
    for ((i=0; i<${#TESTS[@]}; i+=3)); do
        check_output "./parsim-omp" "${TESTS[i]}" "${TESTS[i+1]}" "${TESTS[i+2]}" "omp" "$threads"
    done
done

# Print summary of all speedups
echo -e "\nSpeedup Summary:" >> ../output.txt
for ((i=0; i<${#TESTS[@]}; i+=3)); do
    params="${TESTS[i]}"
    echo "Test [$params]:" >> ../output.txt
    for threads in "${THREADS[@]}"; do
        speedup=$(echo "scale=2; ${serial_times[$params]} / ${omp_times[${params}_${threads}]}" | bc)
        echo "  $threads threads: $speedup" >> ../output.txt
    done
done