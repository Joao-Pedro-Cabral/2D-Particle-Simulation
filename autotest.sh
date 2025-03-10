#!/bin/bash

# Test parameters and expected results (params coords collisions)
TESTS=(
    "1 5000 100 1000000 100" "3899.787 156.291" "163"
    "1 5000 20 1000000 10" "3918.912 143.364" "19"
    "1 1000 3 10000 10000" "287.788 261.446" "31"
    "3 5000 50 1000000 300" "3819.032 25.659" "469"
    "3 5000 50 1000000 500" "3738.436 58.743" "804"
    "-1 1000 30 100000 1000" "575.878 370.663" "1203"
)

declare -A serial_times
declare -A omp_times

function check_output() {
    local cmd=$1
    local params=$2
    local expected_coords=$3
    local expected_collisions=$4
    local version=$5

    echo "$ $cmd $params" >> output.txt
    output=$($cmd $params)
    echo "$output" >> output.txt

    # Read the output
    coords=$(echo "$output" | head -n 1)
    collisions=$(echo "$output" | head -n 2 | tail -n 1)
    execution_time=$(echo "$output" | head -n 3 | tail -n 1 | cut -d's' -f1)

    # Store execution time based on version
    if [ "$version" == "serial" ]; then
        serial_times[$params]=$execution_time
    else
        omp_times[$params]=$execution_time
        # Calculate and print speedup
        serial_time=${serial_times[$params]}
        speedup=$(echo "scale=2; $serial_time / $execution_time" | bc)
        echo "Speedup for test [$params]: $speedup" >> output.txt
    fi

    # Check if output matches expected values
    if [ "$coords" != "$expected_coords" ] || [ "$collisions" != "$expected_collisions" ]; then
        echo "Test failed for $cmd $params!"
        echo "Expected coordinates: $expected_coords, got: $coords"
        echo "Expected collisions: $expected_collisions, got: $collisions"
        exit 1
    fi
}

echo "Serial version" >> output.txt
cd serial
make profile

# Run serial tests
for ((i=0; i<${#TESTS[@]}; i+=3)); do
    check_output "./parsim" "${TESTS[i]}" "${TESTS[i+1]}" "${TESTS[i+2]}" "serial"
done

echo "OMP version" >> output.txt
cd ../omp
make profile

# Run OMP tests and compute speedups
for ((i=0; i<${#TESTS[@]}; i+=3)); do
    check_output "./parsim-omp" "${TESTS[i]}" "${TESTS[i+1]}" "${TESTS[i+2]}" "omp"
done

# Print summary of all speedups
echo -e "\nSpeedup Summary:" >> output.txt
for params in "${!serial_times[@]}"; do
    speedup=$(echo "scale=2; ${serial_times[$params]} / ${omp_times[$params]}" | bc)
    echo "Test [$params]: $speedup" >> output.txt
done