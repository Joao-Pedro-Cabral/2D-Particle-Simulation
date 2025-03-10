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

function check_output() {
    local cmd=$1
    local params=$2
    local expected_coords=$3
    local expected_collisions=$4

    echo "$ $cmd $params" >> output.txt
    output=$($cmd $params)
    echo "$output" >> output.txt

    # Read the first two lines of output
    coords=$(echo "$output" | head -n 1)
    collisions=$(echo "$output" | head -n 2 | tail -n 1)

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
    check_output "./parsim" "${TESTS[i]}" "${TESTS[i+1]}" "${TESTS[i+2]}"
done

echo "OMP version" >> output.txt
cd ../omp
make profile

# Run OMP tests
for ((i=0; i<${#TESTS[@]}; i+=3)); do
    check_output "./parsim-omp" "${TESTS[i]}" "${TESTS[i+1]}" "${TESTS[i+2]}"
done