#!/bin/bash

# Run all serial tests and save their execution times
SERIAL_OUT="results/serial_results.txt"
echo -e "========== Serial Version ==========\n" > "./$SERIAL_OUT"

cd ../serial
make profile

#declare -A serial_times

#TODO: FIX readarray -t TESTS < <(grep -v '^\s*#' ../test-suite-mpi/TESTS.txt)
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

for ((i=0; i<${#TESTS[@]}; i+=3)); do
    params="${TESTS[i]}"
    expected_coords="${TESTS[i+1]}"
    expected_collisions="${TESTS[i+2]}"

    # Append the execution time (stderr) after the two lines output (stdout)
    exec 3>&1 4>&2
    exec 2>temp.err
    output=$(./parsim $params)
    echo "$ ./parsim $params" >>  ../test-suite-mpi/$SERIAL_OUT
    exec 2>&4 4>&-
    exec_time=$(cat temp.err | cut -d's' -f1)
    rm temp.err

    echo "$output" >> ../test-suite-mpi/$SERIAL_OUT
    echo "$exec_time" >> ../test-suite-mpi/$SERIAL_OUT
    coords=$(echo "$output" | head -n 1)
    collisions=$(echo "$output" | head -n 2 | tail -n 1)

   # serial_times[$params]=$execution_time
    
    # Check if output matches expected values
    if [ "$coords" != "$expected_coords" ] || [ "$collisions" != "$expected_collisions" ]; then
        echo "Test failed for $cmd $params "
        echo "Expected coordinates: $expected_coords, got: $coords"
        echo "Expected collisions: $expected_collisions, got: $collisions"
        exit 1
    fi

done

echo "Serial tests completed and saved to test-suite-mpi/$SERIAL_OUT"
