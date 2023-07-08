#!/bin/bash

# Set the process ID to wait for
process_id="130791"

# Function to wait for process to end
wait_for_process() {
    while ps -p "$1" > /dev/null; do
        sleep 1
    done
}

# Wait for the process to end
wait_for_process "$process_id"

mpicc ../examples/"script test.c" -o scri -L../make -l:libphysics.a -fopenmp -lm -lc
mpiexec -n 1 ./scri