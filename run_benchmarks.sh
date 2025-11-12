#!/bin/bash

# --- Configuration ---
# The name of your compiled program
EXECUTABLE="./LocAG"

# The CSV file where time results will be stored
OUTPUT_FILE="benchmark_results.csv"

# A directory to store the detailed console output of each run
LOG_DIR="benchmark_logs"

# --- End of Configuration ---

# Function to run a single benchmark test
# Usage: run_benchmark <config_name> <array_type> <method> <policy>
run_benchmark() {
    local CONFIG=$1
    local TYPE=$2
    local METHOD=$3
    local POLICY=$4

    # Create a unique name for the log file
    local LOG_FILE="${LOG_DIR}/${CONFIG}_${TYPE}_${METHOD}_${POLICY}.log"
    
    # Create a temporary file to store the output of the 'time' command
    local TIME_FILE=$(mktemp)

    echo "Running: Config=$CONFIG, Type=$TYPE, Method=$METHOD, Policy=$POLICY"

    # Run the command:
    # ( ... ) groups the command so we can time it.
    # > $LOG_FILE redirects the program's standard output to our log file.
    # 2> $TIME_FILE redirects the 'time' command's standard error to our temp time file.
    ( time -p $EXECUTABLE $CONFIG $TYPE $METHOD $POLICY ) > $LOG_FILE 2> $TIME_FILE

    # Check if the program finished successfully (exit code 0)
    if [ $? -ne 0 ]; then
        echo "  ERROR: Run failed. Check $LOG_FILE for details."
        # Still log the failure and times
        local REAL_TIME="ERROR"
        local USER_TIME="ERROR"
        local SYS_TIME="ERROR"
    else
        # Parse the real, user, and sys times from the temp file
        local REAL_TIME=$(grep "real" $TIME_FILE | awk '{print $2}')
        local USER_TIME=$(grep "user" $TIME_FILE | awk '{print $2}')
        local SYS_TIME=$(grep "sys" $TIME_FILE | awk '{print $2}')
        
        echo "  Done. Real time: ${REAL_TIME}s. Log saved to $LOG_FILE"
    fi

    # Append the results to our main CSV file
    echo "$CONFIG,$TYPE,$METHOD,$POLICY,$REAL_TIME,$USER_TIME,$SYS_TIME" >> $OUTPUT_FILE

    # Clean up the temporary time file
    rm $TIME_FILE
}

# --- Main execution ---

# 1. First, (re)compile the project
echo "Compiling project with 'make'..."
make
if [ $? -ne 0 ]; then
    echo "Make failed! Fix compilation errors before running benchmarks."
    exit 1
fi
echo "Compile successful."

# 2. Check if the executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    exit 1
fi

# 3. Create the log directory and the CSV header
mkdir -p $LOG_DIR
echo "Config,ArrayType,Method,Policy,RealTime_s,UserTime_s,SysTime_s" > $OUTPUT_FILE

# 4. Define all the jobs you want to run
run_all_benchmarks() {
    # --- Add your desired benchmark runs here ---

    # Note: Your LocAG.cpp currently has a hardcoded loop for t=2 and lambda=1-4.
    # Each 'run_benchmark' call will execute that entire loop.

    # Run 'SPINS' config with all 4 algorithm/policy combinations
    run_benchmark "SPINS" "locating" "ga" "parallel"
    run_benchmark "SPINS" "locating" "density" "serial"

    # Run 'Mobile' config with all 4 algorithm/policy combinations
    run_benchmark "Mobile" "locating" "ga" "parallel"
    run_benchmark "Mobile" "locating" "density" "serial"

    # Run a 'detecting' array example
    run_benchmark "2^30" "detecting" "ga" "serial"
    run_benchmark "2^30" "detecting" "ga" "parallel"

    # --- Add more configs as needed ---
    # run_benchmark "Flex" "locating" "ga" "serial"
    # run_benchmark "TCAS" "locating" "ga" "serial"
    
}

# 5. Run all benchmarks
echo "Starting all benchmarks... Results will be saved to $OUTPUT_FILE"
START_TIME=$(date +%s)

run_all_benchmarks

END_TIME=$(date +%s)
TOTAL_DURATION=$((END_TIME - START_TIME))
echo "All benchmarks finished in $TOTAL_DURATION seconds."
