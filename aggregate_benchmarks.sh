#!/bin/bash

# --- Configuration ---
EXECUTABLE="./LocAG"
OUTPUT_FILE="benchmark_results.csv"
LOG_DIR="benchmark_logs"

# Timeout for the ENTIRE execution (runs lambda 1-4 sequentially)
# If this limit is hit, the program is killed, and subsequent lambdas in the loop won't run.
TIMEOUT_DURATION="45m"

# How many times to repeat the full experiment for statistics
TRIALS=5

# List of configurations defined in your LocAG.cpp
CONFIGS=(
    "SPINS"
    "Mobile"
    "Apache"
    "Bugzilla"
    "Flex"
    "GCC"
    "Make"
    "SPINV"
    "TCAS"
    "Wireless"
    "2^30"
)
# --- End of Configuration ---

# 1. Compile
echo "Compiling project..."
make
if [ $? -ne 0 ]; then
    echo "Make failed! Fix compilation errors before running benchmarks."
    exit 1
fi

# 2. Prepare Output
mkdir -p $LOG_DIR
# CSV Header
echo "Config,ArrayType,Method,Policy,d,t,lambda,Trial,Time_ms,Status" > $OUTPUT_FILE

# Function to parse log and extract specific lambda timings
parse_and_log() {
    local log_file=$1
    local config=$2
    local type=$3
    local method=$4
    local policy=$5
    local trial=$6
    local status=$7

    # awk script to find completed runs in the log
    awk -v config="$config" -v type="$type" -v method="$method" -v policy="$policy" -v trial="$trial" -v status="$status" '
    /------------d=/ {
        # Reset vars for new block
        d_val=""; t_val=""; l_val="";

        # Parse line like: ------------d=1, t=2, lambda=1, filename=...
        gsub(/-/, "", $0);
        split($0, parts, ",");
        for (i in parts) {
            split(parts[i], kv, "=");
            key = kv[1]; sub(/^ /, "", key);
            val = kv[2]; sub(/^ /, "", val);
            if (key == "d") d_val = val;
            if (key == "t") t_val = val;
            if (key == "lambda") l_val = val;
        }
    }
    /Solution 0:/ {
        # If we found a solution line, this specific lambda finished successfully
        match($0, /Time total=[0-9]+/);
        t_str = substr($0, RSTART, RLENGTH);
        split(t_str, t_parts, "=");
        time_ms = t_parts[2];

        # Log this specific lambda result as COMPLETED
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", config, type, method, policy, d_val, t_val, l_val, trial, time_ms, "COMPLETED";
    }
    ' "$log_file" >> "$OUTPUT_FILE"
}

# Function to run a suite
run_suite() {
    local config=$1
    local type=$2
    local method=$3
    local policy=$4

    for (( i=1; i<=TRIALS; i++ ))
    do
        local log_file="${LOG_DIR}/${config}_${type}_${method}_${policy}_run${i}.log"
        echo "  [Run $i/$TRIALS] $config ($method/$policy)..."

        # --- EXECUTION WITH TIMEOUT ---
        # "timeout" runs the command and kills it if it exceeds duration.
        # It returns exit code 124 if it timed out.
        timeout "$TIMEOUT_DURATION" $EXECUTABLE "$config" "$type" "$method" "$policy" > "$log_file"
        local exit_code=$?

        local run_status="COMPLETED"

        if [ $exit_code -eq 124 ]; then
            echo "    WARNING: Run timed out after $TIMEOUT_DURATION!"
            run_status="TIMEOUT"
        elif [ $exit_code -ne 0 ]; then
            echo "    ERROR: Run failed with exit code $exit_code. See $log_file"
            run_status="ERROR"
        fi

        # Parse whatever output made it to the log file before the kill/finish
        parse_and_log "$log_file" "$config" "$type" "$method" "$policy" "$i" "$run_status"

        # Optional: If you want to record the Timeout event itself for the *whole* batch in CSV,
        # you could add a specific echo here, but usually knowing which lambdas finished is enough.
    done
}

# --- Main Execution Loop ---

echo "Starting benchmarks with $TIMEOUT_DURATION timeout..."
START_TIME=$(date +%s)

for conf in "${CONFIGS[@]}"; do
    echo "Processing Configuration: $conf"

    # 1. GA Parallel
    run_suite "$conf" "locating" "ga" "parallel"

    # 2. Greedy Serial
    run_suite "$conf" "locating" "greedy" "serial"

    # 3. CE Serial
    run_suite "$conf" "locating" "ce" "serial"

done

END_TIME=$(date +%s)
TOTAL_DURATION=$((END_TIME - START_TIME))
echo "All benchmarks finished in $TOTAL_DURATION seconds."
