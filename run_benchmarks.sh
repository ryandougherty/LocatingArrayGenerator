#!/bin/bash

# --- Configuration ---
EXECUTABLE="./LocAG"
OUTPUT_FILE="benchmark_results.csv"
LOG_DIR="benchmark_logs"
TIMEOUT_DURATION="45m"
TRIALS=5

# Configurations to test
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

echo "Compiling project..."
make
if [ $? -ne 0 ]; then
    echo "Make failed!"
    exit 1
fi

mkdir -p $LOG_DIR

# --- SAFETY CHECK: Handle Output File ---
if [ -f "$OUTPUT_FILE" ]; then
    echo "Found existing $OUTPUT_FILE. Backing it up to ${OUTPUT_FILE}.bak..."
    cp "$OUTPUT_FILE" "${OUTPUT_FILE}.bak"
    echo "Appending new results to existing $OUTPUT_FILE..."
else
    # Only write the header if the file usually doesn't exist
    echo "Config,ArrayType,Method,Policy,d,t,lambda,Trial,Time_ms,Status" > $OUTPUT_FILE
fi

parse_and_log() {
    local log_file=$1
    local config=$2
    local type=$3
    local method=$4
    local policy=$5
    local trial=$6
    local status=$7

    awk -v config="$config" -v type="$type" -v method="$method" -v policy="$policy" -v trial="$trial" -v status="$status" '
    /------------d=/ {
        d_val=""; t_val=""; l_val="";
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
        match($0, /Time total=[0-9]+/);
        t_str = substr($0, RSTART, RLENGTH);
        split(t_str, t_parts, "=");
        time_ms = t_parts[2];
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", config, type, method, policy, d_val, t_val, l_val, trial, time_ms, "COMPLETED";
    }
    ' "$log_file" >> "$OUTPUT_FILE"
}

run_suite() {
    local config=$1
    local type=$2
    local method=$3
    local policy=$4

    for (( i=1; i<=TRIALS; i++ ))
    do
        local log_file="${LOG_DIR}/${config}_${type}_${method}_${policy}_run${i}.log"

        # Skip if log file exists and is non-empty (Optional: uncomment to resume interrupted runs)
        # if [ -s "$log_file" ]; then echo "Skipping completed run $i..."; continue; fi

        echo "  [Run $i/$TRIALS] $config $type ($method/$policy)..."

        timeout "$TIMEOUT_DURATION" $EXECUTABLE "$config" "$type" "$method" "$policy" > "$log_file"
        local exit_code=$?
        local run_status="COMPLETED"

        if [ $exit_code -eq 124 ]; then
            echo "    WARNING: Run timed out!"
            run_status="TIMEOUT"
        elif [ $exit_code -ne 0 ]; then
            echo "    ERROR: Run failed (exit code $exit_code)."
            run_status="ERROR"
        fi

        parse_and_log "$log_file" "$config" "$type" "$method" "$policy" "$i" "$run_status"
    done
}

echo "Starting benchmarks with $TIMEOUT_DURATION timeout..."
START_TIME=$(date +%s)

for conf in "${CONFIGS[@]}"; do
    echo "=== Configuration: $conf ==="

    # --- LOCATING ARRAYS (Commented out to preserve your 24h run) ---
    echo "Skipping Locating Arrays (already done)..."
    run_suite "$conf" "locating" "ga" "parallel"
    run_suite "$conf" "locating" "greedy" "serial"
    run_suite "$conf" "locating" "ce" "serial"

    # --- DETECTING ARRAYS (Enabled) ---
    run_suite "$conf" "detecting" "ga" "parallel"
    run_suite "$conf" "detecting" "greedy" "serial"
    run_suite "$conf" "detecting" "ce" "serial"
done

END_TIME=$(date +%s)
TOTAL_DURATION=$((END_TIME - START_TIME))
echo "All benchmarks finished in $TOTAL_DURATION seconds."
