import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt

# --- Configuration ---
LOG_DIR = "benchmark_logs"
OUTPUT_IMAGE = "benchmark_array_sizes.png"
# ---

def parse_log_files(log_dir):
    """
    Parses all .log files in the given directory to extract benchmark results.
    """
    parsed_data = []
    
    if not os.path.exists(log_dir):
        print(f"Error: Log directory '{log_dir}' not found.")
        print("Please run './run_benchmarks.sh' first.")
        sys.exit(1)

    # Regex to parse the filename, e.g., "SPINS_locating_ga_serial.log"
    # Groups: 1=Config, 2=Type, 3=Method, 4=Policy
    filename_re = re.compile(r'(.+?)_(locating|detecting)_(ga|density)_(serial|parallel)\.log')
    
    # Regex to find lambda and the *first* "Solution 0" N total for that lambda
    # re.DOTALL makes '.' match newlines, so we can span across multiple lines.
    log_content_re = re.compile(r"lambda=(\d+).*?Solution 0: N total=(\d+)", re.DOTALL)

    for filename in os.listdir(log_dir):
        if not filename.endswith(".log"):
            continue

        match = filename_re.match(filename)
        if not match:
            print(f"Skipping unrecognized log file: {filename}")
            continue
            
        config, array_type, method, policy = match.groups()
        run_type = f"{method}-{policy}"
        
        log_path = os.path.join(log_dir, filename)
        with open(log_path, 'r') as f:
            content = f.read()

        # Find all (lambda, N_total) pairs in this log file
        results = log_content_re.findall(content)
        
        if not results:
            print(f"Warning: No 'Solution 0' found in {filename}.")
            continue

        for lambda_val, n_total in results:
            parsed_data.append({
                'Config': config,
                'ArrayType': array_type,
                'RunType': run_type,
                'Lambda': int(lambda_val),
                'N_total': int(n_total)
            })

    return parsed_data

def visualize_array_sizes(data, output_image):
    """
    Generates and saves a grouped bar chart for array sizes.
    """
    if not data:
        print("No valid data parsed from log files. Nothing to plot.")
        sys.exit(0)

    df = pd.DataFrame(data)

    # Get a list of unique configurations to create one plot for each
    configs = df['Config'].unique()
    num_configs = len(configs)
    
    if num_configs == 0:
        print("No configurations found in parsed data.")
        return

    # Create a subplot for each configuration
    fig, axes = plt.subplots(nrows=num_configs, ncols=1, 
                             figsize=(15, 7 * num_configs), squeeze=False)
    
    for i, config in enumerate(configs):
        ax = axes[i][0]
        config_data = df[df['Config'] == config]

        # Pivot data to get Lambda on x-axis and RunType as bar groups
        try:
            pivot = config_data.pivot_table(index='Lambda', columns='RunType', values='N_total', aggfunc='first')
        except Exception as e:
            print(f"Error pivoting data for config '{config}': {e}")
            continue

        pivot.plot(kind='bar', ax=ax, width=0.8, edgecolor='black')

        # --- Prettify the Plot ---
        ax.set_title(f'Final Array Size (N total) for "{config}"', fontsize=18, pad=20)
        ax.set_ylabel('Total Rows (N)', fontsize=12)
        ax.set_xlabel('Lambda Value', fontsize=12)
        ax.tick_params(axis='x', rotation=0)
        ax.yaxis.grid(True, linestyle='--', alpha=0.7)
        ax.set_axisbelow(True)
        ax.legend(title='Algorithm (Method-Policy)', loc='upper left')

        # Add data labels on top of each bar
        for bar in ax.patches:
            height = bar.get_height()
            if pd.notna(height) and height > 0:
                ax.annotate(
                    f'{int(height)}',  # Format as integer
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center',
                    va='bottom',
                    fontsize=9
                )

    plt.tight_layout(pad=3.0)  # Adjust plot to prevent labels from overlapping

    # --- 5. Save and Show ---
    plt.savefig(output_image)
    print(f"Array size visualization saved to '{output_image}'")

if __name__ == "__main__":
    parsed_data = parse_log_files(LOG_DIR)
    visualize_array_sizes(parsed_data, OUTPUT_IMAGE)