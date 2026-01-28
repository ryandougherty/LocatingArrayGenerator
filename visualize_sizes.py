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
    parsed_data = []

    if not os.path.exists(log_dir):
        print(f"Error: Log directory '{log_dir}' not found.")
        sys.exit(1)

    # UPDATED REGEX: Handles "_run1.log", "_run2.log", etc.
    # Groups: 1=Config, 2=Type, 3=Method, 4=Policy, 5=TrialNumber
    filename_re = re.compile(r'(.+?)_(locating|detecting)_(ga|ce|greedy)_(serial|parallel)_run(\d+)\.log')

    # Regex to find lambda and N total in the logs
    log_content_re = re.compile(r"lambda=(\d+).*?Solution 0: N total=(\d+)", re.DOTALL)

    for filename in os.listdir(log_dir):
        if not filename.endswith(".log"):
            continue

        match = filename_re.match(filename)
        if not match:
            # Try matching the old format just in case
            continue

        config, array_type, method, policy, trial = match.groups()
        run_type = f"{method}-{policy}"

        log_path = os.path.join(log_dir, filename)
        try:
            with open(log_path, 'r') as f:
                content = f.read()
        except Exception as e:
            print(f"Error reading {filename}: {e}")
            continue

        results = log_content_re.findall(content)

        for lambda_val, n_total in results:
            parsed_data.append({
                'Config': config,
                'ArrayType': array_type,
                'RunType': run_type,
                'Trial': int(trial),
                'Lambda': int(lambda_val),
                'N_total': int(n_total)
            })

    return parsed_data

def visualize_array_sizes(data, output_image):
    if not data:
        print("No valid data parsed from log files. Nothing to plot.")
        return

    df = pd.DataFrame(data)

    # --- Aggregation: Average N across Trials ---
    # Since GA is stochastic, size might vary slightly. We take the mean.
    df_agg = df.groupby(['Config', 'Lambda', 'RunType'])['N_total'].mean().reset_index()

    configs = df_agg['Config'].unique()
    num_configs = len(configs)

    if num_configs == 0:
        return

    # Create subplots
    fig, axes = plt.subplots(nrows=num_configs, ncols=1,
                             figsize=(15, 7 * num_configs), squeeze=False)

    for i, config in enumerate(configs):
        ax = axes[i][0]
        config_data = df_agg[df_agg['Config'] == config]

        try:
            pivot = config_data.pivot_table(index='Lambda', columns='RunType', values='N_total')
        except Exception as e:
            continue

        pivot.plot(kind='bar', ax=ax, width=0.8, edgecolor='black', rot=0)

        ax.set_title(f'Average Final Array Size (N) for "{config}"', fontsize=18, pad=20)
        ax.set_ylabel('Total Rows (N)', fontsize=12)
        ax.set_xlabel('Lambda Value', fontsize=12)
        ax.yaxis.grid(True, linestyle='--', alpha=0.7)
        ax.set_axisbelow(True)
        ax.legend(title='Algorithm', loc='upper left')

        # Add labels
        for container in ax.containers:
             ax.bar_label(container, fmt='%d', padding=3, fontsize=9)

    plt.tight_layout(pad=3.0)
    plt.savefig(output_image)
    print(f"Array size visualization saved to '{output_image}'")

if __name__ == "__main__":
    parsed_data = parse_log_files(LOG_DIR)
    visualize_array_sizes(parsed_data, OUTPUT_IMAGE)
