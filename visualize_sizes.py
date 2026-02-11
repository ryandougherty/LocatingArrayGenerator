import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt

# --- Configuration ---
LOG_DIR = "benchmark_logs"
OUTPUT_LOCATING = "benchmark_size_locating.png"
OUTPUT_DETECTING = "benchmark_size_detecting.png"
# ---

def parse_log_files(log_dir):
    parsed_data = []
    if not os.path.exists(log_dir):
        print(f"Error: {log_dir} not found.")
        sys.exit(1)

    filename_re = re.compile(r'(.+?)_(locating|detecting)_(ga|ce|greedy)_(serial|parallel)_run(\d+)\.log')
    log_content_re = re.compile(r"lambda=(\d+).*?Solution 0: N total=(\d+)", re.DOTALL)

    for filename in os.listdir(log_dir):
        if not filename.endswith(".log"): continue
        match = filename_re.match(filename)
        if not match: continue

        config, array_type, method, policy, trial = match.groups()

        with open(os.path.join(log_dir, filename), 'r') as f:
            content = f.read()

        results = log_content_re.findall(content)
        for lambda_val, n_total in results:
            parsed_data.append({
                'Config': config,
                'ArrayType': array_type,
                'RunType': f"{method}-{policy}",
                'Trial': int(trial),
                'Lambda': int(lambda_val),
                'N_total': int(n_total)
            })
    return parsed_data

def plot_sizes(df, title, output_file):
    if df.empty: return

    # Average N across trials
    df_agg = df.groupby(['Config', 'Lambda', 'RunType'])['N_total'].mean().reset_index()
    configs = df_agg['Config'].unique()

    if len(configs) == 0: return

    fig, axes = plt.subplots(nrows=len(configs), ncols=1, figsize=(12, 6 * len(configs)), squeeze=False)

    for i, config in enumerate(configs):
        ax = axes[i][0]
        subset = df_agg[df_agg['Config'] == config]

        pivot = subset.pivot(index='Lambda', columns='RunType', values='N_total')
        pivot.plot(kind='bar', ax=ax, width=0.8, edgecolor='black', rot=0)

        ax.set_title(f'{title}: {config}', fontsize=14)
        ax.set_ylabel('Array Size (N)')
        ax.set_xlabel('Lambda')
        ax.grid(axis='y', linestyle='--', alpha=0.7)

        for container in ax.containers:
             ax.bar_label(container, fmt='%d', padding=3, fontsize=9)

    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Saved {output_file}")
    plt.close()

if __name__ == "__main__":
    data = parse_log_files(LOG_DIR)
    if not data:
        print("No data parsed.")
        sys.exit(0)

    df = pd.DataFrame(data)

    plot_sizes(df[df['ArrayType'] == 'locating'], "Locating Sizes", OUTPUT_LOCATING)
    plot_sizes(df[df['ArrayType'] == 'detecting'], "Detecting Sizes", OUTPUT_DETECTING)
