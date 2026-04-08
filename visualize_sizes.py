import os
import re
import sys
import pandas as pd
import matplotlib.pyplot as plt
import math

# --- Configuration ---
# Change this to "benchmark_logs_subset" if you are running the subset
LOG_DIR = "benchmark_logs"
OUTPUT_PREFIX = "benchmark_sizes_comparison"

# --- FIXED COLOR PALETTE ---
COLOR_MAP = {
    'ga-parallel': '#1f77b4',    # Blue
    'greedy-serial': '#ff7f0e',  # Orange
    'ce-serial': '#2ca02c',      # Green
    'ga-serial': '#9467bd',      # Purple
    'greedy-parallel': '#8c564b',# Brown
    'ce-parallel': '#e377c2'     # Pink
}
# ---------------------

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
        
        try:
            with open(os.path.join(log_dir, filename), 'r') as f:
                content = f.read()
        except Exception as e:
            print(f"Skipping {filename}: {e}")
            continue

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

def generate_comparison_chart(df, configs, filename_suffix):
    if not configs:
        return

    num_rows = len(configs)
    fig, axes = plt.subplots(nrows=num_rows, ncols=2, figsize=(14, 5 * num_rows), squeeze=False)

    for i, config in enumerate(configs):
        config_data = df[df['Config'] == config]
        
        # --- LEFT PLOT: Locating ---
        ax_loc = axes[i][0]
        loc_data = config_data[config_data['ArrayType'] == 'locating']
        
        if not loc_data.empty:
            pivot_loc = loc_data.pivot_table(index='Lambda', columns='RunType', values='N_total')
            
            # Map columns to fixed colors
            colors_loc = [COLOR_MAP.get(col, '#7f7f7f') for col in pivot_loc.columns]
            
            pivot_loc.plot(kind='bar', ax=ax_loc, width=0.8, edgecolor='black', rot=0, color=colors_loc)
            ax_loc.set_title(f"{config} - Locating", fontsize=14, fontweight='bold')
            ax_loc.set_ylabel("Array Size (N)")
            ax_loc.grid(axis='y', linestyle='--', alpha=0.7)
            ax_loc.legend(loc='upper left', fontsize='small')
            
            for container in ax_loc.containers:
                ax_loc.bar_label(container, fmt='%d', padding=3, fontsize=8)
        else:
            ax_loc.text(0.5, 0.5, "No Locating Data", ha='center', va='center')
            ax_loc.set_title(f"{config} - Locating")

        # --- RIGHT PLOT: Detecting ---
        ax_det = axes[i][1]
        det_data = config_data[config_data['ArrayType'] == 'detecting']
        
        if not det_data.empty:
            pivot_det = det_data.pivot_table(index='Lambda', columns='RunType', values='N_total')
            
            # Map columns to fixed colors
            colors_det = [COLOR_MAP.get(col, '#7f7f7f') for col in pivot_det.columns]
            
            pivot_det.plot(kind='bar', ax=ax_det, width=0.8, edgecolor='black', rot=0, color=colors_det)
            ax_det.set_title(f"{config} - Detecting", fontsize=14, fontweight='bold')
            ax_det.set_ylabel("") 
            ax_det.grid(axis='y', linestyle='--', alpha=0.7)
            ax_det.legend(loc='upper left', fontsize='small')

            for container in ax_det.containers:
                ax_det.bar_label(container, fmt='%d', padding=3, fontsize=8)
        else:
            ax_det.text(0.5, 0.5, "No Detecting Data\n(Timeout?)", ha='center', va='center')
            ax_det.set_title(f"{config} - Detecting")

    plt.tight_layout(rect=[0, 0.03, 1, 0.98]) 
    output_file = f"{OUTPUT_PREFIX}_{filename_suffix}.png"
    plt.savefig(output_file)
    print(f"Saved comparison chart to '{output_file}'")
    plt.close()

if __name__ == "__main__":
    data = parse_log_files(LOG_DIR)
    if not data:
        print("No data parsed.")
        sys.exit(0)
        
    df = pd.DataFrame(data)
    
    # Aggregate Trials (Mean N)
    df_agg = df.groupby(['Config', 'ArrayType', 'Lambda', 'RunType'])['N_total'].mean().reset_index()

    # Get unique configs and sort them
    all_configs = sorted(df_agg['Config'].unique())
    
    # Split into chunks
    CHUNK_SIZE = 6
    total_chunks = math.ceil(len(all_configs) / CHUNK_SIZE)
    
    for chunk_idx in range(total_chunks):
        start_i = chunk_idx * CHUNK_SIZE
        end_i = start_i + CHUNK_SIZE
        config_chunk = all_configs[start_i:end_i]
        
        generate_comparison_chart(df_agg, config_chunk, f"part{chunk_idx+1}")