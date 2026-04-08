import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np

# --- Configuration ---
INPUT_FILE = "benchmark_results.csv" # Change to subset if needed
OUTPUT_LOCATING = "benchmark_time_locating_log.png"
OUTPUT_DETECTING = "benchmark_time_detecting_log.png"

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

def plot_data(df, title, output_file):
    if df.empty:
        print(f"No data found for {title}. Skipping plot.")
        return

    # Aggregate: Sum lambdas per Trial -> Average across Trials
    trial_totals = df.groupby(
        ['Config', 'Method', 'Policy', 'Trial']
    )['Time_s'].sum().reset_index()

    final_stats = trial_totals.groupby(
        ['Config', 'Method', 'Policy']
    )['Time_s'].agg(['mean', 'std']).reset_index()

    final_stats['std'] = final_stats['std'].fillna(0)
    final_stats['RunType'] = final_stats['Method'] + '-' + final_stats['Policy']

    # Pivot (Notice we DO NOT fill NaNs with 0 here anymore, because log(0) breaks the chart)
    pivot_mean = final_stats.pivot(index='Config', columns='RunType', values='mean')
    pivot_std = final_stats.pivot(index='Config', columns='RunType', values='std')

    # Map the columns to our fixed colors
    bar_colors = [COLOR_MAP.get(col, '#7f7f7f') for col in pivot_mean.columns]

    # Plot
    ax = pivot_mean.plot(
        kind='bar',
        figsize=(15, 8),
        width=0.8,
        yerr=pivot_std,
        capsize=4,
        edgecolor="black",
        rot=0,
        color=bar_colors
    )

    # --- THE MAGIC: Set Y-Axis to Logarithmic Scale ---
    ax.set_yscale('log')

    ax.set_title(f'{title} Performance (Avg Total Time - Log Scale)', fontsize=16, pad=20)
    ax.set_ylabel('Time (seconds) [Log10 Scale]', fontsize=12)
    ax.set_xlabel('Configuration', fontsize=12)
    
    # Add minor grid lines which look great on log charts
    ax.yaxis.grid(True, linestyle='-', alpha=0.7, which='major')
    ax.yaxis.grid(True, linestyle='--', alpha=0.3, which='minor')
    ax.set_axisbelow(True)
    
    # Add labels to the top of the bars
    for container in ax.containers:
        if isinstance(container[0], plt.Rectangle):
            labels = []
            for v in container.datavalues:
                if pd.isna(v) or v <= 0:
                    labels.append("T/O") # Timeout / No Data
                else:
                    # Format standard numbers cleanly
                    if v < 10:
                        labels.append(f'{v:.2f}')
                    else:
                        labels.append(f'{int(v)}')
            ax.bar_label(container, labels=labels, padding=3, fontsize=8)

    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Saved {output_file}")
    plt.close()

def main():
    try:
        data = pd.read_csv(INPUT_FILE)
    except Exception as e:
        print(f"Error reading {INPUT_FILE}: {e}")
        sys.exit(1)

    # Clean Data
    if 'Status' in data.columns:
        data = data[data['Status'] == 'COMPLETED']
    
    data = data.drop_duplicates(subset=['Config', 'ArrayType', 'Method', 'Policy', 'd', 't', 'lambda', 'Trial'])
    data['Time_s'] = pd.to_numeric(data['Time_ms'], errors='coerce') / 1000.0
    data = data.dropna(subset=['Time_s'])

    # Split by Array Type (No longer splitting by outliers!)
    locating_df = data[data['ArrayType'] == 'locating']
    detecting_df = data[data['ArrayType'] == 'detecting']

    # Plot everything onto just two charts
    plot_data(locating_df, "Locating Array", OUTPUT_LOCATING)
    plot_data(detecting_df, "Detecting Array", OUTPUT_DETECTING)

if __name__ == "__main__":
    main()