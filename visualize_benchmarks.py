import pandas as pd
import matplotlib.pyplot as plt
import sys

# --- Configuration ---
INPUT_FILE = "benchmark_results.csv"
OUTPUT_LOCATING = "benchmark_time_locating.png"
OUTPUT_DETECTING = "benchmark_time_detecting.png"
# ---

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

    # Pivot
    pivot_mean = final_stats.pivot(index='Config', columns='RunType', values='mean').fillna(0)
    pivot_std = final_stats.pivot(index='Config', columns='RunType', values='std').fillna(0)

    # Plot
    ax = pivot_mean.plot(
        kind='bar',
        figsize=(14, 8),
        width=0.8,
        yerr=pivot_std,
        capsize=4,
        edgecolor="black",
        rot=0
    )

    ax.set_title(f'{title} Performance (Avg Total Time)', fontsize=16, pad=20)
    ax.set_ylabel('Time (seconds)', fontsize=12)
    ax.set_xlabel('Configuration', fontsize=12)
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)

    # Add labels
    for container in ax.containers:
        if isinstance(container[0], plt.Rectangle):
            labels = [f'{v:.2f}' if v > 0 else "T/O" for v in container.datavalues]
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

    # 1. Clean Data
    if 'Status' in data.columns:
        data = data[data['Status'] == 'COMPLETED']

    # Drop duplicates just in case
    data = data.drop_duplicates(subset=['Config', 'ArrayType', 'Method', 'Policy', 'd', 't', 'lambda', 'Trial'])

    data['Time_s'] = pd.to_numeric(data['Time_ms'], errors='coerce') / 1000.0
    data = data.dropna(subset=['Time_s'])

    # 2. Split and Plot
    locating_df = data[data['ArrayType'] == 'locating']
    detecting_df = data[data['ArrayType'] == 'detecting']

    plot_data(locating_df, "Locating Array", OUTPUT_LOCATING)
    plot_data(detecting_df, "Detecting Array", OUTPUT_DETECTING)

if __name__ == "__main__":
    main()
