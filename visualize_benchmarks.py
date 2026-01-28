import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# --- Configuration ---
INPUT_FILE = "benchmark_results.csv"
OUTPUT_IMAGE = "benchmark_results.png"
# ---

def visualize_results(input_file, output_image):
    """
    Reads benchmark data from a CSV, aggregates trials, and generates
    a grouped bar chart with error bars.
    """
    # --- 1. Load Data ---
    try:
        data = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{input_file}' is empty.")
        sys.exit(1)

    # Filter for completed runs only
    if 'Status' in data.columns:
        data = data[data['Status'] == 'COMPLETED']

    # Convert Time_ms to seconds
    data['Time_s'] = pd.to_numeric(data['Time_ms'], errors='coerce') / 1000.0
    data = data.dropna(subset=['Time_s'])

    if data.empty:
        print("No valid benchmark data found to plot.")
        sys.exit(0)

    # --- 2. Aggregation Strategy ---
    # Goal: Calculate "Total Time" for the suite (sum of all lambdas) per Trial.
    # Then average those totals across Trials.

    # Group 1: Sum times across all lambdas for each unique Trial
    # (Config + Method + Policy + Trial -> Total Time)
    trial_totals = data.groupby(
        ['Config', 'ArrayType', 'Method', 'Policy', 'Trial']
    )['Time_s'].sum().reset_index()

    # Group 2: Calculate Mean and Std Dev across Trials
    # (Config + Method + Policy -> Mean Time, Std Dev)
    final_stats = trial_totals.groupby(
        ['Config', 'ArrayType', 'Method', 'Policy']
    )['Time_s'].agg(['mean', 'std']).reset_index()

    # Fill NaN std devs (happens if only 1 trial) with 0
    final_stats['std'] = final_stats['std'].fillna(0)

    # Create a label for the bars
    final_stats['RunType'] = final_stats['Method'] + '-' + final_stats['Policy']

    # Pivot for plotting
    # Rows = Config, Columns = RunType, Values = Mean Time
    pivot_mean = final_stats.pivot(index='Config', columns='RunType', values='mean')
    pivot_std = final_stats.pivot(index='Config', columns='RunType', values='std')

    # --- 3. Create the Plot ---
    ax = pivot_mean.plot(
        kind='bar',
        figsize=(16, 9),
        width=0.8,
        yerr=pivot_std,  # Add error bars
        capsize=4,       # Caps on error bars
        edgecolor="black",
        rot=0            # Horizontal x-labels
    )

    # --- 4. Prettify the Plot ---
    ax.set_title(f'LocAG Benchmark Performance (Avg Total Time over Trials)', fontsize=18, pad=20)
    ax.set_ylabel('Wall-Clock Time (seconds)', fontsize=12)
    ax.set_xlabel('Benchmark Configuration', fontsize=12)

    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)

    # Add labels
    for container in ax.containers:
        # Check if this container is for the bars (and not the error bars)
        if isinstance(container[0], plt.Rectangle):
            ax.bar_label(container, fmt='%.2f', padding=3, fontsize=9)

    ax.legend(title='Run Type', loc='upper left')
    plt.tight_layout()

    # --- 5. Save ---
    plt.savefig(output_image)
    print(f"Benchmark visualization saved to '{output_image}'")

if __name__ == "__main__":
    visualize_results(INPUT_FILE, OUTPUT_IMAGE)
