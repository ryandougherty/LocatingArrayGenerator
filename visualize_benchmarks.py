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
    Reads benchmark data from a CSV and generates a grouped bar chart.
    """
    # --- 1. Load Data ---
    try:
        data = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
        print("Please run the './run_benchmarks.sh' script first to generate the data.")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: The file '{input_file}' is empty.")
        print("Please ensure the benchmark script ran correctly.")
        sys.exit(1)

    # Convert time to numeric, replacing 'ERROR' strings with NaN (Not a Number)
    data['RealTime_s'] = pd.to_numeric(data['RealTime_s'], errors='coerce')

    # Drop any rows that had errors
    data = data.dropna(subset=['RealTime_s'])

    if data.empty:
        print("No valid benchmark data found to plot.")
        sys.exit(0)

    # --- 2. Prepare Data for Plotting ---
    
    # Create a combined column for the bar groups (e.g., "ga-serial", "ga-parallel")
    data['RunType'] = data['Method'] + '-' + data['Policy']
    
    # Pivot the data to get Configs as rows and RunTypes as columns
    # This is the ideal format for a grouped bar chart
    try:
        # pivot_data = data.pivot(index='Config', columns='RunType', values='RealTime_s')
        pivot_data = data.pivot(index=['Config', 'ArrayType'], columns='RunType', values='RealTime_s')
    except Exception as e:
        print(f"Error pivoting data: {e}")
        print("Please check your benchmark_results.csv file for correct formatting.")
        sys.exit(1)

    # --- 3. Create the Plot ---
    
    # Use pandas' built-in plotting (which uses matplotlib)
    ax = pivot_data.plot(
        kind='bar',
        figsize=(16, 9),  # Wider figure size for readability
        width=0.8,        # Bar width
        edgecolor="black"
    )

    # --- 4. Prettify the Plot ---
    ax.set_title(f'LocAG Benchmark Performance (Real Time)', fontsize=18, pad=20)
    ax.set_ylabel('Wall-Clock Time (seconds)', fontsize=12)
    ax.set_xlabel('Benchmark Configuration', fontsize=12)
    
    # Make x-axis labels horizontal
    ax.tick_params(axis='x', rotation=0)

    # Add a grid for easier reading
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)

    # Add data labels on top of each bar
    for bar in ax.patches:
        height = bar.get_height()
        if pd.notna(height) and height > 0:
            ax.annotate(
                f'{height:.2f}s',  # Format to 2 decimal places
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points",
                ha='center',
                va='bottom',
                fontsize=9
            )

    # Improve legend
    ax.legend(title='Run Type (Method-Policy)', loc='upper left')

    plt.tight_layout()  # Adjust plot to prevent labels from overlapping

    # --- 5. Save and Show ---
    plt.savefig(output_image)
    print(f"Benchmark visualization saved to '{output_image}'")
    
    # Uncomment the line below if you want the script to open the plot in a window
    # plt.show()

if __name__ == "__main__":
    visualize_results(INPUT_FILE, OUTPUT_IMAGE)