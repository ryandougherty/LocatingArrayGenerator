import pandas as pd
import sys
import numpy as np

# --- Configuration ---
INPUT_FILE = "benchmark_results.csv"
OUTPUT_FILE = "benchmark_table.tex"
# ---------------------

def generate_latex(input_file, output_file):
    # 1. Load Data
    try:
        data = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: '{input_file}' not found.")
        sys.exit(1)

    if data.empty:
        print("Error: CSV file is empty.")
        sys.exit(1)

    # 2. Filter and Convert
    # Only use successfully completed runs
    if 'Status' in data.columns:
        data = data[data['Status'] == 'COMPLETED']

    # Convert ms to seconds
    data['Time_s'] = pd.to_numeric(data['Time_ms'], errors='coerce') / 1000.0
    data = data.dropna(subset=['Time_s'])

    # 3. Aggregation
    # Step A: Sum times for all lambdas within a single Trial
    # (Config + RunType + Trial) -> Total Time for that Trial
    trial_totals = data.groupby(
        ['Config', 'ArrayType', 'Method', 'Policy', 'Trial']
    )['Time_s'].sum().reset_index()

    # Step B: Average across Trials and calculate Std Dev
    stats = trial_totals.groupby(
        ['Config', 'ArrayType', 'Method', 'Policy']
    )['Time_s'].agg(['mean', 'std']).reset_index()

    # Handle NaN std dev (if only 1 trial exists)
    stats['std'] = stats['std'].fillna(0)

    # 4. Format for LaTeX
    # Create a string "Mean ± Std"
    # We use a helper function to format cleanly
    def format_cell(row):
        return f"{row['mean']:.2f} $\\pm$ {row['std']:.2f}"

    stats['Result'] = stats.apply(format_cell, axis=1)

    # Create readable column headers
    stats['RunType'] = stats['Method'] + ' (' + stats['Policy'] + ')'

    # 5. Pivot
    # Rows: Config
    # Columns: RunType
    # Values: Result string
    pivot_table = stats.pivot(
        index=['Config', 'ArrayType'],
        columns='RunType',
        values='Result'
    )

    # 6. Generate LaTeX
    latex_code = pivot_table.to_latex(
        escape=False, # Allow the $\pm$ math symbols to render
        caption="Benchmark Execution Times (seconds, Mean $\\pm$ Std Dev)",
        label="tab:benchmark_results",
        position="htbp",
        column_format="l" * (len(pivot_table.columns) + 1) # simple left alignment
    )

    # 7. Save
    with open(output_file, "w") as f:
        f.write(latex_code)

    print(f"LaTeX table saved to '{output_file}'")
    print("-" * 30)
    print(latex_code)

if __name__ == "__main__":
    generate_latex(INPUT_FILE, OUTPUT_FILE)
