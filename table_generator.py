import pandas as pd
import sys

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

    # 2. Pre-process Data
    # Convert 'RealTime_s' to numeric, forcing errors to NaN
    data['RealTime_s'] = pd.to_numeric(data['RealTime_s'], errors='coerce')
    
    # Create the column header for the table (Method + Policy)
    # e.g., "ga (parallel)", "ce (serial)"
    data['RunType'] = data['Method'] + ' (' + data['Policy'] + ')'

    # 3. Pivot the Data
    # Rows: Config + ArrayType
    # Columns: RunType
    # Values: RealTime_s
    try:
        pivot_data = data.pivot(
            index=['Config', 'ArrayType'], 
            columns='RunType', 
            values='RealTime_s'
        )
    except ValueError:
        print("Error: Duplicate entries found. Check if benchmark_results.csv has duplicates.")
        sys.exit(1)

    # 4. Generate LaTeX
    # format_float: ensures 2 decimal places
    latex_code = pivot_data.to_latex(
        float_format="%.2f",
        na_rep="-",    # Symbol for missing data
        caption="Benchmark Execution Times (seconds)",
        label="tab:benchmark_results",
        position="htbp"
    )

    # 5. Save to file
    with open(output_file, "w") as f:
        f.write(latex_code)
    
    print(f"LaTeX table saved to '{output_file}'")
    print("-" * 30)
    print(latex_code) # Also print to console for quick copy-paste

if __name__ == "__main__":
    generate_latex(INPUT_FILE, OUTPUT_FILE)