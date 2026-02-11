import pandas as pd
import sys

# --- Configuration ---
INPUT_FILE = "benchmark_results.csv"
OUTPUT_FILE = "benchmark_table.tex"
# ---------------------

def generate_latex_table(dataframe, title, label):
    if dataframe.empty:
        return f"% No data for {title}\n"

    pivot = dataframe.pivot(index='Config', columns='RunType', values='Result')
    pivot = pivot.fillna("T/O") # Timeout

    latex = pivot.to_latex(
        escape=False,
        caption=f"{title} Execution Times (seconds)",
        label=label,
        position="htbp",
        column_format="l" + "c" * len(pivot.columns)
    )
    return latex + "\n\\vspace{1em}\n"

def generate_report(input_file, output_file):
    try:
        data = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"Error: {input_file} not found.")
        sys.exit(1)

    if 'Status' in data.columns:
        data = data[data['Status'] == 'COMPLETED']

    # --- SAFETY: Remove duplicates from appended runs ---
    # We identify a unique run by Config, Type, Method, Policy, Trial, and specific Lambda
    data = data.drop_duplicates(subset=['Config', 'ArrayType', 'Method', 'Policy', 'd', 't', 'lambda', 'Trial'])

    data['Time_s'] = pd.to_numeric(data['Time_ms'], errors='coerce') / 1000.0
    data = data.dropna(subset=['Time_s'])

    # Sum lambdas per trial
    trial_totals = data.groupby(
        ['Config', 'ArrayType', 'Method', 'Policy', 'Trial']
    )['Time_s'].sum().reset_index()

    # Average trials
    stats = trial_totals.groupby(
        ['Config', 'ArrayType', 'Method', 'Policy']
    )['Time_s'].agg(['mean', 'std']).reset_index()

    stats['std'] = stats['std'].fillna(0)
    stats['Result'] = stats.apply(lambda r: f"{r['mean']:.2f} $\\pm$ {r['std']:.2f}", axis=1)
    stats['RunType'] = stats['Method'] + ' (' + stats['Policy'] + ')'

    with open(output_file, "w") as f:
        f.write("% --- Locating Table ---\n")
        f.write(generate_latex_table(stats[stats['ArrayType'] == 'locating'], "Locating Array", "tab:locating"))
        f.write("\n% --- Detecting Table ---\n")
        f.write(generate_latex_table(stats[stats['ArrayType'] == 'detecting'], "Detecting Array", "tab:detecting"))

    print(f"Tables saved to '{output_file}'")

if __name__ == "__main__":
    generate_report(INPUT_FILE, OUTPUT_FILE)
