import pandas as pd
import glob
import re

# Find all files matching the pattern.
files = glob.glob("tsv_output/psi11_*_neisize_*_*.tsv")

# Dictionary to store DataFrames keyed by (psi11, neisize)
data_by_param = {}

# Regex to extract psi11 value (float), neisize value (int), and run number (int)
pattern = r'psi11_([0-9.]+)_neisize_([0-9]+)_opt_([0-9.]+)_([0-9]+)\.tsv'

for file in files:
    match = re.search(pattern, file)
    if match:
        psi_val = float(match.group(1))
        neisize_val = int(match.group(2))
        opt_val = float(match.group(3))  # This is the optimization value, not used in key
        key = (psi_val, neisize_val, opt_val)
        # Read the TSV file into a DataFrame
        df = pd.read_csv(file, sep='\t')
        data_by_param.setdefault(key, []).append(df)

# Compute the average DataFrame for each parameter set
averaged_data = {}
for key, dfs in data_by_param.items():
    # Concatenate DataFrames vertically (assuming same row order)
    combined_df = pd.concat(dfs)
    # Group by the row index to calculate the mean
    avg_df = combined_df.groupby(combined_df.index).mean()
    averaged_data[key] = avg_df

# Write each averaged DataFrame into a separate TSV file
for key, avg_df in averaged_data.items():
    psi_val, neisize_val, opt_val = key
    # Define a filename, e.g., psi11_0.5_neisize_10_averaged.tsv
    output_filename = f"tsv_output/psi11_{psi_val}_neisize_{neisize_val}_opt_{opt_val}_averaged.tsv"
    avg_df.to_csv(output_filename, sep="\t", index=False)
    print(f"Output written to {output_filename}")
