import numpy as np
import matplotlib
matplotlib.use('Agg')  # For backend-safe PDF output
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import argparse
import concurrent.futures
import pandas as pd
from matplotlib import ticker

matplotlib.rcParams['pdf.fonttype'] = 42

# Argument parser for command-line input
parser = argparse.ArgumentParser(description="Compute and plot difference heatmap with multi-threading")
parser.add_argument("--input1", required=True, help="First input matrix file (WT)")
parser.add_argument("--input2", required=True, help="Second input matrix file (KO)")
parser.add_argument("--output", required=True, help="Output heatmap PDF file")
parser.add_argument("--min", type=float, default=-5, help="Minimum value for color scale")
parser.add_argument("--max", type=float, default=5, help="Maximum value for color scale")
parser.add_argument("--dpi", type=int, default=150, help="DPI for output PDF")
parser.add_argument("--cpus", type=int, default=8, help="Number of CPU cores to use")

args = parser.parse_args()

# File paths
wt_matrix_file = args.input1
ko_matrix_file = args.input2
output_pdf_file = args.output
zmin = args.min
zmax = args.max
dpi = args.dpi
num_cpus = args.cpus

# Function to read matrix as pandas Series: index = (chr, start, end, peak_id), value = numpy array
def read_matrix_with_keys(file_path):
    coords = []
    matrix = []
    with gzip.open(file_path, "rt") as f:
        for line in f:
            if line.startswith("@"):
                continue
            fields = line.strip().split()
            key = tuple(fields[0:4])  # chr, start, end, peak_id
            values = np.array(fields[6:], dtype=np.float32)
            coords.append(key)
            matrix.append(values)
    return pd.Series(matrix, index=coords)

# Read both matrices
print("[INFO] Reading matrices...")
wt_series = read_matrix_with_keys(args.input1)
ko_series = read_matrix_with_keys(args.input2)

print(f"[INFO] Raw WT rows: {len(wt_series)}")
print(f"[INFO] Raw KO rows: {len(ko_series)}")

# Align by shared keys, maintaining WT order
common_keys = [key for key in wt_series.index if key in ko_series.index]
print(f"[INFO] Matched rows: {len(common_keys)}")

# Extract aligned matrices
wt_matrix = np.stack(wt_series.loc[common_keys].values)
ko_matrix = np.stack(ko_series.loc[common_keys].values)

print(f"[INFO] Final aligned WT matrix shape: {wt_matrix.shape}")
print(f"[INFO] Final aligned KO matrix shape: {ko_matrix.shape}")

# Parallel computation of KO - WT
def compute_difference(wt_chunk, ko_chunk):
    #return ko_chunk - wt_chunk  # Element-wise subtraction
    return np.log2((ko_chunk + 1e-9) / (wt_chunk + 1e-9))

# Split data into chunks for parallel processing
wt_chunks = np.array_split(wt_matrix, num_cpus)
ko_chunks = np.array_split(ko_matrix, num_cpus)

# Compute difference in parallel
with concurrent.futures.ThreadPoolExecutor(max_workers=num_cpus) as executor:
    diff_chunks = list(executor.map(compute_difference, wt_chunks, ko_chunks))

# Merge results
diff_matrix = np.vstack(diff_chunks)

# Function to plot heatmap
def plot_heatmap(matrix, output_file, title, vmin, vmax, dpi):
    plt.figure(figsize=(5, 20))  # Longer and narrower plot
    ax = sns.heatmap(matrix, cmap="coolwarm", vmin=vmin, vmax=vmax, xticklabels=False, yticklabels=False)

    # Rasterize heatmap to reduce PDF size
    img = ax.collections[0]
    img.set_rasterized(True)

    # Reduce colorbar labels for efficiency
    cbar = ax.collections[0].colorbar
    cbar.set_label("Difference (KO - WT)")
    cbar.locator = ticker.MaxNLocator(nbins=3)
    cbar.update_ticks()

    plt.title(title)
    plt.xlabel("Bins")
    plt.ylabel("Regions")

    # Save as optimized PDF
    plt.savefig(output_file, dpi=dpi, format="pdf", bbox_inches="tight", pad_inches=0.1)
    plt.close()
    print(f"Heatmap saved as optimized PDF: {output_file}")

# Plot the heatmap
plot_heatmap(diff_matrix, output_pdf_file, "Difference Heatmap (KO - WT)", zmin, zmax, dpi)
