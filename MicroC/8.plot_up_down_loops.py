#### step1: merge loops at different resolutions and keep the overlaped ones at finest resolution
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import straw
import csv
import random
import math
from tqdm import tqdm
import concurrent.futures

# Define Hi-C files for WT and KO
hic_wt = "WT.hic"
hic_ko = "KO.hic"

# Path to the BEDPE file and output CSV
bedpe_file = "hic_loop_all_sorted.txt"
csv_file = "hic_loop_all_sorted_with_hic_signal.csv"

# Sample size: Set to 100 for testing
sample_size = None  ### None or a number

# Function to determine the appropriate resolution based on loop size
def get_resolution(start1, end1, start2, end2):
    loop_size = max(end1 - start1, end2 - start2)
    if loop_size >= 10000:
        return 10000
    elif loop_size >= 5000:
        return 5000
    elif loop_size >= 2000:
        return 2000
    elif loop_size >= 1000:
        return 1000
    else:
        return 500

# Log2 transformation with handling for zero values
def safe_log2(value):
    return math.log2(value) if value > 0 else 0

# Function to calculate the sum of interaction values
def calculate_sum(result):
    return sum(record.counts for record in result)

# Function to query Hi-C data with appropriate resolution
def query_hic_data(hic_file, chrom, start1, end1, start2, end2):
    resolution = get_resolution(start1, end1, start2, end2)
    region1 = f"{chrom}:{start1}:{end1}"
    region2 = f"{chrom}:{start2}:{end2}"
    result = straw.straw("oe", "KR", hic_file, region1, region2, "BP", resolution)
    individual_values = [(record.binX, record.binY, record.counts) for record in result]
    total_sum = calculate_sum(result)
    return individual_values, total_sum

# Task function for multiprocessing
def task(index, hic_file, region):
    individual_values, total_sum = query_hic_data(hic_file, *region)
    return (index, individual_values, total_sum)  # Return as tuple

# Parallel processing with batch handling to avoid resource overload
def process_in_batches(hic_file, regions, desc, batch_size=1000, max_workers=16):
    all_results = []
    for i in range(0, len(regions), batch_size):
        batch = regions[i:i + batch_size]
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(task, j, hic_file, region): j for j, region in enumerate(batch)
            }
            results = [None] * len(batch)  # Pre-allocate results list
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc=f"{desc} (Batch {i // batch_size + 1})"):
                try:
                    index, individual_values, total_sum = future.result()
                    results[index] = (individual_values, total_sum)
                except Exception as e:
                    print(f"Error processing region: {e}")
                    results[index] = ([], 0)  # Handle error gracefully
        all_results.extend(results)
    return all_results

# Read BEDPE file and extract regions
def read_bedpe_file(file_path, sample_size=None):
    regions = []
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split()
            chrom = parts[0].replace('chr', '')  # Remove 'chr' prefix
            start1, end1 = int(parts[1]), int(parts[2])
            start2, end2 = int(parts[4]), int(parts[5])
            regions.append((chrom, start1, end1, start2, end2))
    return random.sample(regions, sample_size) if sample_size else regions

# Main execution
regions = read_bedpe_file(bedpe_file, sample_size)

# Process WT and KO regions in batches
wt_results = process_in_batches(hic_wt, regions, "Processing WT")
ko_results = process_in_batches(hic_ko, regions, "Processing KO")

# Save results to CSV
with open(csv_file, mode="w", newline="") as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow([
        "chr1", "start1", "end1", "chr2", "start2", "end2",
        "wt_individual_values", "wt_sum_value", "log2_wt_sum",
        "ko_individual_values", "ko_sum_value", "log2_ko_sum"
    ])
    for i, region in enumerate(regions):
        wt_values, wt_sum = wt_results[i]
        ko_values, ko_sum = ko_results[i]
        log2_wt = safe_log2(wt_sum)
        log2_ko = safe_log2(ko_sum)
        writer.writerow([
            f"chr{region[0]}", region[1], region[2], f"chr{region[0]}", region[3], region[4],
            wt_values, wt_sum, log2_wt,
            ko_values, ko_sum, log2_ko
        ])

# Load the CSV file for plotting
data = pd.read_csv(csv_file, delimiter='\t')

# Extract relevant columns
log2_wt_values = data['log2_wt_sum']
log2_ko_values = data['log2_ko_sum']

# Calculate total loops and loops beyond threshold lines
total_loops = len(data)
above_line_05 = sum(log2_ko_values > log2_wt_values + 0.5)
below_line_05 = sum(log2_ko_values < log2_wt_values - 0.5)

# Plotting the dot plot with log2-transformed values
plt.figure(figsize=(8, 8))
plt.scatter(log2_wt_values, log2_ko_values, color="blue", alpha=0.5)

# Add separation lines
x_vals = np.linspace(min(log2_wt_values.min(), log2_ko_values.min()), 
                     max(log2_wt_values.max(), log2_ko_values.max()), 100)
plt.plot(x_vals, x_vals, 'k--', alpha=0.8, label='y = x')
plt.plot(x_vals, x_vals + 0.5, 'r--', alpha=0.6, label='y = x + 0.5')
plt.plot(x_vals, x_vals - 0.5, 'g--', alpha=0.6, label='y = x - 0.5')

# Add loop counts to the plot
plt.text(0.95, 0.05, 
         f'Total Loops: {total_loops}\n'
         f'KO > WT + 0.5: {above_line_05}\n'
         f'WT > KO + 0.5: {below_line_05}', 
         transform=plt.gca().transAxes, fontsize=12, 
         verticalalignment='bottom', horizontalalignment='right', 
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.6))

# Set labels and title
plt.xlabel("log2(WT Sum Value)")
plt.ylabel("log2(KO Sum Value)")
plt.title("log2(WT) vs log2(KO) Interaction Strength (Sampled Regions)")

# Add legend
plt.legend()

# Save the plot as PDF
plt.tight_layout()
plt.savefig("hic_wt_vs_ko_interaction_strength_separated.pdf", format='pdf')

# Show the plot
plt.show()


### STEP2: assign each anchor to genomic annotations(CTCF, Enhancer, Promoter, etc)
import pandas as pd
import pyranges as pr

# File paths
csv_file = "hic_loop_all_sorted_with_hic_signal.csv"
bed_file = "hct_all_region_ctcf_enhancer_sort.bed"  ### annotated regions from USCS genome browser + CTCF + Promoter + Enhancer

# Read the CSV file
columns = [
    "chr1", "start1", "end1", "chr2", "start2", "end2",
    "wt_individual_values", "wt_sum_value", "log2_wt_sum",
    "ko_individual_values", "ko_sum_value", "log2_ko_sum"
]
df = pd.read_csv(csv_file, sep="\t", header=0, names=columns)

# Load the BED file
bed_df = pd.read_csv(bed_file, sep="\t", header=None)
num_bed_columns = bed_df.shape[1]
if num_bed_columns >= 4:
    bed_columns = ["Chromosome", "Start", "End", "Feature"] + [
        f"Extra_{i}" for i in range(5, num_bed_columns + 1)
    ]
else:
    raise ValueError("The BED file must have at least 4 columns.")
bed_df.columns = bed_columns[:num_bed_columns]

# Normalize Chromosome names to ensure they match
df["chr1"] = df["chr1"].apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else x)
df["chr2"] = df["chr2"].apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else x)
bed_df["Chromosome"] = bed_df["Chromosome"].apply(lambda x: f"chr{x}" if not str(x).startswith("chr") else x)

# Initialize columns for intersected classes
df["left_classes"] = None
df["right_classes"] = None

# Function to collect intersected classes for a subset of the BED file
def collect_classes(loop_regions, bed_subset):
    loop_overlap = loop_regions.join(bed_subset)
    overlap_df = loop_overlap.as_df()
    grouped = overlap_df.groupby(["Chromosome", "Start", "End"]).agg({"Feature": lambda x: ",".join(sorted(x.unique()))}).reset_index()
    return grouped.rename(columns={"Feature": "Intersected_Classes"})

# Process chromosome by chromosome
unique_chromosomes = df["chr1"].unique()

for chrom in unique_chromosomes:
    print(f"Processing chromosome: {chrom}")

    # Subset BED file for the current chromosome
    bed_subset = pr.PyRanges(bed_df[bed_df["Chromosome"] == chrom])

    # Subset loops for the current chromosome
    left_subset = pr.PyRanges(df[df["chr1"] == chrom][["chr1", "start1", "end1"]].rename(columns={
        "chr1": "Chromosome", "start1": "Start", "end1": "End"
    }))
    right_subset = pr.PyRanges(df[df["chr2"] == chrom][["chr2", "start2", "end2"]].rename(columns={
        "chr2": "Chromosome", "start2": "Start", "end2": "End"
    }))

    # Intersect left and right regions
    left_classes = collect_classes(left_subset, bed_subset)
    right_classes = collect_classes(right_subset, bed_subset)

    # Merge results back to the main DataFrame
    df = df.merge(left_classes, how="left", left_on=["chr1", "start1", "end1"], right_on=["Chromosome", "Start", "End"]).rename(columns={"Intersected_Classes": "temp_left_classes"})
    df["left_classes"] = df["left_classes"].combine_first(df["temp_left_classes"])
    df = df.drop(columns=["Chromosome", "Start", "End", "temp_left_classes"])

    df = df.merge(right_classes, how="left", left_on=["chr2", "start2", "end2"], right_on=["Chromosome", "Start", "End"]).rename(columns={"Intersected_Classes": "temp_right_classes"})
    df["right_classes"] = df["right_classes"].combine_first(df["temp_right_classes"])
    df = df.drop(columns=["Chromosome", "Start", "End", "temp_right_classes"])

# Save the updated DataFrame to a file
output_file = "HCT_wt_ko_original_with_classes_by_annotation.csv"
df.to_csv(output_file, sep="\t", index=False)

print(f"File with original matrix and new classes saved to: {output_file}")

# Display a sample of the final result
print(df.head())


### step3: prioritize the loop categories
# Load the CSV file with the classes already calculated
csv_file = "HCT_wt_ko_original_with_classes_by_annotation.csv"
df = pd.read_csv(csv_file, sep="\t")

def prioritize_classes(classes):
    if pd.isna(classes):  # If no classes, return 'none'
        return "none"
    class_list = classes.split(",")  # Split comma-separated classes
    class_str = ",".join(class_list)  # Combine into a single string for substring search
    if "ctcf" in class_str:  # Priority 1: CTCF
        return "CTCF"
    elif "tss" in class_str:  # Priority 2: tss
        return "P"  # Rename tss to P
    elif any(sub in class_str for sub in ["Stitch", "Super"]):  # Priority 3: Stitch or Super
        return "E"
    else:  # If no relevant substrings found
        return "none"


# Apply the prioritization function to left_classes and right_classes
df["left_priority"] = df["left_classes"].apply(prioritize_classes)
df["right_priority"] = df["right_classes"].apply(prioritize_classes)

# Save the updated DataFrame with the full matrix
output_file = "HEK_wt_ko_with_prioritized_annotation_classes_full_matrix.csv"
df.to_csv(output_file, sep="\t", index=False)

print(f"File with prioritized classes and full matrix saved to: {output_file}")

# Display a sample of the result
print(df.head())



### step4: plotting J bound CTCF loops, annotate anchors with interested proteins(J or Z)

import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde
import matplotlib as mpl
# Set PDF font type for editable text
mpl.rcParams['pdf.fonttype'] = 42

# File paths
csv_file = "HCT_wt_ko_with_prioritized_annotation_classes_full_matrix.csv"
j_bed_file = "HCT116_J_chip_85_87_sort_dedup_q005_peaks.narrowPeak"

# Read the CSV file
df = pd.read_csv(csv_file, sep="\t")

# Load the J BED file as a PyRanges object
j_bed = pr.read_bed(j_bed_file)

# Create PyRanges objects for "left" and "right" regions
left_regions = pr.PyRanges(df[["chr1", "start1", "end1"]].rename(columns={
    "chr1": "Chromosome", "start1": "Start", "end1": "End"
}))
right_regions = pr.PyRanges(df[["chr2", "start2", "end2"]].rename(columns={
    "chr2": "Chromosome", "start2": "Start", "end2": "End"
}))

# Check overlaps for "left" and "right"
left_overlap = left_regions.intersect(j_bed).df
right_overlap = right_regions.intersect(j_bed).df

# Create a set of overlapping indices for "left" and "right"
left_indices = set(left_overlap.index)
right_indices = set(right_overlap.index)

# Add J overlap results to the DataFrame
df["left_J"] = ["yes" if i in left_indices else "no" for i in df.index]
df["right_J"] = ["yes" if i in right_indices else "no" for i in df.index]

# Filter loops with "J" on either side
j_loops = df[(df["left_J"] == "yes") | (df["right_J"] == "yes")]

# Update total_loops to reflect only J-binding loops
total_loops = len(j_loops)

# Prepare for plotting
log2_wt_values = j_loops["log2_wt_sum"]
log2_ko_values = j_loops["log2_ko_sum"]

# Create masks
is_ctcf_loop = (j_loops["left_priority"] == "CTCF") | (j_loops["right_priority"] == "CTCF")
ctcf_wt = log2_wt_values[is_ctcf_loop]
ctcf_ko = log2_ko_values[is_ctcf_loop]

# Calculate density for CTCF loops
ctcf_xy = np.vstack([ctcf_wt, ctcf_ko])
ctcf_density = gaussian_kde(ctcf_xy)(ctcf_xy)

# Sort points by density
idx = ctcf_density.argsort()
ctcf_wt, ctcf_ko, ctcf_density = ctcf_wt.iloc[idx], ctcf_ko.iloc[idx], ctcf_density[idx]

# Define colormap
colors = ["#1f78b4", "#00a9cf", "#5aae61", "#fee08b", "#fdae61", "#d73027"]
cmap = LinearSegmentedColormap.from_list("blue_cyan_green_yellow_orange_red", colors, N=256)

# Define a function to calculate upregulated and downregulated loops for a given mask
def calculate_up_down_loops(mask):
    up = (log2_ko_values[mask] > log2_wt_values[mask] + 0.5).sum()
    down = (log2_ko_values[mask] < log2_wt_values[mask] - 0.5).sum()
    return up, down

# Calculate upregulated and downregulated loops for EP loops
ep_up_loops, ep_down_loops = calculate_up_down_loops(is_ctcf_loop)

# Filter out points where x=0 or y=0
# Prepare for plotting from j_loops
log2_wt_values_plot = j_loops["log2_wt_sum"]
log2_ko_values_plot = j_loops["log2_ko_sum"]
valid_indices = (log2_wt_values_plot != 0) & (log2_ko_values_plot != 0)
log2_wt_values_plot_filtered = log2_wt_values_plot[valid_indices]
log2_ko_values_plot_filtered = log2_ko_values_plot[valid_indices]

# Also filter CTCF points
log2_wt_ctcf = log2_wt_values[is_ctcf_loop]
log2_ko_ctcf = log2_ko_values[is_ctcf_loop]

# density
ctcf_xy = np.vstack([log2_wt_ctcf, log2_ko_ctcf])
ctcf_z = gaussian_kde(ctcf_xy)(ctcf_xy)

# sort index
ctcf_idx = ctcf_z.argsort()
log2_wt_ctcf = log2_wt_ctcf.iloc[ctcf_idx]
log2_ko_ctcf = log2_ko_ctcf.iloc[ctcf_idx]
ctcf_z = ctcf_z[ctcf_idx]

valid_ctcf_indices = (log2_wt_ctcf != 0) & (log2_ko_ctcf != 0)
log2_wt_ctcf_filtered = log2_wt_ctcf[valid_ctcf_indices]
log2_ko_ctcf_filtered = log2_ko_ctcf[valid_ctcf_indices]
ctcf_z_filtered = ctcf_z[valid_ctcf_indices]

# define colormap
colors = ["#1f78b4", "#00a9cf", "#5aae61", "#fee08b", "#fdae61", "#d73027"]
cmap = LinearSegmentedColormap.from_list("blue_cyan_green_yellow_orange_red", colors, N=256)

# Update counting variables
total_loops_filtered = valid_indices.sum()
ep_loops_filtered = valid_ctcf_indices.sum()
ep_up_loops_filtered = ((log2_ko_ctcf_filtered - log2_wt_ctcf_filtered) > 0.5).sum()
ep_down_loops_filtered = ((log2_wt_ctcf_filtered - log2_ko_ctcf_filtered) > 0.5).sum()

fig, ax = plt.subplots(figsize=(8, 8))  # Ensure the figure itself is square

# Plot J-bound loops with density coloring (gray scale)
scatter_all = ax.scatter(
    log2_wt_values_plot_filtered,
    log2_ko_values_plot_filtered,
    color="grey",
    alpha=0.5,
    s=18,
    label="J-Bound Loops (Density)",
     edgecolors='none'
)

scatter_all.set_rasterized(True)  # Rasterize the "all" points
# Highlight CTCF loops with density coloring
scatter_ctcf = ax.scatter(
    log2_wt_ctcf_filtered,
    log2_ko_ctcf_filtered,
    c=ctcf_z_filtered,
    cmap=cmap,
    alpha=0.8,
    s=30,
    label="CTCF Loops",
     edgecolors='none'
)
scatter_ctcf.set_rasterized(True)  # Rasterize the "CTCF" points




# Add colorbar for density of J-bound loops
cbar_all = plt.colorbar(scatter_ctcf)
cbar_all.set_label("Density (J-Bound CTCF Loops)", fontsize=12)  # Label for the colorbar
cbar_all.ax.tick_params(labelsize=10)  # Adjust tick label size

# seperate lines
x_vals = np.linspace(-2, 12, 100)  # can adjust the values
plt.plot(x_vals, x_vals, color='black', linestyle='--', linewidth=0.8, alpha=1, label='y = x')
plt.plot(x_vals, x_vals + 0.5, color='black', linestyle=':', linewidth=0.7, alpha=1, label='y = x + 0.5')
plt.plot(x_vals, x_vals - 0.5, color='black', linestyle=':', linewidth=0.7, alpha=1, label='y = x - 0.5')


plt.xlim(-0.2, 10)
plt.ylim(-2, 10)


# number panel
plt.text(0.95, 0.05, 
         f'Total J Loops: {total_loops_filtered}\n'
         f'CTCF Loops: {ep_loops_filtered}\n'
         f'CTCF Up Loops: {ep_up_loops_filtered}\n'
         f'CTCF Down Loops: {ep_down_loops_filtered}', 
         transform=plt.gca().transAxes, fontsize=12, 
         verticalalignment='bottom', horizontalalignment='right', 
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.6), color='black')

# add title and axis
plt.xlabel("log2(WT Sum Value)")
plt.ylabel("log2(KO Sum Value)")
plt.title("log2(WT) vs log2(KO) Interaction Strength (J-Bound Loops with CTCF Highlighted)")

# Add legend
ax.legend()

# Adjust aspect ratio so physical lengths of axes are equal
aspect_ratio = (10 - 0) /(10 - (-2))   # y_range / x_range
fig.set_figheight(fig.get_figwidth() * aspect_ratio)  # Adjust height relative to width

# save
plt.tight_layout()
plt.savefig("/hek_J_1B_wt_vs_ko_CTCF_density_fixed.pdf", format="pdf")
plt.show()

### step5: plotting EP loops, annotate anchors with interested proteins(J or Z)
import pandas as pd
import pyranges as pr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import gaussian_kde
import matplotlib as mpl
# Set PDF font type for editable text
mpl.rcParams['pdf.fonttype'] = 42

# File paths
csv_file = "HEK_znf_1B_wt_ko_with_prioritized_classes_full_matrix.csv"
j_bed_file = "HCT116_J_chip_85_87_sort_dedup_q005_peaks.narrowPeak"

# Read the CSV file
df = pd.read_csv(csv_file, sep="\t")

# Load the J BED file as a PyRanges object
j_bed = pr.read_bed(j_bed_file)

# Create PyRanges objects for "left" and "right" regions
left_regions = pr.PyRanges(df[["chr1", "start1", "end1"]].rename(columns={
    "chr1": "Chromosome", "start1": "Start", "end1": "End"
}))
right_regions = pr.PyRanges(df[["chr2", "start2", "end2"]].rename(columns={
    "chr2": "Chromosome", "start2": "Start", "end2": "End"
}))

# Check overlaps for "left" and "right"
left_overlap = left_regions.intersect(j_bed).df
right_overlap = right_regions.intersect(j_bed).df

# Create a set of overlapping indices for "left" and "right"
left_indices = set(left_overlap.index)
right_indices = set(right_overlap.index)

# Add J overlap results to the DataFrame
df["left_J"] = ["yes" if i in left_indices else "no" for i in df.index]
df["right_J"] = ["yes" if i in right_indices else "no" for i in df.index]

# Filter loops with "J" on either side
j_loops = df[(df["left_J"] == "yes") | (df["right_J"] == "yes")]

# Update total_loops to reflect only J-binding loops
total_loops = len(j_loops)

# Prepare for plotting
log2_wt_values = j_loops["log2_wt_sum"]
log2_ko_values = j_loops["log2_ko_sum"]

# Create masks
is_EP_loop = (
    (j_loops["left_priority"].isin(["E", "P"])) | (j_loops["right_priority"].isin(["E", "P"]))
) & ~((j_loops["left_priority"] == "CTCF") | (j_loops["right_priority"] == "CTCF"))

EP_wt = log2_wt_values[is_EP_loop]
EP_ko = log2_ko_values[is_EP_loop]

# Calculate density for EP loops
EP_xy = np.vstack([EP_wt, EP_ko])
EP_density = gaussian_kde(EP_xy)(EP_xy)

# Sort points by density
idx = EP_density.argsort()
EP_wt, EP_ko, EP_density = EP_wt.iloc[idx], EP_ko.iloc[idx], EP_density[idx]

# Define colormap
colors = ["#0461cf", "#b2326c", "#5aae61", "#fee08b", "#fdae61", "#d73027"]
cmap = LinearSegmentedColormap.from_list("blue_cyan_green_yellow_orange_red", colors, N=256)

# Define a function to calculate upregulated and downregulated loops for a given mask
def calculate_up_down_loops(mask):
    up = (log2_ko_values[mask] > log2_wt_values[mask] + 0.5).sum()
    down = (log2_ko_values[mask] < log2_wt_values[mask] - 0.5).sum()
    return up, down

# Calculate upregulated and downregulated loops for EP loops
ep_up_loops, ep_down_loops = calculate_up_down_loops(is_EP_loop)

# Filter out points where x=0 or y=0
# Prepare for plotting from j_loops
log2_wt_values_plot = j_loops["log2_wt_sum"]
log2_ko_values_plot = j_loops["log2_ko_sum"]
valid_indices = (log2_wt_values_plot != 0) & (log2_ko_values_plot != 0)
log2_wt_values_plot_filtered = log2_wt_values_plot[valid_indices]
log2_ko_values_plot_filtered = log2_ko_values_plot[valid_indices]

# Also filter EP points
log2_wt_EP = log2_wt_values[is_EP_loop]
log2_ko_EP = log2_ko_values[is_EP_loop]


EP_xy = np.vstack([log2_wt_EP, log2_ko_EP])
EP_z = gaussian_kde(EP_xy)(EP_xy)


EP_idx = EP_z.argsort()
log2_wt_EP = log2_wt_EP.iloc[EP_idx]
log2_ko_EP = log2_ko_EP.iloc[EP_idx]
EP_z = EP_z[EP_idx]

valid_EP_indices = (log2_wt_EP != 0) & (log2_ko_EP != 0)
log2_wt_EP_filtered = log2_wt_EP[valid_EP_indices]
log2_ko_EP_filtered = log2_ko_EP[valid_EP_indices]
EP_z_filtered = EP_z[valid_EP_indices]


colors = ["#2b388a", "#8babf1", "#5aae61", "#fee08b", "#fdae61", "#d73027"]

cmap = LinearSegmentedColormap.from_list("blue_cyan_green_yellow_orange_red", colors, N=256)

# Update counting variables
total_loops_filtered = valid_indices.sum()
ep_loops_filtered = valid_EP_indices.sum()
ep_up_loops_filtered = ((log2_ko_EP_filtered - log2_wt_EP_filtered) > 0.5).sum()
ep_down_loops_filtered = ((log2_wt_EP_filtered - log2_ko_EP_filtered) > 0.5).sum()


fig, ax = plt.subplots(figsize=(8, 8))  # Ensure the figure itself is square
# Plot J-bound loops with density coloring (gray scale)
scatter_all = ax.scatter(
    log2_wt_values_plot_filtered,
    log2_ko_values_plot_filtered,
    color="grey",
    alpha=0.5,
    s=18,
    label="J-Bound Loops (Density)",
     edgecolors='none'
)

scatter_all.set_rasterized(True)  # Rasterize the "all" points
# Highlight EP loops with density coloring
scatter_EP = ax.scatter(
    log2_wt_EP_filtered,
    log2_ko_EP_filtered,
    c=EP_z_filtered,
    cmap=cmap,
    alpha=0.8,
    s=30,
    label="EP Loops",
     edgecolors='none'
)
scatter_EP.set_rasterized(True)  # Rasterize the "EP" points
# Add colorbar for density of J-bound loops
cbar_all = plt.colorbar(scatter_EP)
cbar_all.set_label("Density (J-Bound EP Loops)", fontsize=12)  # Label for the colorbar
cbar_all.ax.tick_params(labelsize=10)  # Adjust tick label size

x_vals = np.linspace(-2, 12, 100)  # 根据需要调整 x 和 y 的范围
plt.plot(x_vals, x_vals, color='black', linestyle='--', linewidth=0.8, alpha=1, label='y = x')
plt.plot(x_vals, x_vals + 0.5, color='black', linestyle=':', linewidth=0.7, alpha=1, label='y = x + 0.5')
plt.plot(x_vals, x_vals - 0.5, color='black', linestyle=':', linewidth=0.7, alpha=1, label='y = x - 0.5')


plt.xlim(-0.2, 10)
plt.ylim(-2, 10)


plt.text(0.95, 0.05, 
         f'Total J Loops: {total_loops_filtered}\n'
         f'EP Loops: {ep_loops_filtered}\n'
         f'EP Up Loops: {ep_up_loops_filtered}\n'
         f'EP Down Loops: {ep_down_loops_filtered}', 
         transform=plt.gca().transAxes, fontsize=12, 
         verticalalignment='bottom', horizontalalignment='right', 
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.6), color='black')


plt.xlabel("log2(WT Sum Value)")
plt.ylabel("log2(KO Sum Value)")
plt.title("log2(WT) vs log2(KO) Interaction Strength (J-Bound Loops with EP Highlighted)")
ax.legend()

# Adjust aspect ratio so physical lengths of axes are equal
aspect_ratio = (10 - 0) /(10 - (-2))   # y_range / x_range
fig.set_figheight(fig.get_figwidth() * aspect_ratio)  # Adjust height relative to width

# save
plt.tight_layout()
plt.savefig("hct_J_1B_wt_vs_ko_J_loops_EP_density_fixed.pdf", format="pdf")
plt.show()


### step6: export significantly upregulated/downregulated loops

import pandas as pd
import pyranges as pr
import numpy as np
from scipy.stats import gaussian_kde

csv_file = "HEK_znf_1B_wt_ko_with_prioritized_classes_full_matrix.csv"
j_bed_file = "HCT116_J_chip_85_87_sort_dedup_q005_peaks.narrowPeak"
output_dir = "./"

df = pd.read_csv(csv_file, sep="\t")
j_bed = pr.read_bed(j_bed_file)

left_regions = pr.PyRanges(df[["chr1", "start1", "end1"]].rename(columns={"chr1": "Chromosome", "start1": "Start", "end1": "End"}))
right_regions = pr.PyRanges(df[["chr2", "start2", "end2"]].rename(columns={"chr2": "Chromosome", "start2": "Start", "end2": "End"}))


left_overlap = left_regions.intersect(j_bed).df
right_overlap = right_regions.intersect(j_bed).df
left_indices = set(left_overlap.index)
right_indices = set(right_overlap.index)

df["left_J"] = ["yes" if i in left_indices else "no" for i in df.index]
df["right_J"] = ["yes" if i in right_indices else "no" for i in df.index]


j_loops = df[(df["left_J"] == "yes") | (df["right_J"] == "yes")]

is_ctcf_loop = (j_loops["left_priority"] == "CTCF") | (j_loops["right_priority"] == "CTCF")
ctcf_loops_df = j_loops[is_ctcf_loop].copy()

log2_wt_ctcf = ctcf_loops_df["log2_wt_sum"]
log2_ko_ctcf = ctcf_loops_df["log2_ko_sum"]
valid_ctcf_indices = (log2_wt_ctcf != 0) & (log2_ko_ctcf != 0)
ctcf_loops_filtered_df = ctcf_loops_df[valid_ctcf_indices].copy()
log2_wt_ctcf_filtered = log2_wt_ctcf[valid_ctcf_indices]
log2_ko_ctcf_filtered = log2_ko_ctcf[valid_ctcf_indices]


up_mask = (log2_ko_ctcf_filtered - log2_wt_ctcf_filtered) > 0.5
down_mask = (log2_wt_ctcf_filtered - log2_ko_ctcf_filtered) > 0.5

upregulated_loops = ctcf_loops_filtered_df[up_mask].copy()
downregulated_loops = ctcf_loops_filtered_df[down_mask].copy()


up_path = f"{output_dir}/upregulated_ctcf_loops.bedpe"
down_path = f"{output_dir}/downregulated_ctcf_loops.bedpe"

upregulated_loops[["chr1", "start1", "end1", "chr2", "start2", "end2"]].to_csv(up_path, sep="\t", index=False)
downregulated_loops[["chr1", "start1", "end1", "chr2", "start2", "end2"]].to_csv(down_path, sep="\t", index=False)

print(f"upregulated CTCF loops exported to: {up_path}")
print(f"downregulated CTCF loops exported to: {down_path}")
