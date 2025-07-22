import cooler
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import seaborn as sns
import time  # Import the time module
import sys
import subprocess
from multiprocessing import Pool, cpu_count
import os

# Start timing
start_time = time.time()
# Function to compute interaction sums for a given row index
def compute_interaction_sum(i, window_size, matrix):
    start_index = max(0, i - window_size)
    end_index = min(len(matrix), i + window_size + 1)
    # Exclude the diagonal value by summing the upper and lower parts separately
    sum_upper = np.nansum(matrix[start_index:i-2, i])
    sum_lower = np.nansum(matrix[i+2:end_index, i])
    strand = '+' if sum_upper > sum_lower else '-'
    return sum_upper + sum_lower, sum_upper, sum_lower, strand

import cooler




cool_file_path = sys.argv[1]

clr = cooler.Cooler(cool_file_path)
argv_1 = sys.argv[1]
resolution_str = argv_1.split('resolutions/')[1]
resolution_number = int(resolution_str)



# window_size = 1000  # 5kb on each side
window_size1 = sys.argv[2]
prefix = sys.argv[3]
# prefix = bed_filename.split('.')
threshold = sys.argv[4]
region = sys.argv[5] if len(sys.argv) > 5 else None  # Optional region argument
window_size1 = int(window_size1)
threshold = float(threshold)
# Fetch the matrix for the specified region
# matrix = clr.matrix(balance=False).fetch(region)


# Load chromosome sizes
chrom_sizes = pd.read_csv('hg38.clean.chrom.sizes', sep='\t', header=None, names=['chrom', 'size'])


# print(chrom_sizes)

# Set sliding window parameters
window_size = 20000000  # 20 million bp
step_size = 20000000  # 10 million bp
bed_data = []


def process_chromosome(row, window_size, step_size, threshold, prefix, clr, resolution_number):
    """
    Function to process a single chromosome for stripe calculations.
    """
    chromosome = row['chrom']
    chrom_size = row['size']
    start = 0
    while start < chrom_size:
        end = min(start + window_size, chrom_size)
        region = f"{chromosome}:{start}-{end}"
        print(region)
        chromosome, positions = region.split(':')
        prefix1 = prefix + "_" + region + ".bed1"
        prefix2 = prefix + "_" + region + ".bed2"
        matrix = clr.matrix(balance=True).fetch(region)  ### change this if no hic matrix balance

        # Calculate interaction sums and strands separately
        interaction_sums_and_strands = [compute_interaction_sum(i, window_size, matrix) for i in range(len(matrix))]
        interaction_sums = np.array([item[0] for item in interaction_sums_and_strands])  # Numeric data for statistical calculations
        strands = [item[3] for item in interaction_sums_and_strands]
        sum_uppers = np.array([item[1] for item in interaction_sums_and_strands])
        sum_lowers = np.array([item[2] for item in interaction_sums_and_strands])

        if interaction_sums.size > 0:
            mean_signal = np.mean(interaction_sums)
            std_dev_signal = np.std(interaction_sums)
            prominence_value = mean_signal + std_dev_signal * threshold

            # Calculate genomic coordinates and output the detailed information for each bin
            region_start = int(positions.split('-')[0])
            bins = range(len(matrix))
            chromosome_data = [chromosome] * len(matrix)
            bin_starts = [region_start + (bin * clr.binsize) for bin in bins]
            bin_ends = [min(start + resolution_number, chrom_size) for start in bin_starts]
            sum_uppers_for_bins = [sum_upper for sum_upper in sum_uppers]
            sum_lowers_for_bins = [sum_lower for sum_lower in sum_lowers]

            # Create and save BED format for all bins
            bed_df_all_bins = pd.DataFrame({
                'chrom': chromosome_data,
                'start': bin_starts,
                'end': bin_ends,
                'sum_upper': sum_uppers_for_bins,
                'sum_lower': sum_lowers_for_bins
            })
            bed_df_all_bins.to_csv(prefix2, sep='\t', index=False, header=False)

            # Find peaks
            peaks, _ = find_peaks(interaction_sums, prominence=prominence_value)
            peak_coordinates = [region_start + (peak * clr.binsize) for peak in peaks]
            peak_strands = [strands[peak] for peak in peaks]
            uper = [sum_uppers[peak] for peak in peaks]
            lower = [sum_lowers[peak] for peak in peaks]

            # Create and save BED format for peaks
            bed_df = pd.DataFrame({
                'chrom': [chromosome] * len(peak_coordinates),
                'start': [coord + int(resolution_number * 0.1) for coord in peak_coordinates],
                'end': [coord + int(resolution_number * 1.1) for coord in peak_coordinates],
                'strand': peak_strands,
                'sum_upper': uper,
                'sum_lower': lower
            })
            bed_df.to_csv(prefix1, sep='\t', index=False, header=False)
        else:
            print(f"No interaction sums calculated for {region}. Check your matrix and window size.")
        
        start += step_size


# Main execution logic
if len(sys.argv) > 5 and sys.argv[5]:
    # Specific chromosome
    chromosome_data = chrom_sizes[chrom_sizes['chrom'] == sys.argv[5]]
else:
    # Whole genome
    chromosome_data = chrom_sizes

for index, row in chromosome_data.iterrows():
    process_chromosome(row, window_size, step_size, threshold, prefix, clr, resolution_number)


# # Plot the heatmap
# plt.figure(figsize=(10, 10))
# ax = sns.heatmap(matrix, cmap='Reds', vmax=0.75, square=True)

# # Set title and labels for the heatmap
# plt.title(f"Heatmap for region {region}")
# plt.xlabel('Position along chromosome 19')
# plt.ylabel('Position along chromosome 19')

# # Calculate tick positions and labels
# step_size = 5  # Bin step size
# binsize = 200  # 200bp per bin, adjust as necessary
# tick_positions = np.arange(0, len(matrix), step_size)
# tick_labels = [f'{pos}' for pos in tick_positions]

# # Set the tick positions and labels on both axes
# ax.set_xticks(tick_positions)
# ax.set_xticklabels(tick_labels, rotation=90)
# ax.set_yticks(tick_positions)
# ax.set_yticklabels(tick_labels)

# # Save the plot to a PDF file
# pdf_filename = './heatmap_region1.pdf'
# plt.savefig(pdf_filename, bbox_inches='tight')
# plt.close()  # Close the plot to free up memory

# print(f"Heatmap saved to {pdf_filename}.")

# The number of bins corresponding to 5kb on either side of the diagonal


# # Calculate interaction sums
# interaction_sums = [compute_interaction_sum(i, window_size, matrix) for i in range(len(matrix))]
prefix3 = prefix+"_bw.bed"
prefix4 = prefix+"_peak.bed"
prefix5 = prefix+"_bw_plus.bed"
prefix6 = prefix+"_bw_minus.bed"


command = f"cat *.bed1 >> {prefix4}"
result = subprocess.run(command, shell=True, text=True, check=True)


command1 = f"rm *.bed1"
result = subprocess.run(command1, shell=True, text=True, check=True)

command2 = f"cat *.bed2 >> {prefix3}"
result = subprocess.run(command2, shell=True, text=True, check=True)

command3 = f"rm *.bed2"
result = subprocess.run(command3, shell=True, text=True, check=True)


# Remove duplicate lines in prefix3
command_uniq = f"sort -k1,1 -k2,2n {prefix3} | uniq > {prefix3}_uniq"
result = subprocess.run(command_uniq, shell=True, text=True, check=True)

# Process positive and negative signal separately
command4 = f"awk -v OFS='\\t' '{{print $1,$2,$3,$4}}' {prefix3}_uniq > {prefix5}"
command5 = f"awk -v OFS='\\t' '{{print $1,$2,$3,-1*$5}}' {prefix3}_uniq > {prefix6}"
result = subprocess.run(command4, shell=True, text=True, check=True)
result = subprocess.run(command5, shell=True, text=True, check=True)

# # Sort the resulting bed files for proper conversion
# command6 = f"sort -k1,1 -k2,2n {prefix5} > {prefix}_bw_up_sort.bed"
# command7 = f"sort -k1,1 -k2,2n {prefix6} > {prefix}_bw_down_sort.bed"
# result = subprocess.run(command6, shell=True, text=True, check=True)
# result = subprocess.run(command7, shell=True, text=True, check=True)

# Convert the sorted BED files to BigWig format using bedGraphToBigWig
chrom_sizes = "hg38.clean.chrom.sizes"
command8 = f"bedGraphToBigWig {prefix5} {chrom_sizes} {prefix}_bw_up_sort.bw"
command9 = f"bedGraphToBigWig {prefix6} {chrom_sizes} {prefix}_bw_down_sort.bw"
result = subprocess.run(command8, shell=True, text=True, check=True)
result = subprocess.run(command9, shell=True, text=True, check=True)

# Sum column 4 and 5 in {prefix3}_uniq to make {prefix3}_uniq_updown.bed
command_sum = f"awk -v OFS='\\t' '{{print $1,$2,$3,$4+$5}}' {prefix3}_uniq > {prefix3}_uniq_updown.bed"
result = subprocess.run(command_sum, shell=True, text=True, check=True)


print(f"BigWig files generated: {prefix}_bw_up_sort.bw and {prefix}_bw_down_sort.bw")


### zscore normalization
# Define z-score window size
window_zscore = resolution_number * 500

# File paths
hek_bw_path = f"{prefix3}_uniq_updown.bed"
chr_sizes_path = "hg38.clean.chrom.sizes"
output_file = f"{prefix}_zscore_window_{window_zscore}.bed"
output_file1 = f"{prefix}_zscore_window_{window_zscore}_4col.bed"
output_bw = f"{prefix}_zscore_window_{window_zscore}.bw"

# Load data
hek_bw = pd.read_csv(hek_bw_path, sep="\t", names=["chr", "start", "end", "score"])
chr_sizes = pd.read_csv(chr_sizes_path, sep="\t", names=["chr", "size"])
chr_sizes = dict(zip(chr_sizes["chr"], chr_sizes["size"]))  # Convert to dictionary

# Define z-score window and calculation logic
def calculate_window(row, bw_data, chr_size, window_size):
    peak_center = (row["start"] + row["end"]) // 2
    window_start = max(0, peak_center - window_size // 2)
    window_end = min(chr_size, peak_center + window_size // 2)

    # Extract scores within the window
    window_scores = bw_data[(bw_data["start"] < window_end) & (bw_data["end"] > window_start)]["score"]
    window_mean = window_scores.mean()
    window_std = window_scores.std()
    new_score = (row["score"] - window_mean) / window_std if window_std > 0 else 0

    return {
        "chr": row["chr"],
        "start": row["start"],
        "end": row["end"],
        "new_score": new_score,
        "window_start": window_start,
        "window_end": window_end,
        "peak_score": row["score"],
        "window_mean": window_mean,
        "window_std": window_std,
    }

# Parallel processing by chromosome
def process_chromosome(chromosome):
    chr_size = chr_sizes[chromosome]
    bw_data_chr = hek_bw[hek_bw["chr"] == chromosome]
    results = [calculate_window(row, bw_data_chr, chr_size, window_zscore) for _, row in bw_data_chr.iterrows()]
    return results

# Run parallel jobs
n_jobs = min(cpu_count(), 8)
chromosomes = hek_bw["chr"].unique()
print(f"Processing chromosomes using {n_jobs} CPUs...")

with Pool(processes=n_jobs) as pool:
    all_results = pool.map(process_chromosome, chromosomes)

# Flatten and save
all_results_flat = [item for sublist in all_results for item in sublist]
all_results_df = pd.DataFrame(all_results_flat)
all_results_df.to_csv(output_file, sep="\t", index=False, header=False)

print(f"Z-score BED file saved to: {output_file}")

command9 = f"awk -v OFS='\\t' '{{print $1,$2,$3,$4}}' {output_file} > {output_file1}"
command10 = f"bedGraphToBigWig {output_file1} {chrom_sizes} {output_bw}"
result = subprocess.run(command9, shell=True, text=True, check=True)
result = subprocess.run(command10, shell=True, text=True, check=True)



# Print out the running time
end_time = time.time()
print(f"Script running time: {end_time - start_time} seconds")
