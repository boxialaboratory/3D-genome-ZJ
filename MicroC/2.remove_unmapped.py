import cooler
import numpy as np
import pandas as pd

# Load the cool file at 2 kb resolution
cool_path = ".mcool::/resolutions/5000"
clr = cooler.Cooler(cool_path)

# Define thresholds based on your criteria
empty_threshold = 2              # Threshold for truly empty bins
noise_lower_threshold = 10       # Lower bound of noisy region
noise_upper_threshold = 50       # Upper bound of noisy region

# Define a reasonable chunk size to avoid excessive noise
chunk_size_bp = 120_000_000  # Larger chunk size to reduce noise in filtering

# Process each chromosome
for chrom in clr.chromnames:
    print(f"Processing {chrom}")
    bins = clr.bins().fetch(chrom)
    chrom_length = bins['end'].iloc[-1]  # Length of the chromosome
    num_chunks = (chrom_length // chunk_size_bp) + 1

    # Create a temporary file to store results for this chromosome
    with open(f"{chrom}_low_coverage_noisy_regions.bed", "w") as chrom_outfile:
        # Process each chunk
        for chunk in range(num_chunks):
            start_bp = chunk * chunk_size_bp
            end_bp = min((chunk + 1) * chunk_size_bp, chrom_length)
            
            # Select bins within the current chunk
            bins_chunk = bins[(bins['start'] >= start_bp) & (bins['end'] <= end_bp)]
            matrix_chunk = clr.matrix(balance=False).fetch((chrom, start_bp, end_bp))
            
            # Calculate coverage for the chunk
            coverage_chunk = np.array(matrix_chunk.sum(axis=1)).flatten()
            bins_chunk = bins_chunk.copy()
            bins_chunk['coverage'] = coverage_chunk

            # Apply the filters for empty and noisy bins
            filtered_bins_chunk = bins_chunk[
                (bins_chunk['coverage'] < empty_threshold) | 
                ((bins_chunk['coverage'] >= noise_lower_threshold) & (bins_chunk['coverage'] <= noise_upper_threshold))
            ]

            # Write the filtered bins to the chromosome-specific output file
            filtered_bins_chunk[['chrom', 'start', 'end']].to_csv(chrom_outfile, sep="\t", index=False, header=False)

# Combine all chromosome BED files into a single BED file with post-processing to merge adjacent regions
with open("genome_low_coverage_noisy_regions.bed", "w") as outfile:
    for chrom in clr.chromnames:
        df = pd.read_csv(f"{chrom}_low_coverage_noisy_regions.bed", sep="\t", header=None, names=["chrom", "start", "end"])
        
        # Merge adjacent low-coverage regions
        merged_df = df.sort_values(by="start").reset_index(drop=True)
        merged_intervals = []
        current_start, current_end = merged_df.iloc[0]['start'], merged_df.iloc[0]['end']
        
        for i in range(1, len(merged_df)):
            start, end = merged_df.iloc[i]['start'], merged_df.iloc[i]['end']
            if start <= current_end:  # Adjacent or overlapping interval
                current_end = max(current_end, end)
            else:
                merged_intervals.append([chrom, current_start, current_end])
                current_start, current_end = start, end
        merged_intervals.append([chrom, current_start, current_end])  # Add last interval
        
        # Write merged intervals to final output
        for interval in merged_intervals:
            outfile.write("\t".join(map(str, interval)) + "\n")
