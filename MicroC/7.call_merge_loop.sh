## call loop at different resolution
python3 ~/mustache.py -f AA.mcool -r 10kb -pt 0.1 -st 0.88 -o hic_10k_res -p 8 --octaves 2 

### merge loops at different resolutions
cat hic_10k_res.txt hic_5k_res.txt hic_2k_res.txt hic_1k_res.txt >> merge_loop_10k_5k_2k_1k.bed
less merge_loop_10k_5k_2k_1k.bed|awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,"loop_"NR}' > hic_loop_all.txt
bedtools sort -i hic_loop_all.txt > hic_loop_all_sorted.txt


###select the finest loop.py
import pandas as pd
from intervaltree import IntervalTree

# Load the BEDPE file
input_file = "hic_wt_out_10k_5k_2k_1k_real_0.1_0.88_all_sorted.txt"
columns = ['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'loop_id']
df = pd.read_csv(input_file, sep='\s+', header=None, names=columns)

# Calculate loop size for prioritization
df['size'] = (df['end1'] - df['start1']) + (df['end2'] - df['start2'])

# Initialize IntervalTrees for both ends
tree_start = IntervalTree()
tree_end = IntervalTree()

# Add loops to the IntervalTrees based on both ends
for _, row in df.iterrows():
    tree_start[row['start1']:row['end1']] = row.to_dict()
    tree_end[row['start2']:row['end2']] = row.to_dict()

# Function to merge loops, keeping only the smallest when both ends overlap
def merge_loops():
    merged_loops = []
    visited = set()  # Track visited loops to avoid duplication

    # Sort by size to prioritize smaller loops
    intervals = sorted(tree_start, key=lambda x: x.data['size'])

    for interval in intervals:
        loop_data = interval.data
        loop_id = loop_data['loop_id']

        if loop_id in visited:
            continue  # Skip if already processed

        # Find overlaps at both ends
        overlaps_start = tree_start.overlap(interval.begin, interval.end)
        overlaps_end = tree_end.overlap(loop_data['start2'], loop_data['end2'])

        # Collect loops that overlap at both ends
        overlapping_loops = [
            o.data for o in overlaps_start if o.data['loop_id'] in 
            [oe.data['loop_id'] for oe in overlaps_end]
        ]

        if overlapping_loops:
            # Keep only the smallest loop among candidates
            smallest_loop = min(overlapping_loops, key=lambda x: x['size'])

            # Mark all overlapping loops as visited
            for overlap in overlapping_loops:
                visited.add(overlap['loop_id'])

            # Add the smallest loop to the final result
            merged_loops.append(smallest_loop)
        else:
            # If only one end overlaps, retain the current loop
            merged_loops.append(loop_data)

    return pd.DataFrame(merged_loops)

# Apply the merging logic
merged_df = merge_loops()

# Save the result to a new file
output_file = "merged_loops_corrected_final.bed1"
merged_df.to_csv(output_file, sep='\t', header=False, index=False)

print(f"Merged loops saved to {output_file}")
