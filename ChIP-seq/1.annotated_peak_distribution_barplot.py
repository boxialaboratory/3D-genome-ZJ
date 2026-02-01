import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from collections import Counter

# Configure matplotlib to use editable fonts in PDFs
mpl.rcParams['pdf.fonttype'] = 42

# Define the order of regions with priority
region_priority = {
    "ctcf": 1,
    "tss": 2,
    "enhancer": 3,
    "cds": 4,
    "utr5": 5,
    "utr3": 6,
    "intron": 7,
    "intergenic": 8  # Intergenic as the lowest priority
}

# Define colors for each region type
region_colors = {
    "ctcf": "#0485d1",   # CTCF: Blue
    "tss": "#8e82fe",    # TSS: Purple
    "enhancer": "#960056",  # Enhancer: Magenta
    "cds": "#ff7f0e",    # CDS: Orange
    "utr5": "#2ca02c",   # 5'UTR: Green
    "utr3": "#d62728",   # 3'UTR: Red
    "intron": "#9467bd", # Intron: Violet
    "intergenic": "#999999" # Intergenic: Grey
}

# Function to get the region type from the 13th column
def get_region_type(annotation):
    annotation = annotation.lower()
    if "ctcf" in annotation:
        return "ctcf"
    elif "tss" in annotation:
        return "tss"
    elif "enhancer" in annotation or "stitch" in annotation or "super" in annotation:
        return "enhancer"
    elif "cds" in annotation:
        return "cds"
    elif "utr5" in annotation or "_utr5_" in annotation:
        return "utr5"
    elif "utr3" in annotation or "_utr3_" in annotation:
        return "utr3"
    elif "intron" in annotation:
        return "intron"
    else:
        return "intergenic"

# Load the BED files into DataFrames
file_path_aa = "./aa.bed"  # Replace with the path to the annotated regions file
file_path_peaks = "./HEK_H7J_noblack.bed"  # Replace with the path to your peak file

df = pd.read_csv(file_path_aa, sep='\t', header=None)
df_peaks = pd.read_csv(file_path_peaks, sep='\t', header=None)

# Extract the peak and annotated region information
df['peak'] = df[3]  # Assuming the peak is in the 4th column (index 3)
df['region_type'] = df[13].apply(get_region_type)  # Annotated region is in the 14th column (index 13)

# Get the highest priority region type for each unique peak
peak_region_priority = (
    df.groupby('peak')['region_type']
    .apply(lambda x: min(x, key=lambda region: region_priority.get(region, float('inf'))))
)

# Count the occurrences of each region type
region_counts = Counter(peak_region_priority)

# Calculate the total number of peaks from HEK_H7J_noblack.bed
total_peaks = df_peaks.shape[0]
print(f"Total number of peaks in HEK_H7J_noblack.bed: {total_peaks}")

# Calculate the number of intergenic peaks
intergenic_count = total_peaks - len(peak_region_priority)
print(f"Number of intergenic peaks: {intergenic_count}")

# Add the intergenic count to the region counts
region_counts["intergenic"] = intergenic_count

# Sort the region counts based on the defined priority
sorted_regions = sorted(region_counts.items(), key=lambda x: region_priority.get(x[0], float('inf')))
labels, sizes = zip(*sorted_regions)  # Unpack sorted regions into labels and sizes
colors = [region_colors[label] for label in labels]  # Use the color mapping for sections

# Plot the stacked bar chart
fig, ax = plt.subplots(figsize=(10, 2))  # Adjust the height to make it look like a single bar

left = 0  # Track the starting point for each section in the bar
for size, color, label in zip(sizes, colors, labels):
    ax.barh(
        "Annotated Regions", size, left=left, color=color, edgecolor='black',
        label=f"{label} ({size})"
    )
    left += size  # Update the starting point for the next section

# Add the total number of peaks as the title
plt.title(f"Distribution of Annotated Regions Overlapping Peaks\nTotal Peaks: {total_peaks}", fontsize=14)

# Add a legend below the plot
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=len(labels), fontsize=10, frameon=False)

# Remove axes and grid lines for a cleaner look
ax.axis('off')

# Save the bar chart as a PDF
plt.savefig("region_overlap_stacked_bar_ordered.pdf", format='pdf', bbox_inches='tight')

# Display the chart
plt.show()
