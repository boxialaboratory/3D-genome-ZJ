import os
import subprocess

# Input data: filenames and downsampling values
input_data = [
    ("contrl1", 0.7850275072),
    ("sample1", 1),
    ("control2", 1),
    ("sample2", 0.7925966902),
    ("control3",1),
    ("sample3",0.848749704)
]

# Directory for the downsampled files and bigWig outputs
output_dir = "downsample"
bigwig_dir = "bigwig"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(bigwig_dir, exist_ok=True)

# Suffix for BAM files
suffix = ".UniqMapped_sorted_rmdup.bam"

# Loop through the input data
for filename, fraction in input_data:
    bam_file = f"{filename}{suffix}"
    output_bam = os.path.join(output_dir, f"{filename}_downsampled.bam")
    output_bw = os.path.join(bigwig_dir, f"{filename}_hg38_sort_cpm.rmdup.bw")

    if fraction == 1:
        # If the fraction is 1, just copy the BAM file
        print(f"Copying {bam_file} to {output_bam}")
        subprocess.run(["cp", bam_file, output_bam])
    else:
        # Otherwise, downsample the BAM file
        print(f"Downsampling {bam_file} with fraction {fraction}")
        subprocess.run([
            "sambamba", "view", "-h", "-t", "56", "-f", "bam",
            "--subsampling-seed=123", "-s", str(fraction), bam_file,
            "-o", output_bam
        ])

    # Generate index for the downsampled BAM file
    print(f"Indexing {output_bam}")
    subprocess.run(["sambamba", "index", "-t", "56", output_bam])

    # Run bamCoverage
    print(f"Running bamCoverage for {output_bam}")
    subprocess.run([
        "bamCoverage", "--bam", output_bam,
        "-o", output_bw,
        "--binSize", "25",
        "--smoothLength", "75",
        "--numberOfProcessors", "48",
        "--centerReads",
        "--extendReads"
    ])

print("Processing complete. Downsampled BAM and bigWig files are in the 'downsample' and 'bigwig' directories.")
