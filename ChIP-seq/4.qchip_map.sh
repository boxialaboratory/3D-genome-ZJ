### map sequencing reads to combined human hg38 and drosophila dm6 genome.
#!/bin/bash
if [[ ! -d "1.trim" ]]; then
        mkdir 1.trim
fi

if [[ ! -d "2.mapping" ]]
then
                mkdir 2.mapping
fi
if [[ ! -d "3.track_peak" ]]
then
                mkdir 3.track_peak
fi
# Set the reference genome index path
GENOME_INDEX="/home/gcpuser/sky_workdir/homo_droso/hg38_dm6_merged"

# Create a temporary directory for sorting
mkdir -p ./tmp

for file in *_R1_001.fastq.gz; do
    # Extract the prefix before "_R1_001.fastq.gz"
    sample="${file%%_R1_001.fastq.gz}"

    # Print the sample name (prefix)
    echo "Sample: $sample"

    # Optionally, construct the R1 and R2 file names using the sample prefix
    r1_file="./1.trim/${sample}_R1_001_val_1.fq.gz"
    r2_file="./1.trim/${sample}_R2_001_val_2.fq.gz"
    trim_galore -q 25 --paired --phred33 -e 0.1 --fastqc --clip_R1 10 --clip_R2 10 --length 36 --stringency 3 -j 48 -o 1.trim ${sample}_R1_001.fastq.gz ${sample}_R2_001.fastq.gz
    # Run Bowtie2 alignment
    bowtie2 -p 56 --no-mixed --no-discordant -x $GENOME_INDEX -1 "$r1_file" -2 "$r2_file" | \
        grep -v XS: - | samtools view -@ 56 -bhS -F4 - > "./2.mapping/${sample}_UniqMapped.bam"

    # Sort BAM file
    sambamba sort --tmpdir ./tmp/ -t 56 -o "./2.mapping/${sample}_UniqMapped_sort.bam" "./2.mapping/${sample}_UniqMapped.bam"

    # Remove duplicates
    sambamba markdup --tmpdir ./tmp/ -r -t 56 "./2.mapping/${sample}_UniqMapped_sort.bam" "./2.mapping/${sample}_UniqMapped_sort_rmdup.bam"

    # Extract uniquely mapped reads to hg38
    samtools view -@ 56 -h "./2.mapping/${sample}_UniqMapped_sort_rmdup.bam" | grep -v dm6 | sed 's/hg38_chr/chr/g' | \
        samtools view -@ 56 -bhS - > "./2.mapping/${sample}_hg38.UniqMapped_sorted_rmdup.bam"
  
    # Extract uniquely mapped reads to dm6
    samtools view -@ 56 -h "./2.mapping/${sample}_UniqMapped_sort_rmdup.bam" | grep -v hg38 | sed 's/dm6_chr/chr/g' | \
        samtools view -@ 56 -bhS - > "./2.mapping/${sample}_dm6.UniqMapped_sorted_rmdup.bam"

    # Remove intermediate files
    rm -f "./2.mapping/${sample}_UniqMapped.bam" "./2.mapping/${sample}_UniqMapped_sort.bam" "./2.mapping/${sample}_UniqMapped_sort_rmdup.bam"
    rm -f "$r1_file" "$r2_file"
done
        

### get the read number in each file.

# Input directory containing BAM files
bam_dir="./"
output_csv="mapped_reads_summary.csv"

# Create or overwrite the CSV header
echo "BAM File,Total Mapped Reads" > $output_csv

# Loop through each BAM file in the directory
for bam_file in "$bam_dir"/*.bam; do
  # Get the prefix of the BAM file
  prefix=$(basename "$bam_file" .UniqMapped_sorted_rmdup.bam)

  # Generate index using samtools
  samtools index -@ 48 "$bam_file"

  # Run idxstats and calculate total mapped reads
  samtools idxstats "$bam_file" > "${prefix}_idxstats.txt"
  total_reads=$(awk '{sum += $3} END {print sum}' "${prefix}_idxstats.txt")

  # Append the result to the CSV file
  echo "${prefix},${total_reads}" >> $output_csv

  # Optional: Clean up the idxstats file if no longer needed
  # rm "${prefix}_idxstats.txt"
done

echo "Processing complete. Results saved to $output_csv"
