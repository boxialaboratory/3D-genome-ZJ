ChIP-seq analysis

This directory contains scripts used for ChIP-seq preprocessing, spike-in normalization, peak annotation, and feature embedding, as used in the study.

The workflow is organized in a fixed execution order, indicated by the numeric prefixes of each script.

Workflow overview

Annotate ChIP-seq peaks by genomic context

Embed multi-factor ChIP-seq signals using UMAP

Compute and visualize differential ChIP-seq signal matrices

Map reads to combined hg38 + dm6 genome for spike-in normalization

Downsample hg38 reads based on dm6 spike-in counts (scaling factors computed externally)

Scripts
1.annotated_peak_distribution_barplot.py

Purpose
Annotate ChIP-seq peaks by genomic category and visualize their distribution.

Input

Annotated BED file containing overlapping genomic features

Peak BED file (e.g. filtered ChIP-seq peaks)

Method

Assigns each peak to a single genomic category using a priority rule
(CTCF > TSS > enhancer > CDS > UTR > intron > intergenic)

Generates a single stacked horizontal bar plot showing peak distribution

Output

region_overlap_stacked_bar_ordered.pdf

2.UMAP.py

Purpose
Embed ChIP-seq signal patterns across multiple factors into a low-dimensional space.

Input

Merged peak BED file

Multiple ChIP-seq / ATAC-seq / insulation bigWig tracks

Method

Extracts mean signal per peak from each bigWig file

Log-transform and z-score normalize signals

Performs PCA, neighbor graph construction, UMAP embedding, and Leiden clustering

Output

UMAP plots colored by cluster or ChIP-seq signal

Table of UMAP coordinates and cluster labels
(*_adata_cell_cluster_umap.bed)

3.draw_diff.py

Purpose
Compute and visualize differential ChIP-seq signal matrices between two conditions.

Input

Two gzipped matrix files (e.g. WT vs KO), sharing region identifiers

Method

Aligns matrices by genomic coordinates

Computes log₂(KO / WT) difference

Uses multi-threading for matrix computation

Plots a rasterized heatmap optimized for large PDFs

Output

Differential heatmap PDF

4.qchip_map.sh

Purpose
Map ChIP-seq reads to a combined hg38 + dm6 reference genome for spike-in normalization.

Input

Paired-end FASTQ files

Method

Adapter trimming using Trim Galore

Alignment to merged hg38–dm6 genome using Bowtie2

Removal of duplicates

Separation of uniquely mapped hg38 and dm6 reads

Counting mapped reads per BAM file

Output

hg38- and dm6-specific BAM files

mapped_reads_summary.csv containing total mapped reads

5.qchip_downsample.py

Purpose
Downsample hg38-aligned ChIP-seq reads using precomputed scaling factors and generate bigWig tracks.

Important note
Downsampling fractions are not computed in this script.
They are calculated externally (e.g. in Excel) according to the spike-in normalization formula described in the Methods section.

Spike-in normalization formula

For a pairwise comparison (A vs B):

Downsampled hg38 reads
=
min
⁡
(
dm6
𝐴
,
dm6
𝐵
)
dm6
sample
×
hg38
sample
Downsampled hg38 reads=
dm6
sample
	​

min(dm6
A
	​

,dm6
B
	​

)
	​

×hg38
sample
	​


where:

dm6 = number of dm6-aligned reads

hg38 = number of hg38-aligned reads

The resulting fraction is then manually inserted into this script.

Method

Downsamples BAM files using sambamba view

Indexes downsampled BAMs

Generates CPM-normalized bigWig files using bamCoverage

Output

Downsampled BAM files

hg38 bigWig tracks for downstream analysis
