ChIP-seq

This directory contains scripts used for ChIP-seq preprocessing, quantitative spike-in normalization, peak annotation, and feature embedding, as used in this study.

Scripts are intended to be executed in numeric order.

Requirements

Core dependencies used by the scripts in this directory:

Trim Galore ≥ 0.6

Bowtie2 ≥ 2.4

samtools ≥ 1.13

sambamba ≥ 0.8

deepTools ≥ 3.5

MACS2 ≥ 2.2

bedtools ≥ 2.27

Python ≥ 3.8

pandas

numpy

matplotlib

seaborn

pyBigWig

scanpy / anndata

Workflow overview

Peak genomic annotation and summary

Multi-factor ChIP-seq signal embedding (UMAP)

Differential signal heatmap visualization

Spike-in aware read mapping (hg38 + dm6)

Spike-in–normalized downsampling and track generation

Scripts
1.annotated_peak_distribution_barplot.py

Annotates ChIP-seq peaks by genomic context and visualizes their distribution.

Input

Annotated BED file (feature overlaps)

Peak BED file

Method

Each peak is assigned to a single category using a priority rule:
CTCF > TSS > Enhancer > CDS > 5′UTR > 3′UTR > Intron > Intergenic

Outputs a stacked horizontal bar plot

Output

region_overlap_stacked_bar_ordered.pdf

2.UMAP.py

Embeds ChIP-seq and ATAC-seq signal patterns across merged peaks.

Input

Merged peak BED file

Multiple ChIP-seq / ATAC-seq / insulation bigWig tracks

Method

Mean signal extraction per peak (pyBigWig)

log1p + z-score normalization

PCA → UMAP → Leiden clustering

Output

UMAP plots colored by factor enrichment

UMAP coordinates and cluster labels
(*_adata_cell_cluster_umap.bed)

3.draw_diff.py

Computes and visualizes differential ChIP-seq signal matrices between two conditions.

Input

Two gzipped matrices with shared region identifiers (e.g. WT vs KO)

Method

Aligns matrices by genomic coordinates

Computes log₂(KO / WT)

Multi-threaded computation

Rasterized heatmap rendering for large matrices

Output

Differential heatmap PDF

4.qchip_map.sh

Maps ChIP-seq reads to a concatenated hg38 + dm6 genome for spike-in normalization.

Method

Adapter trimming (Trim Galore)

Alignment to merged genome (Bowtie2)

Removal of secondary alignments and duplicates

Separation of hg38 and dm6 reads

Read counting per BAM file

Output

hg38- and dm6-specific BAM files

mapped_reads_summary.csv

5.qchip_downsample.py

Downsamples hg38-aligned reads using precomputed spike-in scaling factors and generates coverage tracks.

Important
Downsampling fractions are computed externally (e.g. Excel), using dm6 spike-in read counts as described in the Methods.

Method

Random downsampling with sambamba view

BAM indexing

bigWig generation using bamCoverage

No CPM/RPKM normalization applied

Output

Downsampled BAM files

hg38 bigWig tracks for quantitative analysis
