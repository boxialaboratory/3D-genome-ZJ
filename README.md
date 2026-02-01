# 3D-genome-ZJ

This repository contains analysis scripts used for **ChIP-seq and Micro-C data processing and visualization** in this study.

The code is organized by data modality and follows a **numeric execution order** within each directory.

---


## ChIP-seq

Scripts for ChIP-seq preprocessing, quantitative spike-in normalization, peak annotation, and feature embedding.

---

### Requirements

**Core command-line tools**

- Trim Galore ≥ 0.6  
- Bowtie2 ≥ 2.4  
- samtools ≥ 1.13  
- sambamba ≥ 0.8  
- deepTools ≥ 3.5  
- MACS2 ≥ 2.2  
- bedtools ≥ 2.27  

**Python (≥ 3.8)**

- pandas  
- numpy  
- matplotlib  
- seaborn  
- pyBigWig  
- scanpy  
- anndata  

---

### Workflow overview

1. Annotate ChIP-seq peaks by genomic context  
2. Embed, cluster and visualize multi-factor ChIP-seq signals using UMAP  
3. Compute differential ChIP-seq signal matrices  
4. Map reads to hg38 + dm6 genome for spike-in normalization  
5. Downsample hg38 reads based on dm6 spike-in levels  

---

### Scripts

#### 1.annotated_peak_distribution_barplot.py

Annotates ChIP-seq peaks by genomic feature and visualizes their distribution.

- Assigns each peak to a single category using a priority rule  
  (CTCF > TSS > Enhancer > CDS > 5′UTR > 3′UTR > Intron > Intergenic)
- Outputs a stacked horizontal bar plot

**Output**
- `region_overlap_stacked_bar_ordered.pdf`

---

#### 2.UMAP.py

Embeds ChIP-seq and ATAC-seq signal patterns across merged peaks.

- Extracts mean signal per peak from bigWig files
- Applies log-transformation and z-score normalization
- Performs PCA, UMAP embedding, and Leiden clustering

**Output**
- UMAP plots
- UMAP coordinates and cluster labels  
  (`*_adata_cell_cluster_umap.bed`)

---

#### 3.draw_diff.py

Computes and visualizes differential ChIP-seq signal matrices between two conditions.

- Aligns matrices by genomic coordinates
- Computes log₂(KO / WT)
- Multi-threaded computation
- Rasterized heatmap output for large matrices

**Output**
- Differential heatmap PDF

---

#### 4.qchip_map.sh

Maps ChIP-seq reads to a **concatenated hg38 + dm6 genome** for spike-in normalization.

- Adapter trimming
- Alignment to merged genome
- Duplicate removal
- Separation of hg38 and dm6 reads
- Read counting per BAM file

**Output**
- hg38- and dm6-specific BAM files  
- `mapped_reads_summary.csv`

---

#### 5.qchip_downsample.py

Downsamples hg38-aligned reads using externally computed spike-in normalization factors and generates coverage tracks.

Spike-in downsampling calculation

For each pairwise comparison (e.g. sample A vs sample B), Drosophila dm6–aligned reads were used as an internal reference.
The sample with the lower dm6 spike-in count was used as the normalization baseline, and human (hg38) reads in the other sample were downsampled accordingly.

Downsampling formula:

hg38_downsampled =
min(dm6_A, dm6_B) / dm6_sample × hg38_sample


where:

dm6_A, dm6_B are the numbers of dm6-aligned reads in samples A and B

dm6_sample is the number of dm6-aligned reads in the sample being scaled

hg38_sample is the number of hg38-aligned reads in the corresponding sample

The resulting scaling factor is calculated externally (e.g. in Excel) and manually provided to this script.

**Notes**

- Random downsampling is performed using `sambamba view`
- No CPM/RPKM normalization is applied during track generation

---
## Micro-C

Scripts for Micro-C data processing, stripe detection, loop calling, differential analysis, and integrated visualization with ChIP-seq tracks.

Scripts are intended to be executed in numeric order.  
They are modular analysis utilities rather than a fully automated pipeline.

---

### Requirements

**Core command-line tools**

- cooler ≥ 0.9  
- bedtools ≥ 2.27  
- bedGraphToBigWig  
- mustache  
- coolpuppy  

**Python (≥ 3.8)**

- numpy  
- pandas  
- scipy  
- matplotlib  
- seaborn  
- cooler  
- cooltools  
- bioframe  
- straw  
- pyranges  

---

### Workflow overview

1. Summarize Micro-C mapping statistics  
2. Identify and remove low-coverage or noisy genomic bins  
3. Detect stripe features using z-score–normalized interaction profiles  
4. Visualize Micro-C contact maps and differential matrices  
5. Quantify aggregate peak analysis (APA) for ChIP-seq–defined anchors  
6. Visualize Micro-C contact maps with overlaid ChIP-seq tracks  
7. Call and merge chromatin loops across resolutions  
8. Quantify and classify loop strength changes between conditions  

---

### Scripts

#### 1.get_stat.py

Summarizes Micro-C mapping statistics from `.stat` files.

- Extracts total reads, duplicate reads, and deduplicated reads  
- Computes cis/trans interaction ratios  
- Outputs per-sample statistics in column-wise format

**Output**
- Printed summary tables for quality control and comparison

---

#### 2.remove_unmapped.py

Identifies low-coverage or noisy genomic regions from Micro-C contact maps.

- Computes per-bin interaction coverage  
- Flags bins with extremely low coverage or excessive noise  
- Merges adjacent low-quality bins into continuous genomic intervals

**Output**
- `genome_low_coverage_noisy_regions.bed`

---

#### 3.stripe_calling_zscore_norm.py

Detects stripe features along the diagonal of balanced Micro-C contact maps.

- Computes directional interaction sums around each genomic bin  
- Identifies stripe peaks based on prominence thresholds  
- Separates positive and negative stripe signals  
- Generates z-score–normalized stripe tracks using local genomic windows

**Output**
- Stripe peak BED files  
- Strand-specific stripe BigWig tracks  
- Z-score–normalized stripe BigWig files

---

#### 4.plot_microc_diff_map.py

Visualizes Micro-C contact maps and differential interaction matrices between conditions.

- Extracts balanced contact matrices from `.mcool` files  
- Plots WT, KO, and KO–WT difference maps in triangular format  
- Performs statistical comparisons (Wilcoxon and Fisher’s exact tests)

**Output**
- Contact matrix PDFs  
- Differential contact map PDFs

---

#### 5.APA_for_pairwise_chip_peak.py

Performs aggregate peak analysis (APA) for pairs of ChIP-seq–defined genomic anchors.

- Generates BEDPE interaction files from ChIP-seq peak sets  
- Computes APA profiles using `coolpuppy`  
- Automates pairwise combinations via a generated shell script

**Output**
- BEDPE interaction files  
- APA matrices and summary text files  
- APA visualization PDFs (via `plotpup.py`)

---

#### 6.microC_visualization_with_ChIP_tracks.py

Generates integrated visualizations of Micro-C contact maps and ChIP-seq signal tracks.

- Plots Micro-C heatmaps at specified resolution  
- Overlays multiple ChIP-seq BigWig tracks  
- Annotates genomic features and motif orientations  
- Supports strand-aware and condition-aware visualization

**Output**
- Multi-track Micro-C + ChIP-seq visualization PDFs

---

#### 7.call_merge_loop.sh

Calls chromatin loops at multiple resolutions and merges results.

- Loop detection performed using `mustache`  
- Loops from different resolutions are merged  
- Final loop set is sorted and indexed

**Output**
- `hic_loop_all_sorted.txt`

---

#### 8.plot_up_down_loops.py

Quantifies loop strength changes between conditions and classifies loop categories.

- Extracts loop interaction strengths from Hi-C / Micro-C data  
- Computes log₂(WT) vs log₂(KO) interaction differences  
- Classifies loops by genomic annotation (CTCF, promoter, enhancer)  
- Identifies upregulated and downregulated loops  
- Exports loop sets and visualization plots

**Output**
- Loop strength comparison plots  
- Upregulated and downregulated loop BEDPE files  
- Annotated loop tables

---

## RNA-seq

Bulk RNA-seq analysis associated with this study was performed following the workflow described in the repository below:

https://github.com/boxialaboratory/Tail-Loss-Shawn-RNA-seq/tree/main/05_bulk_RNA_seq

---

## Citation

If you use this code, please cite the associated manuscript:

> Bai J., Lyu Q., Tan J., Issner R., Tischer B., Liu H., Goel V., Ling X., Bernstein B.E., Tsirigos A., Hansen A.S., Xia B.  
> **High-throughput in silico screen uncovers key regulators of 3D genome architecture**  
> *bioRxiv* 2025.12.09.693120  
> https://doi.org/10.64898/2025.12.09.693120
