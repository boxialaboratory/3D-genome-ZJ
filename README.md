# 3D-genome-ZJ

This repository contains analysis scripts used for **ChIP-seq and Micro-C data processing and visualization** in this study.

The code is organized by data modality and follows a **numeric execution order** within each directory.

---

## Repository structure

3D-genome-ZJ/
├── ChIP-seq/
│ ├── 1.annotated_peak_distribution_barplot.py
│ ├── 2.UMAP.py
│ ├── 3.draw_diff.py
│ ├── 4.qchip_map.sh
│ └── 5.qchip_downsample.py
├── MicroC/
│ ├── 1.get_stat.py
│ ├── 2.remove_unmapped.py
│ ├── 3.stripe_calling_zscore_norm.py
│ ├── 4.plot_microc_diff_map.py
│ ├── 5.APA_for_pairwise_chip_peak.py
│ ├── 6.microC_visualization_with_ChIP_tracks.py
│ ├── 8.call_merge_loop.sh
│ └── 9.plot_up_down_loops.py
└── README.md


---

## General notes

- Scripts are **not a one-click pipeline**
- Paths and parameters are expected to be adjusted to local environments
- Numeric prefixes indicate **recommended execution order**
- Scripts were written for **analysis reproducibility**, not general-purpose packaging

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
2. Embed multi-factor ChIP-seq signals using UMAP  
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

##### Spike-in downsampling calculation

For each pairwise comparison (e.g. sample A vs sample B), Drosophila dm6–aligned reads were used as an internal reference.  
The sample with the lower dm6 spike-in count was used as the normalization baseline, and human (hg38) reads in the other sample were downsampled accordingly.

\[
\mathrm{hg38}_{\mathrm{downsampled}} =
\frac{\min(\mathrm{dm6}_A,\ \mathrm{dm6}_B)}{\mathrm{dm6}_{\mathrm{sample}}}
\times \mathrm{hg38}_{\mathrm{sample}}
\]

where:

- \(\mathrm{dm6}_A\), \(\mathrm{dm6}_B\) are the numbers of dm6-aligned reads in samples A and B  
- \(\mathrm{hg38}_{\mathrm{sample}}\) is the number of hg38-aligned reads in the corresponding sample  

The resulting scaling factor is calculated externally (e.g. in Excel) and manually provided to this script.

**Notes**

- Random downsampling is performed using `sambamba view`
- No CPM/RPKM normalization is applied during track generation

---
