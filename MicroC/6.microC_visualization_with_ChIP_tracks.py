# ---------- Standard library ----------
import os
import h5py

# ---------- Scientific computing ----------
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, fisher_exact

# ---------- Plotting ----------
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

# ---------- Hi-C & genomics ----------
import cooler
import cooltools
from cooltools import expected_cis, expected_trans, insulation
from cooltools.lib import plotting
import bioframe

# ---------- CoolBox ----------
import coolbox
from coolbox.api import *
from coolbox.utilities import GenomeRange

# ---------- pyGenomeTracks ----------
from pygenometracks import readBed, readGtf
from pygenometracks.tracks import BedTrack


# Define plotting aesthetics as global variables
HIGHLIGHT_PARAMS = {"color":"gray",
                   "alpha":0.05,
                   "border_line":False}



SPACER_HEIGHT = 0.4
DEFAULT_COLOR_LIST = ["#FF0000", "#00FA00", "#0000FF", "#A0F000"]
PRO_SEQ_TRACK_HEIGHT = 1.2


ANNOTATIONS_PARAMS = [{"height":4},
                      {"file": "gencode_colored.bed12", # PATH TO BED FILE OF GENE ANNOTATIONS
                       "max_labels":1000,
                       "merge_transcripts":True, 
                        "gene_style": "flybase",
                       "itemRgb":True,
                       "gene_rows":6}]

def cat_prefix_to_list(prefix, file_list):
    """
    Helper function that concatenates a prefix to every element
    in a list. 
    """
    return [os.path.join(prefix, file) for file in file_list]

def make_region_plot(region, resolution, highlight_region_list, microc_file, microc_title,
                     chip_list_list, forward_proseq_list, reverse_proseq_list, condition_order, 
                     bed_files, bigwig_bins=3600, highlight_buffer=1000, pro_bins=None):
    
    """
    Function that makes a plot of a region with Micro-C, ChIP-seq.
    
    region: (str) string that specifies region for plotting.
    resolution: (int) resolution of micro-c map used
    highlight_region_list: List of regions of format [chr, start, end] to highlight
    microc_file: (str) file path to micro-c data
    microc_title: (str) title for micro-c
    chip_list_list: List of lists containing paths to bigwigs
    condition_order: List indicating order of conditions to plot. 
    """
    
    heatmap = Cool(microc_file, resolution=resolution, **MICROC_PLOTTING_PARAMS) + Title(microc_title)
    frame = heatmap
    if pro_bins is None:
        pro_bins = bigwig_bins
    def iterate_over_datasets(frame):
        i = 0 
        for chip_list in chip_list_list:
            bigwig_list = make_bigwig_list(chip_list, region, condition_order, 
                                           DEFAULT_COLOR_LIST[i%len(DEFAULT_COLOR_LIST)], bw_bins=bigwig_bins)
            i += 1
        
            for bigwig in bigwig_list:
                frame += bigwig + Spacer(SPACER_HEIGHT)
# ###pro-seq
#         f_bigwig_list = make_bigwig_list(forward_proseq_list, region, condition_order_pro, 
#                                        "red", bw_bins=pro_bins,  
#                                          track_height=PRO_SEQ_TRACK_HEIGHT)
#         r_bigwig_list = make_bigwig_list(reverse_proseq_list, region, condition_order_pro, len(reverse_proseq_list)*[""], 
#                                "blue", bw_bins=pro_bins, 
#                                          track_height=PRO_SEQ_TRACK_HEIGHT)
#         r_bigwig_list = [r_bigwig + Inverted() for r_bigwig in r_bigwig_list]
#         for f_bigwig, r_bigwig in zip(f_bigwig_list, r_bigwig_list):
#                 frame += f_bigwig + r_bigwig + Spacer(SPACER_HEIGHT)



        return frame
    
    # check if it needs a highlight
    if len(highlight_region_list) > 0:
        highlight_regions = []
        for highlight_region in highlight_region_list:
            highlight_region[1] -= highlight_buffer
            highlight_region[2] += highlight_buffer 
            highlight_regions.append(tuple(highlight_region))

        highlights = HighLights(highlight_regions, **HIGHLIGHT_PARAMS)
        with highlights:
            frame = iterate_over_datasets(frame)
    else:
        frame = iterate_over_datasets(frame)
   ### add bed motif 
    for bed_file in bed_files:
        bed_track = StrandMarkerTrack(bed_file)
        bed_track.fetch_data(region)
        frame += bed_track + Spacer(0.2)
    frame += XAxis()
    annotations = NewBed(*ANNOTATIONS_PARAMS)
    frame += annotations
    return frame

def make_bigwig_list(bigwigs, region, condition_order, color, bw_bins, track_height=3.5, autoscale=False, 
                     y_max=None,threshold_color="blue",threshold=0):
    
    """
    Helper function for `make_region_plot` that converts file paths to bigwig plotting objects.
    """
    assert len(bigwigs) % len(condition_order) == 0
    bigwig_list = []

    for bigwig, title, color, y_max, y_min in zip(bigwigs, condition_order, DEFAULT_COLOR_LIST, ymax, ymin):
        if title == "insulation":  # Adjust color for insulation track
            bigwig_track = BigWig(bigwig, number_of_bins=bw_bins, threshold_color="blue", threshold=0)
        elif title == "stripe":  # Adjust color for insulation track
            bigwig_track = BigWig(bigwig, number_of_bins=bw_bins, threshold_color="#ba76a5", threshold=1.5)
        else:
            bigwig_track = BigWig(bigwig, number_of_bins=bw_bins)
        
        track = bigwig_track + Title(title) + Color(color) + TrackHeight(track_height) + MaxValue(y_max) + MinValue(y_min)
        bigwig_list.append(track)
    if autoscale:
        bigwig_list = auto_scale_bigwigs(bigwig_list, region, y_max=y_max)
    return bigwig_list

def auto_scale_bigwigs(bigwig_list, region, y_max=None, y_min=0):
    """
    Autoscales y axis of bigwigs in bigwig list based off values in `region`.
    """
    def get_max_y_value(bigwig_list, region):
        max_y_values = []
        for bigwig in bigwig_list:
            max_y_values.append(np.amax(bigwig.fetch_plot_data(GenomeRange(region))))

        y_max = round(max(max_y_values) * 1.05)
        return y_max
    
    if y_max is None:
        if type(region) == list:
            y_max_list = [get_max_y_value(bigwig_list, single_region) for single_region in region]
            y_max = max(y_max_list)
        else:
            y_max = get_max_y_value(bigwig_list, region)
                
    return [bigwig + MaxValue(y_max) + MinValue(y_min) for bigwig in bigwig_list]

def read_cooler(cooler_path, resolution=250):
    """
    Simple function wrapper to read cooler file at a particular resolution.
    """
    if clr.fileops.is_cooler(cooler_path):
        return clr.Cooler(cooler_path)
    elif clr.fileops.is_multires_file(cooler_path):
        return clr.Cooler(cooler_path + "::resolutions/" + str(int(resolution)))

# define new bed format
class NewBed(Track):
    """
    Custom pygenometracks class that allows for control over the aesthetics of gene annotations in
    region plots. 
    """
    def __init__(self, coolbox_prop_dict, pygenometracks_prop_dict):
        super().__init__(coolbox_prop_dict)  # init 
        self.pygenometracks_object = BedTrack(pygenometracks_prop_dict)

    def fetch_data(self, gr, **kwargs):
        pass

    def plot(self, ax, gr, **kwargs):
#         x = gr.start + self.properties['offset'] * (gr.end - gr.start)
#         ax.text(x, 0, gr.chrom, fontsize=self.properties['fontsize'])
#         ax.set_xlim(gr.start, gr.end)
        self.pygenometracks_object.plot(ax, gr.chrom, gr.start, gr.end)
        ax.set_xlim([gr.start, gr.end])

        
        ##region_str = ("chr19:27,638,001-30,000,000")
# define StrandMarkerTrack for motif orientation
class StrandMarkerTrack(Track):
    def __init__(self, bed_file, **kwargs):
        super().__init__(**kwargs)
        self.bed_file = bed_file
        self.data = []

    def fetch_data(self, gr, **kwargs):
        """get region from bed"""
        self.data = []
        with open(self.bed_file, "r") as f:
            for line in f:
                fields = line.strip().split("\t")
                chrom, start, end, strand = fields[0], int(fields[1]), int(fields[2]), fields[3]
                if chrom == gr.chrom and start < gr.end and end > gr.start:
                    self.data.append((start, end, strand))
        return self.data

    def plot(self, ax, gr, **kwargs):
        """plot arrow"""
        for start, end, strand in self.data:
            mid = (start + end) / 2
            marker = ">" if strand == "+" else "<"
            ax.text(mid, 0, marker, ha="center", va="center", fontsize=10, color="blue")
        ax.set_xlim(gr.start, gr.end)
        ax.set_ylim(-1, 1)
        ax.axis("off") 


### use parameters(example):
mcool_files = ["hct_wt_merge_all_1.48B_UU_dedup.pairs_50bp.mcool"]
chip_condition_bigwigs = [
    "HEK293_J_HG38_SRR1015994.bigwig",
    "Hela_J_hg38_SRR1016019.bigwig",
"HEK293_CTCF_hg38_SRR299276.bigwig",
    "Hela_CTCF_hg38_SRR6338460.bigwig",
    "23022FL-03-01-28_S28_L001_hg38_sort_cpm.rmdup.bigwig",
    "Hela_RAD21_hg38_SRR6701633.bigwig",
        "H7_2min_6min_merge_sort_dedup_CPM.bigwig",
    "hek_znf_merge_sort_dedup.bigwig",
]

# ctcf motif bed files
bed_files = [
    "hek_CTCF_motifs_fcr0.99_m.bed",
    "hela_CTCF_motifs_fcr099.bed"
]
#"#1fa774","#0014A8","#0485d1","#fd3c06"
condition_order = ["hek_J", "hela_H7J", "hek_ctct", "hela_ctcf", "hek_rad21", "hela_rad21","hek_H7J", "heK_ZNF"]
DEFAULT_COLOR_LIST = ["#1fa774", "#1fa774", "#0485d1", "#0485d1", "#fd3c06", "#fd3c06","#1fa774","#0014A8"]
ymax = [ 0.4, 0.4, 5, 5, 3, 3,0.5,1]
ymin = [0, 0, 0, 0, 0, 0,0,0]
proseq_path = "~/PRO-seq/HCT" # provide path to proseq files
proseq_forward_bigwigs = ["HCT116-PRO_merge_hg38_plus_sort.bw","HCT116-PRO_merge_hg38_plus_sort.bw"]
proseq_reverse_bigwigs = ["HCT116-PRO_merge_hg38_minus_sort.bw","HCT116-PRO_merge_hg38_plus_sort.bw"]



bw_bins=7000
track_height=0.5


region_str = ("17:70469719-72839718")

chrom, positions = region_str.split(":")
start, end = map(int, positions.replace(",", "").split("-"))
region_gr = GenomeRange(f"chr{chrom}", start, end)  

MICROC_PLOTTING_PARAMS = {"style":"triangular",
                          "depth_ratio":"full",
                          "balance":True,
                          "cmap":'fall',
                          "max_value":-6,
                         "min_value":-8}

region_plot_params = {"region": region_gr,
                    "resolution": 5000,
                    "highlight_region_list": [["chr1",154_955_735,154_994_524]], 
                    "microc_file": mcool_files[0],
                    "microc_title": "Micro-C",
                    "chip_list_list": [chip_condition_bigwigs],
                    "forward_proseq_list": proseq_forward_bigwigs,
                    "reverse_proseq_list": proseq_reverse_bigwigs,
                     "condition_order": condition_order,
                     "bed_files":bed_files}

frame = make_region_plot(**region_plot_params)
fig1 = frame.plot(region_str)
fig1
fig1.savefig("hek_hela_sox9_withchip.pdf", format='pdf')



#### to plot hic/Microc bachground in blue
def extract_paired_reads(wt_path, ko_path, region, resolution="5000"):
    clr_wt = cooler.Cooler(f"{wt_path}::resolutions/{resolution}")
    clr_ko = cooler.Cooler(f"{ko_path}::resolutions/{resolution}")

    mat_wt = np.array(clr_wt.matrix(balance=True).fetch(region))
    mat_ko = np.array(clr_ko.matrix(balance=True).fetch(region))

    mask = ~np.isnan(mat_wt) & ~np.isnan(mat_ko)
    return mat_wt, mat_ko, mat_wt[mask], mat_ko[mask]


def plot_contact_matrix_triangular(matrix, title, outname, vmax, cmap):
    matrix = np.where(matrix == 0, np.nan, matrix)  
    masked = np.ma.masked_invalid(matrix)

#     tri_mask = np.tri(*masked.shape, k=0, dtype=bool)
#     masked_tri = np.ma.masked_where(tri_mask, masked)

    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)  
    im = ax.imshow(
        masked,
        cmap=cmap,
        origin="upper",
        aspect="equal",
        vmax=vmax,
        interpolation="none" 
    )
    ax.set_title(title, fontsize=14)
    ax.axis("off")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Balanced contact")
    plt.tight_layout()
    plt.savefig(outname, dpi=300, bbox_inches="tight", transparent=True)  
    plt.close()


def main():
    wt_mcool = 'hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool'
    ko_mcool = 'hek_Jko_merge_all_1.7B_UU_dedup_50bp.mcool'
    resolution = 1000
    region = "chr19:23902824-24143247"

    mat_wt, mat_ko, wt_reads, ko_reads = extract_paired_reads(wt_mcool, ko_mcool, region, resolution)

    all_vals = np.concatenate([mat_wt.flatten(), mat_ko.flatten()])
    nonzero_vals = all_vals[(~np.isnan(all_vals)) & (all_vals > 0)]
#     vmax = np.percentile(nonzero_vals, 95)
    vmax=0.003

    # use fall colormap，set NaN to light blue
    custom_cmap = plt.get_cmap('fall').copy()
    custom_cmap.set_bad((0.78, 0.84, 0.93))   ##RGB


    # save PDF
    plot_contact_matrix_triangular(mat_wt, "WT contact matrix", "wt_contact_matrix.pdf", vmax, custom_cmap)
    plot_contact_matrix_triangular(mat_ko, "KO contact matrix", "ko_contact_matrix.pdf", vmax, custom_cmap)

    # diff map
    diff = mat_ko - mat_wt
    vmax_diff = np.nanpercentile(np.abs(diff), 98)
    masked_diff = np.ma.masked_invalid(diff)
    cmap_diff = plt.cm.seismic

    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    im = ax.imshow(masked_diff, cmap=cmap_diff, origin="upper", aspect="equal",
                   vmin=-vmax_diff, vmax=vmax_diff, interpolation='none')
    ax.set_title("KO - WT difference", fontsize=14)
    ax.axis("off")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    plt.savefig("diff_contact_matrix.pdf", dpi=300, bbox_inches="tight", transparent=True)
    plt.close()


if __name__ == "__main__":
    main()
