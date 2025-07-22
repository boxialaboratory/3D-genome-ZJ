import cooler
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.stats import wilcoxon, fisher_exact
from matplotlib.ticker import EngFormatter
import itertools

plt.rcParams['font.size'] = 12
bp_formatter = EngFormatter('b')


# ---------- 1. plot function ----------
def pcolormesh_45deg(ax, matrix_c, start=0, resolution=1, *args, **kwargs):
    start_pos_vector = [start + resolution * i for i in range(len(matrix_c) + 1)]
    n = matrix_c.shape[0]
    t = np.array([[1, 0.5], [-1, 0.5]])
    matrix_a = np.dot(
        np.array([(i[1], i[0]) for i in itertools.product(start_pos_vector[::-1], start_pos_vector)]),
        t
    )
    x = matrix_a[:, 1].reshape(n + 1, n + 1)
    y = matrix_a[:, 0].reshape(n + 1, n + 1)
    im = ax.pcolormesh(x, y, np.flipud(matrix_c), *args, **kwargs)
    im.set_rasterized(True)
    return im


def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x', rotation=45)


# ---------- 2. extract balanced contact matrix ----------
def extract_paired_reads(wt_path, ko_path, region_tuple, resolution):
    clr_wt = cooler.Cooler(f"{wt_path}::resolutions/{resolution}")
    clr_ko = cooler.Cooler(f"{ko_path}::resolutions/{resolution}")
    mat_wt = np.array(clr_wt.matrix(balance=True).fetch(region_tuple))
    mat_ko = np.array(clr_ko.matrix(balance=True).fetch(region_tuple))
    mask = ~np.isnan(mat_wt) & ~np.isnan(mat_ko)
    return mat_wt, mat_ko, mat_wt[mask], mat_ko[mask]


# ---------- 3. plot diff ----------
def plot_contact_matrix_trio(mat_wt, mat_ko, region_start, resolution, vmax, diff_vmax, cmap, window):
    fig, axes = plt.subplots(nrows=3, figsize=(18, 9), sharex=True)

    diff = mat_ko - mat_wt
    matrices = [mat_wt, mat_ko, diff]
    titles = ["WT contact matrix", "KO contact matrix", "KO - WT difference"]
    cmaps = [cmap, cmap, plt.cm.seismic]
    vmins = [0, 0, -diff_vmax]
    vmaxs = [vmax, vmax, diff_vmax]
    norms = [
        Normalize(vmin=0, vmax=vmax),
        Normalize(vmin=0, vmax=vmax),
        Normalize(vmin=-diff_vmax, vmax=diff_vmax)
    ]

    for ax, mat, title, cm, norm in zip(axes, matrices, titles, cmaps, norms):
        if "difference" not in title:
            mat = np.where(mat == 0, np.nan, mat)  # set 0 to blue

        im = pcolormesh_45deg(
            ax, mat,
            start=region_start,
            resolution=resolution,
            cmap=cm,
            norm=norm
        )
        ax.set_title(title, fontsize=12)
        ax.set_aspect(0.5)
        ax.set_ylim(0, 40 * window)
        ax.xaxis.set_visible(False)
        format_ticks(ax, rotate=False)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="1%", pad=0.05)
        label = "Balanced contact" if "difference" not in title else "KO - WT diff"
        plt.colorbar(im, cax=cax, label=label)

    axes[-1].xaxis.set_visible(True)
    plt.tight_layout()
    plt.savefig("hek_znf_ko_example-nolog.pdf")
    plt.show()


# ---------- 4. main program ----------
def main():
    wt_mcool = "WT.mcool"
    ko_mcool = "KO.mcool"
    resolution = 10000
    region = "chr12:14059779-17038096"

    # Parse
    chrom, coords = region.split(":")
    start, end = map(int, coords.split("-"))
    region_tuple = (chrom, start, end)
    region_start = start

    bin_count = 90
    window = (end - start) // bin_count

    mat_wt, mat_ko, wt_reads, ko_reads = extract_paired_reads(wt_mcool, ko_mcool, region_tuple, resolution)

    # vmax（95 percentile）
    all_vals = np.concatenate([mat_wt.flatten(), mat_ko.flatten()])
    nonzero_vals = all_vals[(~np.isnan(all_vals)) & (all_vals > 0)]
    vmax = np.percentile(nonzero_vals, 95)

    # ✅ use fall colormap，set 0/NaN to blue
    custom_cmap = plt.get_cmap('fall').copy()
    custom_cmap.set_bad('blue')

    
    diff = mat_ko - mat_wt
    vmax_diff = np.nanpercentile(np.abs(diff), 99)

    
    plot_contact_matrix_trio(mat_wt, mat_ko, region_start, resolution, vmax, vmax_diff, custom_cmap, window)

    # statistical test
    pval = wilcoxon(wt_reads, ko_reads).pvalue
    print(f"Wilcoxon P = {pval:.3e}")
    a = np.sum((mat_wt > 0) & ~np.isnan(mat_wt))
    b = np.sum((mat_wt == 0) & ~np.isnan(mat_wt))
    c = np.sum((mat_ko > 0) & ~np.isnan(mat_ko))
    d = np.sum((mat_ko == 0) & ~np.isnan(mat_ko))
    table = np.array([[a, b], [c, d]])
    _, fisher_p = fisher_exact(table)
    print(f"Fisher P = {fisher_p:.3e}, Odds Ratio = {(a*d)/(b*c + 1e-10):.2f}")


if __name__ == "__main__":
    main()
