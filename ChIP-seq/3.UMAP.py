import pathlib
import pandas as pd
import numpy as np
import pyBigWig
from tqdm import tqdm
import ray
import multiprocessing
import anndata
import scanpy as sc

##read merged peaks
gene_bed = pd.read_csv(
    "hct_real_merged_9factor_sorted_noblack.bed",
    sep="\t",
    header=None,
    usecols=[0, 1, 2],  # 只读前三列
    names=["chrom", "start", "end"],
)

gene_bed.index = (
    gene_bed['chrom'].astype(str) + ':' +
    gene_bed['start'].astype(str) + '-' +
    gene_bed['end'].astype(str)
)

use_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX']
gene_bed = gene_bed[gene_bed['chrom'].isin(use_chroms)].copy()
gene_bed.head()

# Define the path to your chip-seq directory
directory_path = pathlib.Path('/broad/.../chip_Umap/bw_files/HCT/fastq/2.track_peak')

# Get all .bw files
bw_paths = list(directory_path.glob('*.bw'))

# Define keywords to filter for the desired paths
keywords = [
    'ChIP-55',
    'ChIP-56',
    'hct116_h3k27ac_encode2012',
    'HCT116_J_chip_85_87',
    'hct_atac_ENCFF121EPT',
    'hct_wt_merge_all_1.48B_UU_dedup.pairs_50bp.mcool_5000',
    'chip260double',
    'hct_h3k4me3',
    'hct_polIIA'
]
# Filter bw_paths to only include files containing the keywords
filtered_bw_paths = [path for path in bw_paths if any(keyword in str(path) for keyword in keywords)]

# Print the filtered paths
print(filtered_bw_paths)

# Define function to extract BigWig statistics
def get_bw_stats(args):
    bw_path, region_bed = args
    with pyBigWig.open(str(bw_path)) as bw:
        vector = []
        for region_id, (chrom, start, end) in region_bed.iterrows():
            try:
                mean = bw.stats(chrom, start, end)
            except Exception as e:
                mean = 0
                print(f"Error at {chrom}:{start}-{end} in {bw_path}, region_id: {region_id}")
                raise e
            vector.append(mean[0])
    return bw_path.name, vector

# Use multiprocessing to parallelize BigWig processing
num_cpus = 8

with multiprocessing.Pool(processes=num_cpus) as pool:
    results = list(tqdm(pool.imap(get_bw_stats, [(bw, gene_bed) for bw in filtered_bw_paths]), total=len(filtered_bw_paths)))

# Process results
bw_names, matrix = zip(*results)

# Convert to DataFrame
df = pd.DataFrame(matrix, columns=gene_bed.index, index=bw_names).T

### save the file for future use
df.to_csv('/broad/.../chip_Umap/hct_peak_withznf_by_bw_pub_ins.csv.gz')
df = pd.read_csv('/broad/.../chip_Umap/hct_peak_withznf_nipbl_by_bw_pub_ins.csv.gz', index_col=0)

#normalize data
import seaborn as sns
from scipy.stats import zscore
df = df.fillna(0)
df_umap = df
insulation_score = df['hct_wt_merge_all_1.48B_UU_dedup.pairs_50bp.mcool_5000_ins.bedgraph.125000.bw']
df_umap = df_umap.drop(columns=['hct_wt_merge_all_1.48B_UU_dedup.pairs_50bp.mcool_5000_ins.bedgraph.125000.bw'])
### zscore normalization
df_norm0 = zscore(np.log1p(df_umap), axis=0)

### Umap calculation
adata = anndata.AnnData(X=df_norm0.values, 
                obs=pd.DataFrame([], index=df_norm0.index), 
                var=pd.DataFrame([], index=df_norm0.columns))

sc.pp.scale(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
print(adata.obsm.keys())
sc.tl.leiden(adata, resolution=1)


### plot
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
from matplotlib.colors import TwoSlopeNorm
mpl.rcParams['pdf.fonttype']=42 ### save txt to changable font

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Plot UMAP and other enrichment files
fig, axes = plt.subplots(figsize=(12, 7), ncols=4, nrows=3, dpi=300)
adata.obs['x'] = adata.obsm['X_umap'][:, 0]
adata.obs['y'] = adata.obsm['X_umap'][:, 1]

axes = axes.ravel()

# First plot: Clusters with leiden
ax = axes[0]
scatter = sns.scatterplot(ax=ax, data=adata.obs, x='x', y='y', hue='leiden', palette='tab10', s=1.5, linewidth=0, legend=None)
clusters = adata.obs['leiden'].unique()
for cluster in clusters:
    cluster_data = adata.obs[adata.obs['leiden'] == cluster]
    x_mean = cluster_data['x'].mean()
    y_mean = cluster_data['y'].mean()
    ax.text(x_mean, y_mean, cluster, fontsize=12, color='black', ha='center', va='center')
# Rasterize the scatter plot
ax.collections[0].set_rasterized(True)

# define color map
balanced_pink_colors_with_labels = [
    (0.0196078431372549, 0.18823529411764706, 0.3803921568627451, 1.0),  
    (0.12941176470588237, 0.4, 0.6745098039215687, 1.0),  
    (0.2627450980392157, 0.5764705882352941, 0.7647058823529411, 1.0),  
    (0.5725490196078431, 0.7725490196078432, 0.8705882352941177, 1.0),  
    (0.8196078431372549, 0.8980392156862745, 0.9411764705882353, 1.0),  
    (0.9686274509803922, 0.9686274509803922, 0.9686274509803922, 1.0), 
    (0.9921568627450981, 0.8, 0.85, 1.0),  
    (0.9411764705882353, 0.5882352941176471, 0.6666666666666666, 1.0),  
    (0.8627450980392157, 0.39215686274509803, 0.47058823529411764, 1.0),  
    (0.6274509803921569, 0.27450980392156865, 0.3137254901960784, 1.0)  
]



# create colormap
custom_cmap = LinearSegmentedColormap.from_list("custom_cmap",balanced_pink_colors_with_labels)
# Define new colormap with white in the center
new_colormap_colors = [
    "#974da3",  # Purple (Start)
    "#ffffff",  # White (Middle)
    "#91c993"   # Teal (End)
]

# Create custom colormap
custom_cmap1 = LinearSegmentedColormap.from_list("custom_cmap", new_colormap_colors)

# Remaining plots: Enrichment of each file with color scale from -5 to 1
for col, ax in zip(adata.var_names, axes[1:]):
    adata.obs[col] = pd.Series(adata[:, col].X.ravel(), index=adata.obs_names)
    points = ax.scatter(
        adata.obs['x'], adata.obs['y'], c=adata.obs[col],
        s=1.5, linewidth=0, cmap=custom_cmap, vmin=-0.5, vmax=0.5
    )
    ax.set_title(col)
    fig.colorbar(points, ax=ax, orientation='vertical')
    # Rasterize the scatter plot
    points.set_rasterized(True)

# Plot insulation score in the next available slot
ax = axes[len(adata.var_names) + 1]
adata.obs['insulation_score'] = insulation_score
norm = TwoSlopeNorm(vmin=-0.5, vcenter=0, vmax=0.5)
insulation_points = ax.scatter(adata.obs['x'], adata.obs['y'], c=adata.obs['insulation_score'], s=1.5,  linewidth=0, cmap=custom_cmap1,  norm=norm)
ax.collections[0].set_rasterized(True)
fig.colorbar(insulation_points, ax=ax, label='Insulation Score')
ax.set_title('Insulation Score')


# # Plot insulation score in the next available slot
# ax = axes[len(adata.var_names) + 2]
# adata.obs['stripe'] = stripe
# insulation_points = ax.scatter(
#     adata.obs['x'], adata.obs['y'],
#     c=adata.obs['insulation_score'], s=1,
#     linewidth=0, cmap='Reds', vmin=0, vmax=1
# )

ax.collections[0].set_rasterized(True)
# fig.colorbar(insulation_points, ax=ax, label='stripe')
# ax.set_title('stripe')


plt.tight_layout()

plt.show()

# extract point position
output_df = adata.obs[['leiden', 'x', 'y']].copy()

# add rowname
output_df['cell'] = output_df.index

# order
output_df = output_df[['cell', 'leiden', 'x', 'y']]

# save for future use
output_df.to_csv("/broad/boxialab/shawn/projects/analysis/chip_Umap/hct_adata_cell_cluster_umap.bed", sep="\t", index=False, header=False)







