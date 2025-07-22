## first step: get peak loops within a range
import pandas as pd
from multiprocessing import Pool
import argparse

def parse_bed(file_path):
    """
    Parse BED file into a DataFrame.
    """
    df = pd.read_csv(file_path, sep="\t", header=None, usecols=[0, 1, 2])
    df.columns = ["chrom", "start", "end"]
    return df

def filter_interactions(df_a, df_b, mindist, maxdist):
    """
    Generate all pairwise interactions between two BED files,
    filtered by the distance between the centers of the intervals.
    """
    merged = pd.merge(df_a, df_b, on="chrom", suffixes=("_a", "_b"))
    merged["dist"] = abs((merged["start_b"] + merged["end_b"]) / 2 -
                         (merged["start_a"] + merged["end_a"]) / 2)
    filtered = merged[(merged["dist"] >= mindist) & (merged["dist"] <= maxdist)]
    return filtered[["chrom", "start_a", "end_a", "chrom","start_b", "end_b", "dist"]]

def process_chrom(chrom, df_a, df_b, mindist, maxdist):
    """
    Process a single chromosome to generate interactions.
    """
    df_a_chrom = df_a[df_a["chrom"] == chrom]
    df_b_chrom = df_b[df_b["chrom"] == chrom]
    return filter_interactions(df_a_chrom, df_b_chrom, mindist, maxdist)

def generate_interactions(bed_a, bed_b, mindist, maxdist, output_file, nproc):
    """
    Generate interactions between two BED files, saving results to output_file.
    """
    df_a = parse_bed(bed_a)
    df_b = parse_bed(bed_b)
    
    common_chroms = set(df_a["chrom"]).intersection(df_b["chrom"])
    
    with Pool(nproc) as pool:
        results = pool.starmap(process_chrom, [(chrom, df_a, df_b, mindist, maxdist) for chrom in common_chroms])
    
    result_df = pd.concat(results)
    result_df.to_csv(output_file, sep="\t", index=False, header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate BEDPE interactions from two BED files.")
    parser.add_argument("--bed_a", required=True, help="Path to first BED file")
    parser.add_argument("--bed_b", required=True, help="Path to second BED file")
    parser.add_argument("--output", required=True, help="Output BEDPE file")
    parser.add_argument("--mindist", type=int, required=True, help="Minimum distance between interactions")
    parser.add_argument("--maxdist", type=int, required=True, help="Maximum distance between interactions")
    parser.add_argument("--nproc", type=int, default=1, help="Number of processes for parallelization")
      
    args = parser.parse_args()
    
    generate_interactions(
        bed_a=args.bed_a,
        bed_b=args.bed_b,
        mindist=args.mindist,
        maxdist=args.maxdist,
        output_file=args.output,
        nproc=args.nproc
    )

## second step: generate a .sh file to compute APA and plotting
import os

# List of input BED files located in ../
files = [
    "../hek_J_on_C_no_ZJR.bed",
    "../hek_J_on_h3k27ac_J_noCRZ.bed",
    "../hek_ctcf_on_CRZJ_descend.bed",
    "../hek_ctcf_on_CR_no_JZ_descend.bed",
    "../hek_ctcf_on_CZR_noJ_descend.bed"
]

# Input COOL file path
mcool_path = "../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000"

# Output directory (current working directory)
output_dir = "."

# Create work.sh file to store commands
with open("work.sh", "w") as f_out:
    f_out.write("#!/bin/bash\n\n")  # Add shebang for the script
    
    # Loop through all pairwise combinations
    for i in range(len(files)):
        for j in range(i, len(files)):
            file_a = files[i]
            file_b = files[j]

            # Generate BEDPE file name (without ../ for output)
            bedpe_name = f"{os.path.basename(file_a)[:-4]}_{os.path.basename(file_b)[:-4]}.bedpe"
            bedpe_path = os.path.join(output_dir, bedpe_name)

            # Generate call_bedpe.py command
            call_bedpe_cmd = (
                f"python ../call_bedpe.py "
                f"--bed_a {file_a} --bed_b {file_b} "
                f"--output {bedpe_path} --mindist 100000 --maxdist 2000000 --nproc 48"
            )
            f_out.write(call_bedpe_cmd + "\n")

            # Generate output TXT file name for coolpuppy
            out_name = f"{os.path.basename(file_a)[:-4]}_{os.path.basename(file_b)[:-4]}_coolpup.txt"
            out_path = os.path.join(output_dir, out_name)
            
            # Generate coolpuppy command
            coolpuppy_cmd = (
                f"coolpup.py {mcool_path} {bedpe_path} "
                f"--outname {out_path} "
                f"--flank 28000 --mindist 0 --n_proc 48"
            )
            f_out.write(coolpuppy_cmd + "\n")

#### the generated sh file:
#!/bin/bash

# python ../call_bedpe.py --bed_a ../hek_J_on_C_no_ZJR.bed --bed_b ../hek_J_on_C_no_ZJR.bed --output ./hek_J_on_C_no_ZJR_hek_J_on_C_no_ZJR.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_J_on_C_no_ZJR_hek_J_on_C_no_ZJR.bedpe --outname ./hek_J_on_C_no_ZJR_hek_J_on_C_no_ZJR_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_J_on_C_no_ZJR.bed --bed_b ../hek_J_on_h3k27ac_J_noCRZ.bed --output ./hek_J_on_C_no_ZJR_hek_J_on_h3k27ac_J_noCRZ.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_J_on_C_no_ZJR_hek_J_on_h3k27ac_J_noCRZ.bedpe --outname ./hek_J_on_C_no_ZJR_hek_J_on_h3k27ac_J_noCRZ_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_J_on_C_no_ZJR.bed --bed_b ../hek_ctcf_on_CRZJ_descend.bed --output ./hek_J_on_C_no_ZJR_hek_ctcf_on_CRZJ_descend.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_J_on_C_no_ZJR_hek_ctcf_on_CRZJ_descend.bedpe --outname ./hek_J_on_C_no_ZJR_hek_ctcf_on_CRZJ_descend_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_J_on_C_no_ZJR.bed --bed_b ../hek_ctcf_on_CR_no_JZ_descend.bed --output ./hek_J_on_C_no_ZJR_hek_ctcf_on_CR_no_JZ_descend.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_J_on_C_no_ZJR_hek_ctcf_on_CR_no_JZ_descend.bedpe --outname ./hek_J_on_C_no_ZJR_hek_ctcf_on_CR_no_JZ_descend_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_J_on_C_no_ZJR.bed --bed_b ../hek_ctcf_on_CZR_noJ_descend.bed --output ./hek_J_on_C_no_ZJR_hek_ctcf_on_CZR_noJ_descend.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_J_on_C_no_ZJR_hek_ctcf_on_CZR_noJ_descend.bedpe --outname ./hek_J_on_C_no_ZJR_hek_ctcf_on_CZR_noJ_descend_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_J_on_h3k27ac_J_noCRZ.bed --bed_b ../hek_J_on_h3k27ac_J_noCRZ.bed --output ./hek_J_on_h3k27ac_J_noCRZ_hek_J_on_h3k27ac_J_noCRZ.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_J_on_h3k27ac_J_noCRZ_hek_J_on_h3k27ac_J_noCRZ.bedpe --outname ./hek_J_on_h3k27ac_J_noCRZ_hek_J_on_h3k27ac_J_noCRZ_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_J_on_h3k27ac_J_noCRZ.bed --bed_b ../hek_ctcf_on_CRZJ_descend.bed --output ./hek_J_on_h3k27ac_J_noCRZ_hek_ctcf_on_CRZJ_descend.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_J_on_h3k27ac_J_noCRZ_hek_ctcf_on_CRZJ_descend.bedpe --outname ./hek_J_on_h3k27ac_J_noCRZ_hek_ctcf_on_CRZJ_descend_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_J_on_h3k27ac_J_noCRZ.bed --bed_b ../hek_ctcf_on_CR_no_JZ_descend.bed --output ./hek_J_on_h3k27ac_J_noCRZ_hek_ctcf_on_CR_no_JZ_descend.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_J_on_h3k27ac_J_noCRZ_hek_ctcf_on_CR_no_JZ_descend.bedpe --outname ./hek_J_on_h3k27ac_J_noCRZ_hek_ctcf_on_CR_no_JZ_descend_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_J_on_h3k27ac_J_noCRZ.bed --bed_b ../hek_ctcf_on_CZR_noJ_descend.bed --output ./hek_J_on_h3k27ac_J_noCRZ_hek_ctcf_on_CZR_noJ_descend.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_J_on_h3k27ac_J_noCRZ_hek_ctcf_on_CZR_noJ_descend.bedpe --outname ./hek_J_on_h3k27ac_J_noCRZ_hek_ctcf_on_CZR_noJ_descend_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_ctcf_on_CRZJ_descend.bed --bed_b ../hek_ctcf_on_CRZJ_descend.bed --output ./hek_ctcf_on_CRZJ_descend_hek_ctcf_on_CRZJ_descend.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_ctcf_on_CRZJ_descend_hek_ctcf_on_CRZJ_descend.bedpe --outname ./hek_ctcf_on_CRZJ_descend_hek_ctcf_on_CRZJ_descend_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_ctcf_on_CRZJ_descend.bed --bed_b ../hek_ctcf_on_CR_no_JZ_descend.bed --output ./hek_ctcf_on_CRZJ_descend_hek_ctcf_on_CR_no_JZ_descend.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_ctcf_on_CRZJ_descend_hek_ctcf_on_CR_no_JZ_descend.bedpe --outname ./hek_ctcf_on_CRZJ_descend_hek_ctcf_on_CR_no_JZ_descend_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_ctcf_on_CRZJ_descend.bed --bed_b ../hek_ctcf_on_CZR_noJ_descend.bed --output ./hek_ctcf_on_CRZJ_descend_hek_ctcf_on_CZR_noJ_descend.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_ctcf_on_CRZJ_descend_hek_ctcf_on_CZR_noJ_descend.bedpe --outname ./hek_ctcf_on_CRZJ_descend_hek_ctcf_on_CZR_noJ_descend_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_ctcf_on_CR_no_JZ_descend.bed --bed_b ../hek_ctcf_on_CR_no_JZ_descend.bed --output ./hek_ctcf_on_CR_no_JZ_descend_hek_ctcf_on_CR_no_JZ_descend.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_ctcf_on_CR_no_JZ_descend_hek_ctcf_on_CR_no_JZ_descend.bedpe --outname ./hek_ctcf_on_CR_no_JZ_descend_hek_ctcf_on_CR_no_JZ_descend_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_ctcf_on_CR_no_JZ_descend.bed --bed_b ../hek_ctcf_on_CZR_noJ_descend.bed --output ./hek_ctcf_on_CR_no_JZ_descend_hek_ctcf_on_CZR_noJ_descend.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_ctcf_on_CR_no_JZ_descend_hek_ctcf_on_CZR_noJ_descend.bedpe --outname ./hek_ctcf_on_CR_no_JZ_descend_hek_ctcf_on_CZR_noJ_descend_coolpup.txt --flank 28000 --mindist 0 --n_proc 48
# python ../call_bedpe.py --bed_a ../hek_ctcf_on_CZR_noJ_descend.bed --bed_b ../hek_ctcf_on_CZR_noJ_descend.bed --output ./hek_ctcf_on_CZR_noJ_descend_hek_ctcf_on_CZR_noJ_descend.bedpe --mindist 100000 --maxdist 2000000 --nproc 48
# coolpup.py ../hek_wt_merge_all_1.7B_UU_dedup_50bp.mcool::resolutions/1000 ./hek_ctcf_on_CZR_noJ_descend_hek_ctcf_on_CZR_noJ_descend.bedpe --outname ./hek_ctcf_on_CZR_noJ_descend_hek_ctcf_on_CZR_noJ_descend_coolpup.txt --flank 28000 --mindist 0 --n_proc 48


### draw coolpuppy:

# for txt in *.txt; do 
#     pdf="${txt%.txt}_plot.pdf"
#     plotpup.py --input_pups "$txt" --not_symmetric --vmin 0.9 --vmax 3.3 --output "$pdf"
# done
