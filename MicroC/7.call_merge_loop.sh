## call loop at different resolution
python3 ~/mustache.py -f AA.mcool -r 10kb -pt 0.1 -st 0.88 -o hic_10k_res -p 8 --octaves 2 

### merge loops at different resolutions
cat hic_10k_res.txt hic_5k_res.txt hic_2k_res.txt hic_1k_res.txt >> merge_loop_10k_5k_2k_1k.bed
less merge_loop_10k_5k_2k_1k.bed|awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,"loop_"NR}' > hic_loop_all.txt
bedtools sort -i hic_loop_all.txt > hic_loop_all_sorted.txt
