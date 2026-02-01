**Spike-in downsampling calculation**

For each pairwise comparison (e.g. sample A vs sample B), Drosophila dm6–aligned reads were used as an internal reference.  
The sample with the lower dm6 spike-in count was used as the normalization baseline, and human (hg38) reads in the other sample were downsampled accordingly.

The number of downsampled human reads was calculated as:

\[
\mathrm{hg38}_{\mathrm{downsampled}} =
\frac{\min(\mathrm{dm6}_A,\ \mathrm{dm6}_B)}{\mathrm{dm6}_{\mathrm{sample}}}
\times \mathrm{hg38}_{\mathrm{sample}}
\]

where  
- \(\mathrm{dm6}_A\), \(\mathrm{dm6}_B\) are the numbers of dm6-aligned reads in samples A and B, respectively  
- \(\mathrm{hg38}_{\mathrm{sample}}\) is the number of hg38-aligned reads in the corresponding sample  

The resulting scaling factor was calculated externally (e.g. in Excel) and manually provided to `5.qchip_downsample.py`.
