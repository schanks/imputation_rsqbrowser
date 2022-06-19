# Effect of variant caller on imputation quality
Commands and scripts used to compare imputed data from the GotCloud and GATK callsets for the METSIM (Finnish) dataset used for Figure S4.

## Input file: 
bi.snv.tab: Merged R2 file for GATK results.  Analogous to metsim.bi.snv.tab downloadable from RsqBrowser

## Calculate proportion of variants with r2>0.8
```
Rscript ../percentileaggr2/percentile_bisnv.R -r bi.snv.tab -o metsim.gatk.bi.snv.perc8.tab
```
