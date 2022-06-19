## Effect of real vs. WGS-based arrays
Commands and scripts used to compare imputed data from real Illumina OmniExpress array to WGS-based OmniExpress array for the METSIM (Finnish) dataset used for Figure S3. 

## Input files:
humanomni.bi.snv.tab: Variant-level Rsq for every variant in METSIM WGS. Not provided due to size.

## Proportion of variants with rsq>0.8:
```
Rscript percentile_metsim.R -r humanomni.bi.snv.tab -o metsim.humanomni.bi.snv.perc8.tab 
```

##  Mean rsq:
```
Rscript mean_metsim.R -r humanomni.bi.snv.tab -o metsim.humanomni.bi.snv.mean.tab
```
