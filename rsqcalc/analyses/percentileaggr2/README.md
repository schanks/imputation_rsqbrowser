# Aggregate (mean) r2 and proportion of well-imputed variants (Rsq>0.8)
Commands and scripts used to calculate the proportion of variants with Rsq>0.8 and the mean Rsq by MAF bins used in Figures 1,S2,S12. These aggregated files are provided.

## Input files:
../../[ancestry].[variant type].tab variant-level Rsq results files that can be downloaded from https://imputationserver.sph.umich.edu/rsq-browser/downloads

## Proportion of variants with rsq>0.8
Biallelic SNVs:
```
for j in inpsyght biome mlof metsim; do echo "Rscript percentile_bisnv.R -r ../../$j.bi.snv.tab -o $j.bi.snv.perc8.tab"; done > prop.sh
```

All variants (TOPMed only):
```
for j in inpsyght biome mlof metsim; do for k in bi.snv bi.indel multi.snv multi.indel; do echo "Rscript percentile_allvar.R -r ../../$j.$k.tab -o $j.$k.top.perc8.tab"; done; done > prop.sh
```

## Mean Rsq:
```
for j in inpsyght biome mlof metsim; do echo "Rscript mean_bisnv.R -r ../../$j.bi.snv.tab -o $j.bi.snv.mean.tab"; done > mean.sh
```








