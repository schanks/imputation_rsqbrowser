# Downsample

Commands and scripts used to downsample individuals for Figure S1. All commands use METSIM dataset as an example. The proportion of variants with Rsq>0.8 and the mean Rsq by MAF bins for each subset are provided.

## Input files:
[ancestry].unrelated.txt: List of unrelated IDs

chr20.[sample].bi.snv.sites: Site lists of variants present in random subsample

all_results.tab: Variant-level Rsq for all subsamples


## Calculate mean Rsq:
METSIM example:
```
Rscript down_mean.R --results all_results.tab --one chr20.100.bi.snv.sites --two chr20.200.bi.snv.sites --three chr20.300.bi.snv.sites --four chr20.400.bi.snv.sites --five chr20.500.bi.snv.sites --six chr20.1000.bi.snv.sites --seven chr20.2000.bi.snv.sites
```
***Repeat for other ancestries to generate [ancestry].all_mean.tab***

## Calculate the proportion of variants with rsq>0.8
METSIM example:
```
Rscript down_perc8.R --results all_results.tab --one chr20.100.bi.snv.sites --two chr20.200.bi.snv.sites --three chr20.300.bi.snv.sites --four chr20.400.bi.snv.sites --five chr20.500.bi.snv.sites --six chr20.1000.bi.snv.sites --seven chr20.2000.bi.snv.sites
```

***Repeat for other ancestries to generate [ancestry].all_perc8.tab***





