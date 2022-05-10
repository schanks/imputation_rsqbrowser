# Imputation quality by predicted variant consequences

Commands used to calculate proportion of well-imputed biallelic SNVs by predicted functional impact (VEP) for Figure S11. Aggregate results files are provided.

## Download
```
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
```

## Input files:
chr*.bi.snv.vcf.gz: WGS VCF of sequenced biallelic SNVs

../../[ancestry].bi.snv.tab: variant-level Rsq results files that can be downloaded from [RsqBrowser](https://imputationserver.sph.umich.edu/rsq-browser/downloads)

## Annotate sequenced variants:
```
for i in {1..22}; do echo "ensembl-vep/vep -i ../chr$i.bi.snv.vcf.gz --cache -o chr$i.vep.txt"; done > command.sh
```

## Extract relevant columns, aggregate chromosomes, remove duplicate rows:
```
for i in {1..22}; do echo "cat chr$i.vep.txt | grep -v \"##\" | cut -f2,3,14 | sed 's/:/\t/g' | sed 's/;/\t/g' | cut -f1,2,3,4 | tail -n +2 > chr$i.impact"; done > command.sh
cat chr1.impact > all.impact
for i in {2..22}; do echo "cat chr$i.impact >> all.impact"; done > command.sh
cat all.impact | sort | uniq -c > all.vep
```

## Calculate the proportion of variants with rsq>0.8
```
Rscript percentile_vep_bisnv.R -r ../../[ancestry].bi.snv.tab -v [ancestry].all.vep -o [ancestry].vep.perc8.tab
```

