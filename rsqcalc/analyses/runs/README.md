# Runs of consecutively well-imputed variants
Script (Switch.R) and commands used to calculate length of runs of consecutively well-imputed variants for Figure 3. Due to size, the output files cannot be provided.

## Input files:
chr$i.bi.snv.sites: Tab delimited file with CHR, POS, REF, ALT, AF for each sequenced biallelic SNV

chr$i.bi.snv (output of Rsq.pl, data can be constructed from RsqBrowser downloads)

centromere: downloaded from UCSC Table Browser (Mapping and Sequencing - Centromeres)

## Calculate length of runs of consecutively well-imputed variants
Common variants (MAF>5%)
```
for i in {1..22}; do echo "Rscript Switch.R -w chr$i.bi.snv.sites -i chr$i.bi.snv -m 0.5 -n 0.05 -t 0.8 -c centromere -o chr$i.common.stretch
```

Low frequency variants (0.5%<MAF<5%)
```
for i in {1..22}; do echo "Rscript Switch.R -w chr$i.bi.snv.sites -i chr$i.bi.snv -m 0.05 -n 0.005 -t 0.8 -c centromere -o chr$i.lowfreq.stretch
```

Rare variants (MAF<5%)
```
for i in {1..22}; do echo "Rscript Switch.R -w chr$i.bi.snv.sites -i chr$i.bi.snv -m 0.005 -n 0.0002 -t 0.8 -c centromere -o chr$i.rare.stretch
```

## Aggregate over all chromosomes:
```
for j in common lowfreq rare; do echo "cat chr1.$j.stretch > all.$j.stretch"; done > command.sh
for j in common lowfreq rare; do for i in {2..22}; do echo "cat chr$i.$j.stretch >> all.$j.stretch"; done; done > command.sh
```


