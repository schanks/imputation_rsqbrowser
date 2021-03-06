# Genotype concordance
This README includes commands used to calculate the concordance between sequenced and imputed 'best-guess' genotypes. 
The individual-level results files cannot be provided due to privacy concerns.

## Input files:
chr$i.dose.vcf.gz: Imputed VCF downloaded from server (TOPMed) and lifted to b38 (1000G + HRC) 

chr$i.common.vcf.gz: WGS VCF of all common (MAF>5%) biallelic SNVs

chr$i.lowfreq.vcf.gz: WGS VCF of all low frequencey (0.5%<MAF<5%) biallelic SNVs

chr$i.rare.vcf.gz: WGS VCF of all rare (MAF<0.5%) biallelic SNVs

## Download bed-diff.pl
https://raw.githubusercontent.com/statgen/gotcloud/master/scripts/bed-diff.pl.

## Calculate concordance by variant MAF category
```
for i in {1..22}; do for j in common lowfreq rare; do echo "perl bed-diff.pl -vcf1 chr$i.$j.vcf.gz -vcf2 chr$i.dose.vcf.gz -minmaj -o chr$i.$j"; done; done > concordance.sh
```

## Combine over all chromosomes for each individual by variant MAF cateogry
```
for j in common lowfreq rare; do echo "cat chr1.$j.ind > all.$j.ind"; done > merge.sh
for i in {2..22}; do for j in common lowfreq rare; do echo "cat chr$i.$j.ind >> all.$j.ind"; done; done > merge.sh
awk '{ b[$1]; for(i=2;i<=NF;i++)a[$1,i]+=$i }END{for( i in b){printf("%s",i);for(j=2;j<=NF;j++){printf("%s%s",OFS,a[i,j])} print ""}}' all.rare.ind | sed 's/\s/\t/g' > tot.rare.ind
awk '{ b[$1]; for(i=2;i<=NF;i++)a[$1,i]+=$i }END{for( i in b){printf("%s",i);for(j=2;j<=NF;j++){printf("%s%s",OFS,a[i,j])} print ""}}' all.lowfreq.ind | sed 's/\s/\t/g' > tot.lowfreq.ind
awk '{ b[$1]; for(i=2;i<=NF;i++)a[$1,i]+=$i }END{for( i in b){printf("%s",i);for(j=2;j<=NF;j++){printf("%s%s",OFS,a[i,j])} print ""}}' all.common.ind | sed 's/\s/\t/g' > tot.common.ind
```




