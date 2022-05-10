
# Microarray Simuluation

This file details how WGS variants were subset to those on each of the 4 arrays. 


## Input files
chr*.wgs.vcf.gz = VCFs of sequenced variants that pass QC (Taliun et al 2021)

## Download strand files (we downloaded on March 23 2020):
**InfiniumCore**
```
https://www.well.ox.ac.uk/~wrayner/strand/InfiniumCore-24v1-1_A1-b38-strand.zip
```
**InfiniumOmniExpress**
```
wget https://www.well.ox.ac.uk/~wrayner/strand/InfiniumOmniExpress-24v1-2_A1-b38-strand.zip
```

**InfiniumOmni2.5M**
```
wget https://www.well.ox.ac.uk/~wrayner/strand/InfiniumOmni2-5-8v1-3_A1-b38-strand.zip
```
**MEGA**
```
wget https://www.well.ox.ac.uk/~wrayner/strand/Multi-EthnicGlobal_A1-b38-strand.zip
```

## Convert strand files to sorted bed files:
```
cat InfiniumCore-24v1-1_A1-b38.strand | awk -v OFS='\t' '{print "chr"$2,$3-1,$3}' > core.unsort.bed
bedtools sort -i core.unsort.bed > core.bed
```
*repeat for other 3 arrays*

## Subset WGS arrays to variants present on each array:
```
for i in {1..22}; do echo "tabix -h chr$i.wgs.vcf.gz -R core.bed | bgzip -c > chr$i.core.unsort.vcf.gz"; done > core.sh
for i in {1..22}; do echo "bcftools sort -Oz chr$i.core.unsort.vcf.gz > chr$i.core.sort.vcf.gz"; done > sort.sh
for i in {1..22}; do echo "tabix -p vcf chr$i.core.sort.vcf.gz"; done > tabix.sh
```
*repeat for other 3 arrays*

## Upload
Sorted arrays (chr*.[array].sort.vcf.gz) are uploaded to [Michigan Imputation Server](https://imputationserver.sph.umich.edu/).


