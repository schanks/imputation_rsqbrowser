#########################################################################################
This README includes the commands used to calculate Rsq between sequenced genotypes and imputed dosages for each combination of genotyping array and reference panel used in the paper and to merge all results for each variant type. 

Merged Rsq results files can be downloaded at https://imputationserver.sph.umich.edu/rsq-browser/downloads. These are the same files that were generated in Step 2 below.

Scripts used to generate the plots for the paper figures are found in 'figures' directory. All other scripts used to perform analyses reported in the paper are found in 'analyses' directory.
#########################################################################################
Input files:
Please note that all 1000G and HRC files were lifted to b38 before Rsq calculation.
chr$i.wgs.bi.snv.vcf.gz: VCF of sequenced biallelic SNVs
chr$i.wgs.bi.indel.vcf.gz: VCF of sequenced biallelic indels
chr$i.wgs.multi.snv.vcf.gz: VCF of sequenced multiallelic SNVs
chr$i.wgs.multi.indel.vcf.gz: VCF of sequenced multiallelic indels
chr$i.dose.vcf.gz: Imputed VCF downloaded from server (TOPMed) and lifted to b38 (1000G + HRC) 
wgs.[variant type].sites: Tab delimited file with CHR, POS, REF, ALT, AF for each sequenced variant

Scripts:
Rsq.pl: Calculates Rsq between sequenced genotypes and imputed dosages for biallelic variants.
Rsq_multi.pl: Calculates Rsq between sequenced genotypes and imputed dosages for multiallelic variants.
rsq_merge_all.R: Merges Rsq results for all combinations of genotyping arrays and reference panels for biallelic SNVs.
rsq_merge_TOPMed.R: Merges Rsq results for all genotyping arrays for TOPMed reference panel (other variant types).

1. Calculate observed imputation r2 for each combination of genotyping array and reference panel:
Biallelic SNVs:
for i in {1..22}; do echo "perl Rsq.pl --dosage-file chr$i.dose.vcf.gz --geno-file chr$i.wgs.bi.snv.vcf.gz > chr$i.bi.snv"; done > bisnv.sh

Biallelic indels (TOPMed only):
for i in {1..22}; do echo "perl Rsq.pl --dosage-file chr$i.dose.vcf.gz --geno-file chr$i.wgs.bi.indel.vcf.gz > chr$i.bi.indel"; done > biindel.sh

Multiallelic variants (TOPMed only):
for i in {1..22}; do echo "perl Rsq_multi.pl --dosage-file chr$i.dose.vcf.gz --geno-file chr$i.wgs.multi.snv.vcf.gz > chr$i.multi.snv"; done > multisnv.sh

for i in {1..22}; do echo "perl Rsq_multi.pl --dosage-file chr$i.dose.vcf.gz --geno-file chr$i.multi.indel.vcf.gz > chr$i.multi.indel"; done > multiindel.sh

Merge all chromosomes together:
for j in bi.snv bi.indel multi.snv multi.indel; do echo "cat chr1.$j | sed 's/^/chr/g' | sed 's/chrCHR/CHR/g' > all.$j"; done > merge.sh
for j in bi.snv bi.indel multi.snv multi.indel; do for i in {2..22}; do echo "cat chr$i.$j | tail -n +2 | sed 's/^/chr/g' >> all.$j"; done; done > merge.sh

2. Merge all combinations of genotyping arrays and reference panels. Assumes that final output file from Step 1 is of the form all.[array].[panel].[variant type]
Biallelic SNVs:
Rscript rsq_merge_all.R -n wgs.bi.snv.sites -a all.Core.1000G.bi.snv -b all.Core.HRC.bi.snv -c all.Core.TOPMed.bi.snv -e all.OmniExp.1000G.bi.snv -f all.OmniExp.HRC.bi.snv -g all.OmniExp.TOPMed.bi.snv -j all.Omni25.1000G.bi.snv -k all.Omni25.HRC.bi.snv -l all.Omni25.TOPMed.bi.snv -m all.MEGA.1000G.bi.snv -o all.MEGA.HRC.bi.snv -p all.MEGA.TOPMed.bi.snv -z bi.snv

Other variants:
for j in bi.indel multi.snv multi.indel; do echo "Rscript rsq_merge_TOPMed.R -n all.$j.sites -c all.Core.TOPMed.$j -g all.OmniExp.TOPMed.$j -l all.Omni25.TOPMed.$j -p all.MEGA.TOPMed.$j -z $j"; done > merge.sh





