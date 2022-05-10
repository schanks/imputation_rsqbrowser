# Analyses
This directory contains all scripts used to generate intermediate results files from merged Rsq results files downloaded at https://imputationserver.sph.umich.edu/rsq-browser/downloads. The intermediate results files are also provided as size and privacy concerns allow.

## Table of contents
1. [concordance](https://github.com/schanks/imputation_rsqbrowser/tree/main/rsqcalc/analyses/concordance): Commands used to calculate concordance between sequenced and imputed 'best-guess' genotypes and to generate intermediate results for Figures 2,S5,S6,S7. Individual-level summary statistics cannot be provided due to privacy concerns. 

2. [downsample](https://github.com/schanks/imputation_rsqbrowser/tree/main/rsqcalc/analyses/downsample): Commands and scripts used to downsample individuals for Figure S1. The proportion of variants with Rsq>0.8 and the mean Rsq by MAF bins for each subset are provided. 

3. [genomicfeatures](https://github.com/schanks/imputation_rsqbrowser/tree/main/rsqcalc/analyses/genomicfeatures):Data sources and scripts used to perform analyses associating Rsq with local genomic features for Figures 4,S9,S10. The output from statistical analyses are provided.

4. [omniexp](https://github.com/schanks/imputation_rsqbrowser/tree/main/rsqcalc/analyses/omniexp): Commands and scripts used to compare imputed data from real Illumina OmniExpress array to WGS-based OmniExpress array for the METSIM (Finnish) dataset used for Figure S3. The proportion of variants with Rsq>0.8 for both array types are provided.

5. [percentileaggr2](https://github.com/schanks/imputation_rsqbrowser/tree/main/rsqcalc/analyses/percentileaggr2): Commands and scripts used to calculate the proportion of variants with Rsq>0.8 and the mean Rsq by MAF bins used in Figures 1,S2,S12. These aggregated files are provided.

6. [runs](https://github.com/schanks/imputation_rsqbrowser/tree/main/rsqcalc/analyses/runs): Scripts used to calculate length of runs of consecutively well-imputed variants for Figure 3. The results files are not provided due to size but can be reproduced if desired.

7. [varcaller](https://github.com/schanks/imputation_rsqbrowser/tree/main/rsqcalc/analyses/varcaller): Commands and scripts used to compare imputed data from the GotCloud and GATK callsets for the METSIM (Finnish) dataset used for Figure S4. The proportion of variants with Rsq>0.8 for both variant callers are provided.

8. [VEP](https://github.com/schanks/imputation_rsqbrowser/tree/main/rsqcalc/analyses/VEP): Commands used to calculate proportion of well-imputed biallelic SNVs by predicted functional impact (VEP) for Figure S11. The proportion of variants with Rsq>0.8 by VEP category are provided.



